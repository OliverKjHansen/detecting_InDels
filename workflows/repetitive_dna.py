import math
from functools import partial
from os.path import expanduser

from gwf import Workflow

from pipeline import Pipeline

# ---------------------------------------------------------------------------
# Repetitive DNA annotation. Annotates tandem repeats per contig for every
# period k, under three models:
#
#   M0 — mreps, exact/pure repeats
#   M1 — TRF, standard scoring from RepeatMasker
#   M2 — mreps, resolution 1: one substitution allowed between adjacent copies
#
# and the excluded group, which neither InDel workflow measures: homopolymer
# runs, and tracts with periods above the models', from TRF.
#
# Configuration is read from config.yaml. Run gwf from the repository root,
# after precompute.py has split the assemblies into contigs:
#
#   gwf -f workflows/precompute.py run
#   gwf -f workflows/repetitive_dna.py run
#
# Paths come from the `layout` templates and the parameters from the `repeats`
# block in config.yaml, through pipeline.py.
# ---------------------------------------------------------------------------

pipeline          = Pipeline()
config            = pipeline.config
assembly_type     = config["assembly_type"]
scripts           = pipeline.scripts
mreps             = expanduser(config["tools"]["mreps"])
trf               = expanduser(config["tools"]["trf"])

# Tandem repeat periods the models annotate. Homopolymers (k=1) and longer
# periods are the excluded group, below.
periods = config["repeats"]["periods"]
periods_str = " ".join(str(k) for k in periods)
max_period = max(periods)

# M0 and M2 — mreps minimum exponent, i.e. at least this many full copies.
mreps_exp = config["repeats"]["mreps_min_exponent"]

# M2 — mreps resolution: substitutions allowed between adjacent copies.
m2_res = config["repeats"]["m2_resolution"]

# M1 — TRF parameters: match mismatch indel PM PI minscore.
trf_params = config["repeats"]["trf_params"]

# M1 — TRF runs on windows of a contig, one job per window, overlapping by
# trf_overlap (see config.yaml).
trf_window  = config["repeats"]["trf_window_size"]
trf_overlap = config["repeats"]["trf_overlap"]

# M1 — TRF aligns a repeat over at most -l million bp (default 2). In a uniform
# array longer than that, e.g. a 4 Mb centromere, it keeps restarting the
# alignment and a window runs for days instead of minutes. So the limit covers
# everything a window job fetches, and no array can exceed it. TRF allocates
# ~0.6 GB per million up front: -l 12 for 10 Mb windows is ~7.5 GB.
trf_max_length = math.ceil((trf_window + 2 * trf_overlap) / 1_000_000)

models = config["repeats"]["models"]

# The excluded group (decisions Q8, Q16): homopolymer runs at least this long,
# and tracts with periods from the first above `periods` up to long_period_max.
homopolymer_min_length = config["repeats"]["homopolymer_min_length"]
long_period_max        = config["repeats"]["long_period_max"]
long_periods           = f"{max_period + 1}-{long_period_max}"

# Resource requests. mreps and TRF are single-threaded, so one core per job.
gwf = Workflow(
    defaults={
        'cores': 1,
        'memory': '10g',
        'walltime': '8:00:00',
        'account': 'mutationalscanning'})


for species, samples in config["samples"].items():
    for sample in samples:
        path        = partial(pipeline.path, species=species, sample=sample)
        assembly    = path("assembly")
        outdir      = pipeline.folder("tandem_repeats", species, sample)
        per_k_dir   = pipeline.folder("repeats_per_k", species, sample)
        bases_dir   = pipeline.folder("repeat_bases", species, sample)
        excluded_dir = pipeline.folder("excluded", species, sample)
        long_period_dir = pipeline.folder("long_periods_window", species, sample)
        resultsdir  = pipeline.folder("repeats_summary", species, sample)
        long_contigs = pipeline.long_contigs(species, sample)

        for contig, length in long_contigs.items():
            contig_fa = path("contig_fasta", contig=contig)

            # Step 1a — one job per model per period k. Each job writes its BED
            # to a .tmp file and renames it only when every command succeeded, so
            # a failed job never leaves a finished-looking output behind.
            for k in periods:

                # M0: mreps exact/pure repeats.
                # -exp N:      at least N full copies (exponent = length/period)
                # -allowsmall: allow tracts shorter than the random-occurrence threshold
                m0_bed = path("repeats_per_k", contig=contig, model="M0", k=k)
                gwf.target(
                    f"M0_mreps_k{k}_{contig}_{species}_{sample}_{assembly_type}",
                    inputs=[contig_fa],
                    outputs=[m0_bed],
                ) << f"""
                set -e
                mkdir -p {per_k_dir}
                WORK=$(mktemp -d)
                trap 'rm -rf "$WORK"' EXIT
                {mreps} -minperiod {k} -maxperiod {k} -exp {mreps_exp} -allowsmall -fasta {contig_fa} > $WORK/mreps.txt
                python {scripts}/parse_mreps_canonical.py $WORK/mreps.txt {contig} {k} {species} {sample} M0 > {m0_bed}.tmp
                mv {m0_bed}.tmp {m0_bed}
                """

                # M2: mreps with resolution m2_res.
                # -res N: allow up to N substitutions between adjacent copies
                m2_bed = path("repeats_per_k", contig=contig, model="M2", k=k)
                gwf.target(
                    f"M2_mreps_k{k}_{contig}_{species}_{sample}_{assembly_type}",
                    inputs=[contig_fa],
                    outputs=[m2_bed],
                ) << f"""
                set -e
                mkdir -p {per_k_dir}
                WORK=$(mktemp -d)
                trap 'rm -rf "$WORK"' EXIT
                {mreps} -minperiod {k} -maxperiod {k} -exp {mreps_exp} -res {m2_res} -allowsmall -fasta {contig_fa} > $WORK/mreps.txt
                python {scripts}/parse_mreps_canonical.py $WORK/mreps.txt {contig} {k} {species} {sample} M2 > {m2_bed}.tmp
                mv {m2_bed}.tmp {m2_bed}
                """

            # Step 1a (M1) — TRF, one job per window of the contig. One run
            # covers every period, so the k-loop below only splits its output.
            #
            # Windows are cut with an overlap on both sides, and a window keeps
            # only the tracts whose start lies in its own core [core_start,
            # core_end), so neighbouring windows never report the same tract.
            # The overlap is what lets a tract crossing a boundary still be
            # found whole by the window that owns it.
            m1_window_beds = []

            for core_start in range(0, length, trf_window):
                core_end    = min(core_start + trf_window, length)
                fetch_start = max(0, core_start - trf_overlap)
                fetch_end   = min(length, core_end + trf_overlap)
                window_beds = [path("repeats_per_k_window", contig=contig, model="M1", k=k, start=core_start)
                               for k in periods]
                m1_window_beds += window_beds
                # The window's BED for the period in the shell loop's $k.
                k_window_bed = path("repeats_per_k_window", contig=contig, model="M1", k="${k}", start=core_start)

                gwf.target(
                    f"M1_trf_{contig}_{core_start}_{species}_{sample}_{assembly_type}",
                    inputs=[assembly, path("assembly_fai")],
                    outputs=window_beds,
                    walltime="8:00:00",
                ) << f"""
                set -e
                mkdir -p {per_k_dir}
                WORK=$(mktemp -d)
                trap 'rm -rf "$WORK"' EXIT
                samtools faidx {assembly} {contig}:{fetch_start + 1}-{fetch_end} > $WORK/window.fa
                # -ngs -> machine-readable stdout, -h -> suppress HTML
                {trf} $WORK/window.fa {trf_params} {max_period} -ngs -h -l {trf_max_length} > $WORK/trf.txt
                for k in {periods_str}; do
                    awk -v k=$k '/^@/{{print}} !/^@/{{if ($3==k) print}}' $WORK/trf.txt > $WORK/trf_${{k}}.txt
                    python {scripts}/parse_trf.py $WORK/trf_${{k}}.txt $k {species} {sample} M1 \
                        {contig} {fetch_start} {core_start} {core_end} \
                        > {k_window_bed}.tmp
                    mv {k_window_bed}.tmp {k_window_bed}
                done
                """

            # Step 1b — per model: join the per-k BEDs of the contig (header kept
            # from the first file only) and resolve overlapping tracts, keeping
            # the smallest k. M1 also has one set of BEDs per window.
            #
            # Resolution drops a tract that overlaps a kept one whole, so the
            # dropped tract's other bases are in no tract of the model, and
            # workflow 4 never sees them. Workflow 3 masks them all the same:
            # its mask is built from every base in a tract before resolution,
            # written here as merged intervals. Those bases are then in neither
            # workflow, like the excluded group.
            k_beds_by_model = {
                "M0": [path("repeats_per_k", contig=contig, model="M0", k=k) for k in periods],
                "M1": m1_window_beds,
                "M2": [path("repeats_per_k", contig=contig, model="M2", k=k) for k in periods],
            }

            for model in models:
                k_beds = k_beds_by_model[model]
                k_bed_list = " ".join(k_beds)
                contig_bed = path("tandem_repeats", contig=contig, model=model)
                bases_bed  = path("repeat_bases", contig=contig, model=model)

                gwf.target(
                    f"{model}_merge_k_{contig}_{species}_{sample}_{assembly_type}",
                    inputs=k_beds,
                    outputs=[contig_bed, bases_bed],
                    memory="16g",
                    walltime="02:00:00",
                ) << f"""
                set -e
                mkdir -p {outdir} {bases_dir}
                WORK=$(mktemp -d)
                trap 'rm -rf "$WORK"' EXIT
                awk 'FNR==1 && NR==1 {{print}} FNR>1 {{print}}' {k_bed_list} > $WORK/all.bed
                python {scripts}/resolve_repeat_overlaps.py $WORK/all.bed > {contig_bed}.tmp
                python {scripts}/merge_regions.py $WORK/all.bed > {bases_bed}.tmp
                mv {contig_bed}.tmp {contig_bed}
                mv {bases_bed}.tmp {bases_bed}
                """

            # Step 3a — tracts with periods above the models', from TRF run as
            # for M1 (same scoring, same windows) but up to long_period_max. M1
            # keeps its own run, so M1 does not change. These tracts are only
            # masked, so a window keeps every tract it reports, not just those
            # starting in its core: an array longer than the overlap is then
            # covered by the windows' pieces, with no hole between them.
            long_period_window_beds = []

            for core_start in range(0, length, trf_window):
                core_end    = min(core_start + trf_window, length)
                fetch_start = max(0, core_start - trf_overlap)
                fetch_end   = min(length, core_end + trf_overlap)
                window_bed  = path("long_periods_window", contig=contig, start=core_start)
                long_period_window_beds.append(window_bed)

                # In satellite arrays TRF takes far longer up to period 2000 than
                # up to 6: some windows ran out of 8 h on the cluster (2026-10-23).
                gwf.target(
                    f"long_period_trf_{contig}_{core_start}_{species}_{sample}_{assembly_type}",
                    inputs=[assembly, path("assembly_fai")],
                    outputs=[window_bed],
                    walltime="48:00:00",
                ) << f"""
                set -e
                mkdir -p {long_period_dir}
                WORK=$(mktemp -d)
                trap 'rm -rf "$WORK"' EXIT
                samtools faidx {assembly} {contig}:{fetch_start + 1}-{fetch_end} > $WORK/window.fa
                {trf} $WORK/window.fa {trf_params} {long_period_max} -ngs -h -l {trf_max_length} > $WORK/trf.txt
                python {scripts}/parse_trf.py $WORK/trf.txt {long_periods} {species} {sample} long_period \\
                    {contig} {fetch_start} > {window_bed}.tmp
                mv {window_bed}.tmp {window_bed}
                """

            # Step 3b — the excluded group of the contig: its homopolymer runs,
            # its long-period tracts joined from the windows of step 3a as
            # merged intervals (neighbouring windows report the tracts in their
            # overlap twice), and the union of the two, which workflow 3 masks
            # and workflow 4 drops the tracts it overlaps.
            homopolymer_bed = path("homopolymers", contig=contig)
            long_period_bed = path("long_periods", contig=contig)
            excluded_bed    = path("excluded", contig=contig)

            gwf.target(
                f"excluded_{contig}_{species}_{sample}_{assembly_type}",
                inputs=[contig_fa] + long_period_window_beds,
                outputs=[homopolymer_bed, long_period_bed, excluded_bed],
                memory="8g",
                walltime="02:00:00",
            ) << f"""
            set -e
            mkdir -p {excluded_dir}
            python {scripts}/find_homopolymers.py {contig_fa} {homopolymer_min_length} > {homopolymer_bed}.tmp
            python {scripts}/merge_regions.py {" ".join(long_period_window_beds)} > {long_period_bed}.tmp
            python {scripts}/merge_regions.py {homopolymer_bed}.tmp {long_period_bed}.tmp > {excluded_bed}.tmp
            mv {homopolymer_bed}.tmp {homopolymer_bed}
            mv {long_period_bed}.tmp {long_period_bed}
            mv {excluded_bed}.tmp {excluded_bed}
            """

        # Step 2 — per model: concatenate the per-contig BEDs into one BED for the
        # sample, and summarise it into a small table to plot locally. The summary
        # counts tracts per (canonical_motif, k, copies, trf_purity, model_purity,
        # model, species, sample), i.e. BED columns 4-11, and goes to the results
        # path; the full BED stays under output.
        #
        # A second table gives what overlap resolution leaves out: the bases in
        # a tract before it, those in the kept tracts, and the difference, the
        # bases of dropped tracts that no kept tract covers. The kept tracts of
        # a model never overlap, so their lengths add up.
        for model in models:
            model_beds = [path("tandem_repeats", contig=contig, model=model) for contig in long_contigs]
            bases_beds = [path("repeat_bases", contig=contig, model=model) for contig in long_contigs]
            bed_list = " ".join(model_beds)
            sample_bed = path("tandem_repeats_sample", model=model)
            summary = path("repeats_summary", model=model)
            dropped_summary = path("repeats_dropped_summary", model=model)

            gwf.target(
                f"cat_{model}_{species}_{sample}_{assembly_type}",
                inputs=model_beds + bases_beds,
                outputs=[sample_bed, summary, dropped_summary],
            ) << f"""
            mkdir -p {outdir} {resultsdir}
            awk 'FNR==1 && NR==1 {{print}} FNR>1 {{print}}' {bed_list} > {sample_bed}

            printf "canonical_motif\\tk\\tcopies\\ttrf_purity\\tmodel_purity\\tmodel\\tspecies\\tsample\\tcount\\n" > {summary}
            awk -F'\\t' -v OFS='\\t' 'NR>1 {{n[$4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11]++}}
                END {{for (key in n) print key, n[key]}}' {sample_bed} \\
                | sort -k2,2n -k1,1 -k3,3g >> {summary}

            before=$(awk '{{bp += $3 - $2}} END {{print bp + 0}}' {" ".join(bases_beds)})
            kept=$(awk 'NR > 1 {{bp += $3 - $2}} END {{print bp + 0}}' {sample_bed})
            printf "model\\tbp_before\\tbp_kept\\tbp_dropped\\tshare_dropped\\tspecies\\tsample\\n" > {dropped_summary}
            awk -v OFS='\\t' -v before=$before -v kept=$kept 'BEGIN {{dropped = before - kept
                print "{model}", before, kept, dropped, (before > 0 ? dropped / before : 0), "{species}", "{sample}"}}' \\
                >> {dropped_summary}
            """

        # Step 3c — the excluded group per sample: for the homopolymer runs,
        # the long-period tracts and their union, the regions (runs, or merged
        # intervals), their bases, and the share of the contigs' bases.
        contigs_bp = sum(long_contigs.values())
        excluded_summary = path("excluded_summary")
        class_beds = {
            "homopolymer": [path("homopolymers", contig=contig) for contig in long_contigs],
            "long_period": [path("long_periods", contig=contig) for contig in long_contigs],
            "excluded":    [path("excluded", contig=contig) for contig in long_contigs],
        }
        class_commands = "\n        ".join(
            f"awk -v OFS='\\t' '{{n++; bp += $3 - $2}} "
            f"END {{print \"{name}\", n + 0, bp + 0, bp / {contigs_bp}, \"{species}\", \"{sample}\"}}' "
            f"{' '.join(beds)} >> {excluded_summary}.tmp"
            for name, beds in class_beds.items())

        gwf.target(
            f"excluded_summary_{species}_{sample}_{assembly_type}",
            inputs=[bed for beds in class_beds.values() for bed in beds],
            outputs=[excluded_summary],
            memory="2g",
            walltime="00:30:00",
        ) << f"""
        set -e
        mkdir -p {resultsdir}
        printf "class\\tregions\\tbp\\tshare\\tspecies\\tsample\\n" > {excluded_summary}.tmp
        {class_commands}
        mv {excluded_summary}.tmp {excluded_summary}
        """
