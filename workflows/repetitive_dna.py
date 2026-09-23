import math
import os
from os.path import expanduser

import yaml
from gwf import Workflow

# ---------------------------------------------------------------------------
# Repetitive DNA annotation. Annotates tandem repeats per contig for every
# period k, under three models:
#
#   M0 — mreps, exact/pure repeats
#   M1 — TRF, standard scoring from RepeatMasker
#   M2 — mreps, resolution 1: one substitution allowed between adjacent copies
#
# Configuration is read from config.yaml. Run gwf from the repository root,
# after precompute.py has split the assemblies into contigs:
#
#   gwf -f workflows/precompute.py run
#   gwf -f workflows/repetitive_dna.py run
#
# output/ and scripts/ are made absolute from there, because gwf runs the
# recipes from workflows/.
# ---------------------------------------------------------------------------

with open("config.yaml") as handle:
    config = yaml.safe_load(handle)

sample_dict       = config["samples"]
assembly_type     = config["assembly_type"]
min_contig_length = config["min_contig_length"]
inpath_fa         = expanduser(config["paths"]["fasta"])
output            = os.path.abspath(config["paths"]["output"])
results           = os.path.abspath(config["paths"]["results"])
scripts           = os.path.abspath("scripts")
mreps             = expanduser(config["tools"]["mreps"])
trf               = expanduser(config["tools"]["trf"])

# Tandem repeat periods to scan. Homopolymers (k=1) are not annotated.
periods = [2, 3, 4, 5, 6]
periods_str = " ".join(str(k) for k in periods)
max_period = max(periods)

# M0 and M2 — mreps minimum exponent, i.e. at least this many full copies.
mreps_exp = 3

# M2 — mreps resolution: substitutions allowed between adjacent copies.
m2_res = 1

# M1 — TRF parameters: match mismatch indel PM PI minscore.
# With match weight 2 a perfect tract needs >= minscore/2 bp: ~25 bp.
trf_params = "2 7 7 80 10 50"

# M1 — TRF runs on windows of a contig, one job per window, because its runtime
# explodes in long satellite arrays and a whole chromosome can run for days.
# Windows overlap by trf_overlap so a tract crossing a boundary is still seen
# whole; a tract longer than the overlap can be cut in two.
trf_window  = 5_000_000
trf_overlap = 1_000_000

# M1 — TRF aligns a repeat over at most -l million bp (default 2). In a uniform
# array longer than that, e.g. a 4 Mb centromere, it keeps restarting the
# alignment and a window runs for days instead of minutes. So the limit covers
# everything a window job fetches, and no array can exceed it. TRF allocates
# ~0.6 GB per million up front: -l 12 for 10 Mb windows is ~7.5 GB.
trf_max_length = math.ceil((trf_window + 2 * trf_overlap) / 1_000_000)

models = ["M0", "M1", "M2"]

# Resource requests. mreps and TRF are single-threaded, so one core per job.
gwf = Workflow(
    defaults={
        'cores': 1,
        'memory': '10g',
        'walltime': '8:00:00',
        'account': 'mutationalscanning'})


for species, samples in sample_dict.items():
    for sample in samples:
        assembly    = f"{inpath_fa}/{species}/{sample}/{assembly_type}_assembly/{sample}_{assembly_type}.fa"
        contigs_dir = f"{output}/precomute/{species}/{sample}/contigs"
        outdir      = f"{output}/repetitive_dna/{species}/{sample}"
        per_k_dir   = f"{outdir}/per_k"
        resultsdir  = f"{results}/repetitive_dna/{species}/{sample}"
        contig_dict = {}

        with open(f"{assembly}.fai") as contigs:
            for line in contigs:
                contig_dict[line.split()[0]] = line.split()[1]

        long_contigs = {contig: int(length) for contig, length in contig_dict.items() if int(length) > min_contig_length}

        for contig, length in long_contigs.items():
            contig_fa = f"{contigs_dir}/{contig}.fa"

            # Step 1a — one job per model per period k. Each job writes its BED
            # to a .tmp file and renames it only when every command succeeded, so
            # a failed job never leaves a finished-looking output behind.
            for k in periods:

                # M0: mreps exact/pure repeats.
                # -exp N:      at least N full copies (exponent = length/period)
                # -allowsmall: allow tracts shorter than the random-occurrence threshold
                m0_bed = f"{per_k_dir}/{contig}_M0_k{k}.bed"
                gwf.target(
                    f"M0_mreps_k{k}_{contig}_{sample}_{assembly_type}",
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
                m2_bed = f"{per_k_dir}/{contig}_M2_k{k}.bed"
                gwf.target(
                    f"M2_mreps_k{k}_{contig}_{sample}_{assembly_type}",
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
                window_beds = [f"{per_k_dir}/{contig}_M1_k{k}_{core_start}.bed" for k in periods]
                m1_window_beds += window_beds

                gwf.target(
                    f"M1_trf_{contig}_{core_start}_{sample}_{assembly_type}",
                    inputs=[assembly, f"{assembly}.fai"],
                    outputs=window_beds,
                    walltime="24:00:00",
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
                        > {per_k_dir}/{contig}_M1_k${{k}}_{core_start}.bed.tmp
                    mv {per_k_dir}/{contig}_M1_k${{k}}_{core_start}.bed.tmp {per_k_dir}/{contig}_M1_k${{k}}_{core_start}.bed
                done
                """

            # Step 1b — per model: join the per-k BEDs of the contig (header kept
            # from the first file only) and resolve overlapping tracts, keeping
            # the smallest k. M1 also has one set of BEDs per window.
            k_beds_by_model = {
                "M0": [f"{per_k_dir}/{contig}_M0_k{k}.bed" for k in periods],
                "M1": m1_window_beds,
                "M2": [f"{per_k_dir}/{contig}_M2_k{k}.bed" for k in periods],
            }

            for model in models:
                k_beds = k_beds_by_model[model]
                k_bed_list = " ".join(k_beds)
                contig_bed = f"{outdir}/{contig}_{model}.bed"

                gwf.target(
                    f"{model}_merge_k_{contig}_{sample}_{assembly_type}",
                    inputs=k_beds,
                    outputs=[contig_bed],
                    memory="16g",
                    walltime="02:00:00",
                ) << f"""
                set -e
                mkdir -p {outdir}
                WORK=$(mktemp -d)
                trap 'rm -rf "$WORK"' EXIT
                awk 'FNR==1 && NR==1 {{print}} FNR>1 {{print}}' {k_bed_list} > $WORK/all.bed
                python {scripts}/resolve_repeat_overlaps.py $WORK/all.bed > {contig_bed}.tmp
                mv {contig_bed}.tmp {contig_bed}
                """

        # Step 2 — per model: concatenate the per-contig BEDs into one BED for the
        # sample, and summarise it into a small table to plot locally. The summary
        # counts tracts per (canonical_motif, k, copies, trf_purity, model_purity,
        # model, species, sample), i.e. BED columns 4-11, and goes to the results
        # path; the full BED stays under output.
        for model in models:
            model_beds = [f"{outdir}/{contig}_{model}.bed" for contig in long_contigs]
            bed_list = " ".join(model_beds)
            sample_bed = f"{outdir}/{model}.bed"
            summary = f"{resultsdir}/{model}_summary.tsv"

            gwf.target(
                f"cat_{model}_{sample}_{assembly_type}",
                inputs=model_beds,
                outputs=[sample_bed, summary],
            ) << f"""
            mkdir -p {outdir} {resultsdir}
            awk 'FNR==1 && NR==1 {{print}} FNR>1 {{print}}' {bed_list} > {sample_bed}

            printf "canonical_motif\\tk\\tcopies\\ttrf_purity\\tmodel_purity\\tmodel\\tspecies\\tsample\\tcount\\n" > {summary}
            awk -F'\\t' -v OFS='\\t' 'NR>1 {{n[$4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11]++}}
                END {{for (key in n) print key, n[key]}}' {sample_bed} \\
                | sort -k2,2n -k1,1 -k3,3g >> {summary}
            """
