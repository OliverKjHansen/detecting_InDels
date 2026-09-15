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

# M0 and M2 — mreps minimum exponent, i.e. at least this many full copies.
mreps_exp = 3

# M2 — mreps resolution: substitutions allowed between adjacent copies.
m2_res = 1

# M1 — TRF parameters: match mismatch indel PM PI minscore.
# With match weight 2 a perfect tract needs >= minscore/2 bp: ~25 bp.
trf_params = "2 7 7 80 10 50"

models = ["M0", "M1", "M2"]

# Resource requests.
gwf = Workflow(
    defaults={
        'cores': 4,
        'memory': '40g',
        'walltime': '12:00:00',
        'account': 'mutationalscanning'})


for species, samples in sample_dict.items():
    for sample in samples:
        assembly    = f"{inpath_fa}/{species}/{sample}/{assembly_type}_assembly/{sample}_{assembly_type}.fa"
        contigs_dir = f"{output}/precomute/{species}/{sample}/contigs"
        outdir      = f"{output}/repetitive_dna/{species}/{sample}"
        resultsdir  = f"{results}/repetitive_dna/{species}/{sample}"
        contig_dict = {}

        with open(f"{assembly}.fai") as contigs:
            for line in contigs:
                contig_dict[line.split()[0]] = line.split()[1]

        long_contigs = [contig for contig, length in contig_dict.items() if int(length) > min_contig_length]

        # Step 1 — one job per model per contig. Each job loops over every k in
        # bash; per-k intermediates live in a scratch dir removed on exit, so
        # only the per-contig BED reaches disk.
        for contig in long_contigs:
            contig_fa = f"{contigs_dir}/{contig}.fa"
            m0_k_beds = " ".join(f"$WORK/{k}_M0.bed" for k in periods)
            m1_k_beds = " ".join(f"$WORK/{k}_M1.bed" for k in periods)
            m2_k_beds = " ".join(f"$WORK/{k}_M2.bed" for k in periods)

            # Step 1a — M0: mreps exact/pure repeats.
            gwf.target(
                f"M0_mreps_{contig}_{sample}_{assembly_type}",
                inputs=[contig_fa],
                outputs=[f"{outdir}/{contig}_M0.bed"],
            ) << f"""
            mkdir -p {outdir}
            WORK=$(mktemp -d)
            trap 'rm -rf "$WORK"' EXIT

            for k in {periods_str}; do
                # -exp N:      at least N full copies (exponent = length/period)
                # -allowsmall: allow tracts shorter than the random-occurrence threshold
                {mreps} -minperiod $k -maxperiod $k -exp {mreps_exp} -allowsmall -fasta {contig_fa} \
                    > $WORK/${{k}}_mreps.txt
                python {scripts}/parse_mreps_canonical.py \
                    $WORK/${{k}}_mreps.txt {contig} $k {species} {sample} M0 \
                    > $WORK/${{k}}_M0.bed
            done

            # Concatenate all k - header retained from the first file only.
            awk 'FNR==1 && NR==1 {{print}} FNR>1 {{print}}' {m0_k_beds} > $WORK/all_M0.bed

            # Overlapping tracts: keep the smallest k, drop the others.
            python {scripts}/resolve_repeat_overlaps.py $WORK/all_M0.bed > {outdir}/{contig}_M0.bed
            """

            # Step 1b — M1: TRF standard parameters.
            gwf.target(
                f"M1_trf_{contig}_{sample}_{assembly_type}",
                inputs=[contig_fa],
                outputs=[f"{outdir}/{contig}_M1.bed"],
            ) << f"""
            mkdir -p {outdir}
            WORK=$(mktemp -d)
            trap 'rm -rf "$WORK"' EXIT

            for k in {periods_str}; do
                # -ngs -> machine-readable stdout, -h -> suppress HTML
                # TRF has no -minperiod, so awk keeps only rows with period == k
                {trf} {contig_fa} {trf_params} $k -ngs -h \
                    | awk -v k=$k '/^@/{{print}} !/^@/{{if ($3==k) print}}' \
                    > $WORK/${{k}}_trf.txt
                python {scripts}/parse_trf.py \
                    $WORK/${{k}}_trf.txt $k {species} {sample} M1 \
                    > $WORK/${{k}}_M1.bed
            done

            # Concatenate all k - header retained from the first file only.
            awk 'FNR==1 && NR==1 {{print}} FNR>1 {{print}}' {m1_k_beds} > $WORK/all_M1.bed

            # Overlapping tracts: keep the smallest k, drop the others.
            python {scripts}/resolve_repeat_overlaps.py $WORK/all_M1.bed > {outdir}/{contig}_M1.bed
            """

            # Step 1c — M2: mreps with resolution m2_res.
            gwf.target(
                f"M2_mreps_{contig}_{sample}_{assembly_type}",
                inputs=[contig_fa],
                outputs=[f"{outdir}/{contig}_M2.bed"],
            ) << f"""
            mkdir -p {outdir}
            WORK=$(mktemp -d)
            trap 'rm -rf "$WORK"' EXIT

            for k in {periods_str}; do
                # -res N: allow up to N substitutions between adjacent copies
                {mreps} -minperiod $k -maxperiod $k -exp {mreps_exp} -res {m2_res} -allowsmall -fasta {contig_fa} \
                    > $WORK/${{k}}_mreps.txt
                python {scripts}/parse_mreps_canonical.py \
                    $WORK/${{k}}_mreps.txt {contig} $k {species} {sample} M2 \
                    > $WORK/${{k}}_M2.bed
            done

            # Concatenate all k - header retained from the first file only.
            awk 'FNR==1 && NR==1 {{print}} FNR>1 {{print}}' {m2_k_beds} > $WORK/all_M2.bed

            # Overlapping tracts: keep the smallest k, drop the others.
            python {scripts}/resolve_repeat_overlaps.py $WORK/all_M2.bed > {outdir}/{contig}_M2.bed
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
