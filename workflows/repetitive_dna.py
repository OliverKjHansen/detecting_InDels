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
#   M2 — mreps, custom parameters (not implemented yet)
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
scripts           = os.path.abspath("scripts")
mreps             = expanduser(config["tools"]["mreps"])
trf               = expanduser(config["tools"]["trf"])

# Tandem repeat periods to scan.
periods = [1, 2, 3, 4, 5, 6]
periods_str = " ".join(str(k) for k in periods)

# M0 — mreps minimum exponent, i.e. at least this many full copies.
mreps_exp = 3

# M1 — TRF parameters: match mismatch indel PM PI minscore.
# With match weight 2 a perfect tract needs >= minscore/2 bp: ~25 bp.
trf_params = "2 7 7 80 10 50"

models = ["M0", "M1"]

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
                # parse_mreps_canonical.py writes the 11-column M0 BED directly
                python {scripts}/parse_mreps_canonical.py \
                    $WORK/${{k}}_mreps.txt {contig} $k {species} {sample} \
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

        # Step 2 — concatenate the per-contig BEDs into one file per model per sample.
        for model in models:
            model_beds = [f"{outdir}/{contig}_{model}.bed" for contig in long_contigs]
            bed_list = " ".join(model_beds)

            gwf.target(
                f"cat_{model}_{sample}_{assembly_type}",
                inputs=model_beds,
                outputs=[f"{outdir}/{model}.bed"],
            ) << f"""
            mkdir -p {outdir}
            awk 'FNR==1 && NR==1 {{print}} FNR>1 {{print}}' {bed_list} > {outdir}/{model}.bed
            """
