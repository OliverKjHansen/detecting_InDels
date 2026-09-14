import os
from os.path import expanduser

import yaml
from gwf import Workflow

# ---------------------------------------------------------------------------
# Configuration, read from config.yaml. Run gwf from the repository root:
#
#   gwf -f workflows/precompute.py run
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
inpath_bam        = expanduser(config["paths"]["bam"])
output            = os.path.abspath(config["paths"]["output"])
scripts           = os.path.abspath("scripts")


# Resource requests.
gwf = Workflow(
    defaults={
        'cores': 4,
        'memory': '40g',
        'walltime': '12:00:00',
        'account': 'mutationalscanning'})


for species, samples in sample_dict.items():
    for sample in samples:
        assembly  = f"{inpath_fa}/{species}/{sample}/{assembly_type}_assembly/{sample}_{assembly_type}.fa"
        alignment = f"{inpath_bam}/{species}/{sample}/kinetics/{sample}_kinetics_{assembly_type}.bam"
        outdir    = f"{output}/precomute/{species}/{sample}"
        coverage  = f"{outdir}/avg_coverage.txt"
        contig_dict = {}

        with open(f"{assembly}.fai") as contigs:
            for line in contigs:
                contig_dict[line.split()[0]] = line.split()[1]

        # Step 1 — split the assembly into per-contig FASTAs.
        contigs_dir = f"{outdir}/contigs"

        gwf.target(
            f"faSplit_{sample}_{assembly_type}",
            inputs=[assembly],
            outputs=[f"{contigs_dir}/{contig}.fa" for contig, length in contig_dict.items() if int(length) > min_contig_length],
        ) << f"""
        mkdir -p {contigs_dir}
        faSplit byname {assembly} {contigs_dir}/
        """

        # Step 2 — average depth across the positions samtools reports.
        long_contigs_bed = f"{outdir}/long_contigs.bed"

        gwf.target(
            f"average_coverage_{sample}",
            inputs=[alignment, f"{assembly}.fai"],
            outputs=[long_contigs_bed, coverage],
            cores=8,
            memory="100g",
            walltime="08:00:00",
        ) << f"""
        mkdir -p {outdir}
        awk -v min={min_contig_length} 'BEGIN {{OFS="\\t"}} $2 > min {{print $1, 0, $2}}' {assembly}.fai > {long_contigs_bed}
        samtools depth -b {long_contigs_bed} {alignment} | awk '{{sum += $3; count++}} END {{print sum / count}}' > {coverage}
        """
        # for contig, length in contig_dict.items():
        #     if int(length) > 10**6:
                
        # Step 3 — unresolved haplotypes, per contig.
        # Commented out while debugging steps 1 and 2.
        # for contig in contigs:
        #     unresolved = f"{outdir}/unresolved_haplotypes/{contig}.bed"
        #
        #     gwf.target(
        #         f"unresolved_haplotype_{sample}_{contig}_{assembly_type}",
        #         inputs=[assembly, alignment, coverage],
        #         outputs=[unresolved],
        #     ) << f"""
        #     mkdir -p {outdir}/unresolved_haplotypes
        #     python {scripts}/resolve_haplotypes.py {contig} {alignment} {assembly} {coverage} > {unresolved}
        #     """
