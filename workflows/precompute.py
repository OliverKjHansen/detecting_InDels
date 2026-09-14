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
tissue_dict       = config["tissue_samples"]
assembly_type     = config["assembly_type"]
min_contig_length = config["min_contig_length"]
inpath_fa         = expanduser(config["paths"]["fasta"])
inpath_bam        = expanduser(config["paths"]["bam"])
inpath_labels     = expanduser(config["paths"]["read_labels"])
output            = os.path.abspath(config["paths"]["output"])
scripts           = os.path.abspath("scripts")

# should properly go to config 
chunk_size = 40_000_000


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

        # Step 3 — unresolved haplotypes, per contig, split into windows.
        #
        # The windows [0, chunk_size), [chunk_size, 2*chunk_size), ... tile the
        # contig, so every position is scanned by exactly one job. A region is
        # built from the full span of its reads, so it can reach past its
        # window; the pieces from neighbouring windows then overlap, and
        # bedtools merge joins them back into one region.
        chunk_dir = f"{outdir}/unresolved_haplotypes/chunks"

        for contig, length in contig_dict.items():
            length = int(length)
            if length <= min_contig_length:
                continue

            unresolved = f"{outdir}/unresolved_haplotypes/{contig}.bed"
            chunk_beds = []

            # Step 3a — one job per window.
            for start in range(0, length, chunk_size):
                end = min(start + chunk_size, length)
                chunk_bed = f"{chunk_dir}/{contig}_{start}_{end}.bed"
                chunk_beds.append(chunk_bed)

                gwf.target(
                    f"unresolved_haplotype_{sample}_{contig}_{start}_{assembly_type}",
                    inputs=[assembly, alignment, coverage],
                    outputs=[chunk_bed],
                ) << f"""
                mkdir -p {chunk_dir}
                python {scripts}/resolve_haplotypes.py {contig} {start} {end} {alignment} {assembly} {coverage} > {chunk_bed}
                """

            # Step 3b — join the windows back into one BED per contig.
            chunk_list = " ".join(chunk_beds)

            gwf.target(
                f"merge_unresolved_{sample}_{contig}_{assembly_type}",
                inputs=chunk_beds,
                outputs=[unresolved],
                cores=1,
                memory="4g",
                walltime="00:30:00",
            ) << f"""
            cat {chunk_list} | sort -k1,1 -k2,2n | bedtools merge -i - > {unresolved}
            """

        # Step 4 — read -> tissue table, sorted by read ID. Only samples listed
        # under tissue_samples in config.yaml have read labels.
        if sample in tissue_dict.get(species, []):
            read_labels = f"{inpath_labels}/{species}/{sample}/{sample}_read_labels.txt"
            tissue      = f"{outdir}/sorted_tissue_classification.tsv"

            gwf.target(
                f"sorting_tissue_classification_{sample}_{assembly_type}",
                inputs=[read_labels],
                outputs=[tissue],
            ) << f"""
            mkdir -p {outdir}
            printf "readID\\ttissue\\n" > {tissue}
            awk -v OFS="\\t" '{{gsub(/[ \\t]+/, "\\t"); print $1, $2}}' {read_labels} | sort -k1,1 >> {tissue}
            """
