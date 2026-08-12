from os.path import expanduser
from gwf import *

sample_dict = {
                'chimp':    ['yoran'],
                'macaque':  ['jako'],
                'marmoset': ['marmoset', 'lambiek']
               }

models = ['M0', 'M1']

assembly_type = "diploid"
inpath_fa  = expanduser("~/mutationalscanning/DerivedData/fasta/HiFi")
inpath_bam = expanduser("~/mutationalscanning/DerivedData/bam/HiFi")
scripts    = expanduser("~/mutationalscanning/Workspaces/oliver/endjoining/scripts")

gwf = Workflow(
    defaults={
        'cores':    4,
        'memory':   '10g',
        'walltime': '48:00:00',
        'account':  'mutationalscanning'})


for species, samples in sample_dict.items():
    for sample in samples:
        index_assembly = f"{inpath_fa}/{species}/{sample}/diploid_assembly/{sample}_diploid.fa.fai"
        assembly       = f"{inpath_fa}/{species}/{sample}/diploid_assembly/{sample}_diploid.fa"
        alignment      = f"{inpath_bam}/{species}/{sample}/kinetics/{sample}_kinetics_diploid.bam"
        avg_coverage_file = expanduser(
            f"~/mutationalscanning/Workspaces/oliver/microsatellites/ratemicrosat/coverage/{sample}_avg_coverage.txt"
        )
        tissue_file = (
            f"/faststorage/project/mutationalscanning/Workspaces/oliver/"
            f"microsatellites/marmosetanalysis/{species}/{sample}_sorted_tissue_classification.tsv"
        )

        contig_dict = {}
        with open(index_assembly) as contigs:
            for line in contigs:
                contig_dict[line.split()[0]] = line.split()[1]

        for model in models:

            for contig, length in contig_dict.items():
                if int(length) <= 10**6:
                    continue

                tr_bed = f"output/tandem_repeats/{species}/{sample}/{contig}_{model}.bed"

                # Step 1: Detect period-multiple InDels inside TR regions
                gwf.target(
                    f"TRIndels_{contig}_{sample}_{model}",
                    inputs=[assembly, alignment, avg_coverage_file,
                            tissue_file, tr_bed],
                    outputs=[f"output/tandem_repeat_indels/{species}/{sample}/{contig}_{model}_tr_indels.tsv"]
                ) << f"""
                mkdir -p output/tandem_repeat_indels/{species}/{sample}
                python {scripts}/tandem_repeat_indels.py \
                    {tr_bed} {contig} {alignment} {assembly} \
                    {avg_coverage_file} {tissue_file} \
                    > output/tandem_repeat_indels/{species}/{sample}/{contig}_{model}_tr_indels.tsv
                """

                # Step 2: Extract mutation-supporting reads as FASTQ for haploid remapping
                gwf.target(
                    f"TRBuildFastq_{contig}_{sample}_{model}",
                    inputs=[f"output/tandem_repeat_indels/{species}/{sample}/{contig}_{model}_tr_indels.tsv"],
                    outputs=[f"output/tandem_repeat_mutation_fastq/{species}/{sample}/{contig}_{model}_remapping.fastq"]
                ) << f"""
                mkdir -p output/tandem_repeat_mutation_fastq/{species}/{sample}
                python {scripts}/build_mutation_fastq.py \
                    output/tandem_repeat_indels/{species}/{sample}/{contig}_{model}_tr_indels.tsv \
                    {alignment} {contig} \
                    output/tandem_repeat_mutation_fastq/{species}/{sample}/{contig}_{model}_remapping.fastq
                """

                # Step 3: Count TR opportunities (parallel with Step 2)
                gwf.target(
                    f"TROpportunities_{contig}_{sample}_{model}",
                    inputs=[assembly, alignment, avg_coverage_file,
                            tissue_file, tr_bed],
                    outputs=[f"output/tandem_repeat_opportunities/{species}/{sample}/{contig}_{model}_tr_opportunities.tsv"]
                ) << f"""
                mkdir -p output/tandem_repeat_opportunities/{species}/{sample}
                python {scripts}/tandem_repeat_opportunities.py \
                    {tr_bed} {contig} {alignment} {assembly} \
                    {avg_coverage_file} {tissue_file} \
                    > output/tandem_repeat_opportunities/{species}/{sample}/{contig}_{model}_tr_opportunities.tsv
                """

            large_contigs = [c for c, l in contig_dict.items() if int(l) > 10**6]

            # Step 4: Remap mutation reads to haploid assembly
            gwf.target(
                f"TRHaploidMapping_{sample}_{model}",
                inputs=[f"output/tandem_repeat_mutation_fastq/{species}/{sample}/{contig}_{model}_remapping.fastq"
                        for contig in large_contigs],
                outputs=[f"output/tandem_repeat_haploid_mapping/{species}_{sample}_{model}_haploid_locations.tsv"],
                cores=8,
                memory="100g",
                walltime='24:00:00',
                account="mutationalscanning"
            ) << f"""
            mkdir -p output/tandem_repeat_haploid_mapping
            cat output/tandem_repeat_mutation_fastq/{species}/{sample}/*_{model}_remapping.fastq \
                > output/tandem_repeat_mutation_fastq/{species}/{sample}/{sample}_{model}_combined.fastq
            pbmm2 align \
                {inpath_fa}/{species}/{sample}/haploid_assembly/{sample}_haploid.fa \
                output/tandem_repeat_mutation_fastq/{species}/{sample}/{sample}_{model}_combined.fastq \
                output/tandem_repeat_haploid_mapping/{species}_{sample}_{model}_remapped.bam \
                --preset HIFI --sort -j 12 -J 4 -m 8G
            samtools view output/tandem_repeat_haploid_mapping/{species}_{sample}_{model}_remapped.bam \
                | awk -v OFS='\\t' '{{print $1,$2,$3,$4,$5,$6}}' \
                | sort -k1,1 \
                > output/tandem_repeat_haploid_mapping/{species}_{sample}_{model}_haploid_locations.tsv
            """

            # Step 5: Haploid concordance filter (per contig)
            for contig, length in contig_dict.items():
                if int(length) <= 10**6:
                    continue

                gwf.target(
                    f"TRHaploidFilter_{contig}_{sample}_{model}",
                    inputs=[f"output/tandem_repeat_haploid_mapping/{species}_{sample}_{model}_haploid_locations.tsv",
                            f"output/tandem_repeat_indels/{species}/{sample}/{contig}_{model}_tr_indels.tsv"],
                    outputs=[f"output/tandem_repeat_indels/{species}/{sample}/{contig}_{model}_tr_indels_haploid_filtered.tsv"]
                ) << f"""
                mkdir -p output/tandem_repeat_indels/{species}/{sample}
                python {scripts}/filtering_haploid_mappping.py \
                    output/tandem_repeat_indels/{species}/{sample}/{contig}_{model}_tr_indels.tsv \
                    output/tandem_repeat_haploid_mapping/{species}_{sample}_{model}_haploid_locations.tsv \
                    > output/tandem_repeat_indels/{species}/{sample}/{contig}_{model}_tr_indels_haploid_filtered_tmp.tsv
                awk 'NR==1 {{print $0 "\\t" "species" "\\t" "sample"}}; NR>1 {{print $0 "\\t" "{species}" "\\t" "{sample}"}}' \
                    output/tandem_repeat_indels/{species}/{sample}/{contig}_{model}_tr_indels_haploid_filtered_tmp.tsv \
                    > output/tandem_repeat_indels/{species}/{sample}/{contig}_{model}_tr_indels_haploid_filtered.tsv
                rm -f output/tandem_repeat_indels/{species}/{sample}/{contig}_{model}_tr_indels_haploid_filtered_tmp.tsv
                """

            # Step 6: Merge and plot
            gwf.target(
                f"TRMergeAndPlot_{sample}_{model}",
                inputs=(
                    [f"output/tandem_repeat_indels/{species}/{sample}/{contig}_{model}_tr_indels_haploid_filtered.tsv"
                     for contig in large_contigs]
                    + [f"output/tandem_repeat_opportunities/{species}/{sample}/{contig}_{model}_tr_opportunities.tsv"
                       for contig in large_contigs]
                ),
                outputs=[
                    f"output/final/{sample}_{model}_tr_haploid_combined.tsv",
                    f"output/final/{sample}_{model}_tr_opportunities_combined.tsv",
                    f"output/plots/{sample}_{model}_tr_indel_plots.pdf"
                ],
                cores=8,
                memory="100g",
                walltime='24:00:00',
                account="mutationalscanning"
            ) << f"""
            mkdir -p output/final output/plots
            awk 'FNR==1 && NR==1 || FNR>1' \
                output/tandem_repeat_indels/{species}/{sample}/*_{model}_tr_indels_haploid_filtered.tsv \
                > output/final/{sample}_{model}_tr_haploid_combined.tsv
            awk 'FNR==1 && NR==1 || FNR>1' \
                output/tandem_repeat_opportunities/{species}/{sample}/*_{model}_tr_opportunities.tsv \
                > output/final/{sample}_{model}_tr_opportunities_combined.tsv
            Rscript {scripts}/plot_indels.R \
                output/final/{sample}_{model}_tr_haploid_combined.tsv \
                output/final/{sample}_{model}_tr_opportunities_combined.tsv \
                output/plots/{sample}_{model}_tr_indel_plots.pdf
            """

            # Step 7: Final filtering
            gwf.target(
                f"TRFinalFiltering_{sample}_{model}",
                inputs=[
                    f"output/final/{sample}_{model}_tr_haploid_combined.tsv",
                    f"output/final/{sample}_{model}_tr_opportunities_combined.tsv"
                ],
                outputs=[
                    f"output/final/{sample}_{model}_tr_combined_filtered.tsv",
                    f"output/final/{sample}_{model}_tr_opportunities_combined_filtered.tsv"
                ],
                cores=4,
                memory="20g",
                walltime='1:00:00',
                account="mutationalscanning"
            ) << f"""
            mkdir -p output/final
            python {scripts}/filtering_small_indels.py \
                output/final/{sample}_{model}_tr_haploid_combined.tsv \
                -o output/final/{sample}_{model}_tr_combined_filtered.tsv \
                --haploid-check True \
                --min-quality 93 \
                --min-read-end-distance 100 \
                --report output/final/{sample}_{model}_tr_filter_report.tsv
            python {scripts}/filtering_opportunities.py \
                output/final/{sample}_{model}_tr_opportunities_combined.tsv \
                -o output/final/{sample}_{model}_tr_opportunities_combined_filtered.tsv \
                --min-kmerq 93
            """
