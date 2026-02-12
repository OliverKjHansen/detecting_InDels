from os.path import expanduser
from gwf import *

sample_dict = {
                'chimp':['yoran'],
                'macaque':['jako'],
                'marmoset':['marmoset', 'lambiek']
               }

assembly_type = "diploid"
inpath_fa = expanduser("~/mutationalscanning/DerivedData/fasta/HiFi")
inpath_bam = expanduser("~/mutationalscanning/DerivedData/bam/HiFi")
scripts = expanduser("~/mutationalscanning/Workspaces/oliver/endjoining/scripts")

gwf = Workflow(
    defaults={
        'cores': 4,
        'memory': '10g',
        'walltime': '48:00:00',
        'account' : 'mutationalscanning'})


for species, samples in sample_dict.items():
    for sample in samples:
        index_assembly = f"{inpath_fa}/{species}/{sample}/diploid_assembly/{sample}_diploid.fa.fai"
        assembly = f"{inpath_fa}/{species}/{sample}/diploid_assembly/{sample}_diploid.fa"
        contig_dict = {}
        alignment = f"{inpath_bam}/{species}/{sample}/kinetics/{sample}_kinetics_diploid.bam"
        avg_coverage_file = expanduser(f"~/mutationalscanning/Workspaces/oliver/microsatellites/ratemicrosat/coverage/{sample}_avg_coverage.txt")
        tissue_file = f"/faststorage/project/mutationalscanning/Workspaces/oliver/microsatellites/marmosetanalysis/{species}/{sample}_sorted_tissue_classification.tsv"

        with open(index_assembly) as contigs:
            for line in contigs:
                contig_dict[line.split()[0]] = line.split()[1]

            for contig, length in contig_dict.items():
                if (int(length) > 10**6):
                    microsatfile = f"/faststorage/project/mutationalscanning/Workspaces/oliver/microsatellites/ratemicrosat/{species}/{sample}/{contig}_trs_tidy.tsv"
                    unresolved_bed = f"/faststorage/project/mutationalscanning/Workspaces/oliver/microsatellites/ratemicrosat/{species}/{sample}/{contig}_unresolved_haplotypes.bed"

                    # Step 1: Build filter BED and call small indels per contig
                    gwf.target(
                        f"SmallIndels_{contig}_{sample}_diploid",
                        inputs=[assembly, alignment, avg_coverage_file, tissue_file, microsatfile],
                        outputs=[f"output/small_indels/{species}/{sample}/{contig}_filters.tsv",
                                f"output/small_indels/{species}/{sample}/{contig}_small_indels.tsv"]
                        ) << f"""
                        mkdir -p output/small_indels/{species}/{sample}
                        cat <(awk 'OFS="\\t" {{print $1,$2,$3}}' {microsatfile}) {unresolved_bed} | sort -k1,1 -k2,2n | bedtools merge -i - > output/small_indels/{species}/{sample}/{contig}_filters.tsv
                        python {scripts}/small_indels.py output/small_indels/{species}/{sample}/{contig}_filters.tsv {contig} {alignment} {assembly} {avg_coverage_file} {tissue_file} > output/small_indels/{species}/{sample}/{contig}_small_indels.tsv
                        """

                    # Step 2: Extract mutation-supporting reads as FASTQ for haploid remapping
                    gwf.target(
                        f"build_mutation_fastq_{contig}_{sample}_diploid",
                        inputs=[f"output/small_indels/{species}/{sample}/{contig}_small_indels.tsv"],
                        outputs=[f"output/mutation_fastq/{species}/{sample}/{contig}_remapping.fastq"]
                        ) << f"""
                        mkdir -p output/mutation_fastq/{species}/{sample}
                        python {scripts}/build_mutation_fastq.py output/small_indels/{species}/{sample}/{contig}_small_indels.tsv {alignment} {contig} "output/mutation_fastq/{species}/{sample}/{contig}_remapping.fastq"
                        """

                    # Step 3: Count indel opportunities per contig
                    gwf.target(
                        f"opportunities_{contig}_{sample}_diploid",
                        inputs=[assembly, alignment, avg_coverage_file, tissue_file,
                               f"output/small_indels/{species}/{sample}/{contig}_filters.tsv"],
                        outputs=[f"output/opportunities/{species}/{sample}/{contig}_opportunities.tsv"]
                        ) << f"""
                        mkdir -p output/opportunities/{species}/{sample}
                        python {scripts}/opportunities.py output/small_indels/{species}/{sample}/{contig}_filters.tsv {contig} {alignment} {assembly} {avg_coverage_file} {tissue_file} > output/opportunities/{species}/{sample}/{contig}_opportunities.tsv
                        """

            # Step 4: Remap mutation reads to haploid assembly
            gwf.target(
                f"haploid_mapping_{sample}_diploid",
                inputs=[f"output/mutation_fastq/{species}/{sample}/{contig}_remapping.fastq" for contig, length in contig_dict.items() if (int(length) > 10**6)],
                outputs=[f"output/haploid_mapping/{species}_{sample}_haploid_locations.tsv"],
                cores=8,
                memory="100g",
                walltime='24:00:00',
                account="mutationalscanning") << f"""
                mkdir -p output/haploid_mapping
                cat output/mutation_fastq/{species}/{sample}/*_remapping.fastq > output/mutation_fastq/{species}/{sample}/{sample}_combined.fastq
                pbmm2 align {inpath_fa}/{species}/{sample}/haploid_assembly/{sample}_haploid.fa output/mutation_fastq/{species}/{sample}/{sample}_combined.fastq output/haploid_mapping/{species}_{sample}_remapped.bam --preset HIFI --sort -j 12 -J 4 -m 8G
                samtools view output/haploid_mapping/{species}_{sample}_remapped.bam | awk -v OFS='\t' '{{print $1,$2,$3,$4,$5,$6}}' | sort -k1,1 > output/haploid_mapping/{species}_{sample}_haploid_locations.tsv
                """

            # Step 5: Filter indels using haploid mapping concordance
            for contig, length in contig_dict.items():
                if (int(length) > 10**6):
                    gwf.target(
                        f"haploid_filter_{contig}_{sample}_diploid",
                        inputs=[f"output/haploid_mapping/{species}_{sample}_haploid_locations.tsv",
                                f"output/small_indels/{species}/{sample}/{contig}_small_indels.tsv"],
                        outputs=[f"output/small_indels/{species}/{sample}/{contig}_small_indels_haploid_filtered.tsv"]
                        ) << f"""
                        mkdir -p output/small_indels/{species}/{sample}
                        python {scripts}/filtering_haploid_mappping.py output/small_indels/{species}/{sample}/{contig}_small_indels.tsv output/haploid_mapping/{species}_{sample}_haploid_locations.tsv > output/small_indels/{species}/{sample}/{contig}_small_indels_haploid_filtered_tmp.tsv
                        awk 'NR==1 {{print $0 "\t" "species" "\t" "sample"}}; NR>1 {{print $0 "\t" "{species}" "\t" "{sample}"}}' output/small_indels/{species}/{sample}/{contig}_small_indels_haploid_filtered_tmp.tsv > output/small_indels/{species}/{sample}/{contig}_small_indels_haploid_filtered.tsv
                        rm -f output/small_indels/{species}/{sample}/{contig}_small_indels_haploid_filtered_tmp.tsv
                        """

            # Step 6: Merge per-contig results and generate plots
            gwf.target(
                f"merge_and_plot_{sample}_diploid",
                inputs=[f"output/small_indels/{species}/{sample}/{contig}_small_indels_haploid_filtered.tsv" for contig, length in contig_dict.items() if (int(length) > 10**6)] +
                       [f"output/opportunities/{species}/{sample}/{contig}_opportunities.tsv" for contig, length in contig_dict.items() if (int(length) > 10**6)],
                outputs=[f"output/final/{sample}_haploid_combined.tsv",
                         f"output/final/{sample}_opportunities_combined.tsv",
                         f"output/plots/{sample}_indel_plots.pdf"]
                ) << f"""
                mkdir -p output/final output/plots
                awk 'FNR==1 && NR==1 || FNR>1' output/small_indels/{species}/{sample}/*_small_indels_haploid_filtered.tsv > output/final/{sample}_haploid_combined.tsv
                awk 'FNR==1 && NR==1 || FNR>1' output/opportunities/{species}/{sample}/*_opportunities.tsv > output/final/{sample}_opportunities_combined.tsv
                Rscript {scripts}/plot_indels.R output/final/{sample}_haploid_combined.tsv output/final/{sample}_opportunities_combined.tsv output/plots/{sample}_indel_plots.pdf
                """

            # Step 7: Apply quality and haploid-concordance filters with reporting
            gwf.target(
                f"initial_filtering_{sample}_diploid",
                inputs=[f"output/final/{sample}_haploid_combined.tsv",
                        f"output/final/{sample}_opportunities_combined.tsv"],
                outputs=[f"output/final/{sample}_combined_filtered.tsv",
                        f"output/final/{sample}_opportunities_combined_filtered.tsv"],
                cores=8,
                memory="100g",
                walltime='1:00:00',
                account="mutationalscanning"
                ) << f"""
                mkdir -p output/final
                python {scripts}/filtering_small_indels.py output/final/{sample}_haploid_combined.tsv -o output/final/{sample}_combined_filtered.tsv --haploid-check True --min-quality 93 --min-read-end-distance 100 --report output/final/{sample}_filter_report.tsv
                python {scripts}/filtering_opportunities.py output/final/{sample}_opportunities_combined.tsv -o output/final/{sample}_opportunities_combined_filtered.tsv --min-kmerq 93
                """