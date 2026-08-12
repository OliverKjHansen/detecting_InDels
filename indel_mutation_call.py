from os.path import expanduser
from gwf import *
gwf = Workflow()
import pysam
import sys
import csv
from collections import defaultdict

sample_dict = {#'human':["chk"],
            #   'human':['ob006', 'ob007', 'ob008', 'ob009', 'da1', 'la4','tsaed','ob2v','ob003h',"aa_sperm","chk","la","ph","tr"],
                'chimp':['yoran'],
            #   'gorilla':['kipenzi'], # 'mbewe','samson'], empty kinestics assembly
            #   'baboon':['sethi'], #'bs5v'contaminated
            #   'gibbon':['tarzan'],
                'macaque':['jako'], #lenny']#,
                'marmoset':['marmoset', 'lambiek'] 
               }

# classified_samples = ['ob006', 'ob007', 'ob008', 'ob009' #'ob2v','ob003h', # humans
#                       #'stephan', 'carl',#'marlock', 'martin', 'mongo', #chimp #these samples do not have methylation data 
#                       #'mbewe','samson',#'kipenzi', # gorilla
#                       ] 

#methylation_classification = expanduser("~/mutationalscanning/Workspaces/vinod/AllDatasets_germFracDynamics_vs_ll_HiFiReadsLL_UpdateAgreed_testLowRAM_update2_final/ToOliver_germReads")
assembly_type = "diploid"
inpath_fa = expanduser("~/mutationalscanning/DerivedData/fasta/HiFi")
inpath_bam = expanduser("~/mutationalscanning/DerivedData/bam/HiFi")
#tandem_reapeat_finder = expanduser("~/mutationalscanning/Workspaces/oliver/microsatellites/development/mreps")
scripts = expanduser("~/mutationalscanning/Workspaces/oliver/endjoining/scripts")

gwf = Workflow(
    defaults={
        'cores': 4,
        'memory': '10g',
        'walltime': '48:00:00',
        'account' : 'mutationalscanning'})


for species, samples in sample_dict.items():
    #print(f"Species: {species}")
    for sample in samples:
        #print(f"{inpath_fa}/{species}/{sample}/diploid_assembly/{sample}_diploid.fa.fai")
        index_assembly = f"{inpath_fa}/{species}/{sample}/diploid_assembly/{sample}_diploid.fa.fai"
        assembly = f"{inpath_fa}/{species}/{sample}/diploid_assembly/{sample}_diploid.fa"
        contig_dict = {}
<<<<<<< HEAD:workflow.py
        alignemt = f"{inpath_bam}/{species}/{sample}/kinetics/{sample}_kinetics_diploid.bam" # new version
        avg_coverage_file = expanduser(f"~/mutationalscanning/Workspaces/oliver/microsatellites/ratemicrosat/coverage/{sample}_avg_coverage.txt")
        tissue_file = f"/faststorage/project/mutationalscanning/Workspaces/oliver/microsatellites/marmosetanalysis/{species}/{sample}_sorted_tissue_classification.tsv"
=======
        alignment = f"{inpath_bam}/{species}/{sample}/kinetics/{sample}_kinetics_diploid.bam"
        avg_coverage_file = expanduser(f"~/mutationalscanning/Workspaces/oliver/microsatellites/ratemicrosat/coverage/{sample}_avg_coverage.txt") #should change location of this file
        tissue_file = f"/faststorage/project/mutationalscanning/Workspaces/oliver/microsatellites/marmosetanalysis/{species}/{sample}_sorted_tissue_classification.tsv" #should change location of this file
>>>>>>> development:indel_mutation_call.py

        with open(index_assembly) as contigs:
            for line in contigs:
                contig_dict[line.split()[0]] = line.split()[1]
            for contig, length in contig_dict.items(): #removed the -d flag? # minscore chaged to 30 # maxperiod is 50 # added -l 10
                if (int(length) > 10**6):
                    #print(f"contig: {contig} length: {length}")
                    microsatfile = f"/faststorage/project/mutationalscanning/Workspaces/oliver/microsatellites/ratemicrosat/{species}/{sample}/{contig}_trs_tidy.tsv"
                    unresolved_bed = f"/faststorage/project/mutationalscanning/Workspaces/oliver/microsatellites/ratemicrosat/{species}/{sample}/{contig}_unresolved_haplotypes.bed"
                    gwf.target(#split fasta by contig so analysis is smooth
                        f"SmallIndels_{contig}_{sample}_diploid",
                        inputs=[assembly,alignemt,avg_coverage_file,tissue_file,microsatfile],
                        outputs=[f"output/small_indels/{species}/{sample}/{contig}_filters.tsv",
                                f"output/small_indels/{species}/{sample}/{contig}_small_indels.tsv"]
                        ) << f"""
                        mkdir -p output/
                        mkdir -p output/small_indels/{species}/
                        mkdir -p output/small_indels/{species}/{sample}/
                        awk 'OFS="\t" {{print $1,$2,$3}}' {microsatfile} | sort -k1,1 -k2,2n | bedtools merge -i - -i {unresolved_bed} > output/small_indels/{species}/{sample}/{contig}_filters.tsv
                        python {scripts}/small_indels.py output/small_indels/{species}/{sample}/{contig}_filters.tsv {contig} {alignemt} {assembly} {avg_coverage_file} {tissue_file} > output/small_indels/{species}/{sample}/{contig}_small_indels.tsv
                        """
                    gwf.target(#split fasta by contig so analysis is smooth
                        f"build_mutation_fastq_{contig}_{sample}_diploid",
                        inputs=[f"output/small_indels/{species}/{sample}/{contig}_small_indels.tsv"],
                        outputs=[f"output/mutation_fastq/{species}/{sample}/{contig}_remapping.fastq"]
                        ) << f"""
                        mkdir -p output/
                        mkdir -p output/mutation_fastq/{species}/
                        mkdir -p output/mutation_fastq/{species}/{sample}/
                        python {scripts}/build_mutation_fastq.py output/small_indels/{species}/{sample}/{contig}_small_indels.tsv {alignemt} {contig} "output/mutation_fastq/{species}/{sample}/{contig}_remapping.fastq"
                        """
                    # gwf.target(#split fasta by contig so analysis is smooth
                    #     f"LargeIndels_{contig}_{sample}_diploid",
                    #     inputs=[assembly,alignemt,avg_coverage_file,tissue_file],
                    #     outputs=[f"output/large_indels/{species}/{sample}/{contig}_large_indels.tsv"]
                    #     ) << f"""
                    #     mkdir -p output/
                    #     mkdir -p output/large_indels/{species}/
                    #     mkdir -p output/large_indels/{species}/{sample}/
                    #     python {scripts}/large_indels.py {contig} {alignemt} {assembly} {avg_coverage_file} {tissue_file} {microsatfile} > output/large_indels/{species}/{sample}/{contig}_large_indels.tsv
                    #     """
                    gwf.target(#split fasta by contig so analysis is smooth
                        f"opportunities_{contig}_{sample}_diploid",
                        inputs=[assembly,alignemt,avg_coverage_file,tissue_file,
                               f"output/small_indels/{species}/{sample}/{contig}_filters.tsv"],
                        outputs=[f"output/opportunities/{species}/{sample}/{contig}_opportunities.tsv"]
                        ) << f"""
                        mkdir -p output/
                        mkdir -p output/opportunities/{species}/
                        mkdir -p output/opportunities/{species}/{sample}/
                        python {scripts}/opportunities.py output/small_indels/{species}/{sample}/{contig}_filters.tsv {contig} {alignemt} {assembly} {avg_coverage_file} {tissue_file} > output/opportunities/{species}/{sample}/{contig}_opportunities.tsv
                        """
                    # gwf.target(#split fasta by contig so analysis is smooth
                    #     f"summary_SmallIndels_{contig}_{sample}_diploid",
                    #     inputs=[assembly,alignemt,avg_coverage_file,tissue_file,microsatfile],
                    #     outputs=[f"output/small_indels/{species}/{sample}/{contig}_small_indels_summary.tsv"]
                    #     ) << f"""
                    #     mkdir -p output/
                    #     mkdir -p output/small_indels/{species}/
                    #     mkdir -p output/small_indels/{species}/{sample}/
                    #     python {scripts}/small_indels_summary.py {microsatfile} {contig} {alignemt} {assembly} {avg_coverage_file} {tissue_file} > output/small_indels/{species}/{sample}/{contig}_small_indels_summary.tsv
                    #     """
                    
            # gwf.target(#split fasta by contig so analysis is smooth
            #     f"MergingIndels_{sample}_diploid",
            #     inputs=[f"output/large_indels/{species}/{sample}/{contig}_large_indels.tsv" for contig, length in contig_dict.items() if int(length)>10**6] +
            #             [f"output/small_indels/{species}/{sample}/{contig}_small_indels.tsv" for contig, length in contig_dict.items() if int(length)>10**6] +
            #             [f"output/small_indels/{species}/{sample}/{contig}_small_indels_summary.tsv" for contig, length in contig_dict.items() if int(length)>10**6] + 
            #             [f"output/opportunities/{species}/{sample}/{contig}_opportunities.tsv" for contig, length in contig_dict.items() if int(length)>10**6],
            #     outputs=[f"output/merged/{sample}_small_indels.tsv", 
            #             f"output/merged/{sample}_large_indels.tsv", 
            #             f"output/merged/{sample}_opportunities.tsv",
            #             f"output/merged/{sample}_small_indels_summary.tsv"],
            #     cores=8,
            #     memory="50g",
            #     walltime='1:00:00',
            #     account="mutationalscanning") << f"""
            #     mkdir -p output/
            #     mkdir -p output/merged/
            #     awk 'FNR==1 && NR==1 || FNR>1' output/small_indels/{species}/{sample}/*_small_indels.tsv > output/merged/{sample}_small_indels.tsv
            #     awk 'FNR==1 && NR==1 || FNR>1' output/small_indels/{species}/{sample}/*_small_indels_summary.tsv > output/merged/{sample}_small_indels_summary.tsv
            #     awk 'FNR==1 && NR==1 || FNR>1' output/large_indels/{species}/{sample}/*_large_indels.tsv > output/merged/{sample}_large_indels.tsv
            #     awk 'FNR==1 && NR==1 || FNR>1' output/opportunities/{species}/{sample}/*_opportunities.tsv > output/merged/{sample}_opportunities.tsv
            #     """
            gwf.target(#split fasta by contig so analysis is smooth
                f"haploid_mapping_{sample}_diploid",
                inputs=[f"output/mutation_fastq/{species}/{sample}/{contig}_remapping.fastq" for contig, length in contig_dict.items() if (int(length) > 10**6)],
                outputs=[f"output/haploid_mapping/{species}_{sample}_haploid_locations.tsv"],
                cores=8,
                memory="100g",
                walltime='24:00:00',
                account="mutationalscanning") << f"""
                mkdir -p output/haploid_mapping/
                cat output/mutation_fastq/{species}/{sample}/*_remapping.fastq > output/mutation_fastq/{species}/{sample}/{sample}_combined.fastq
                pbmm2 align {inpath_fa}/{species}/{sample}/haploid_assembly/{sample}_haploid.fa output/mutation_fastq/{species}/{sample}/{sample}_combined.fastq output/haploid_mapping/{species}_{sample}_remapped.bam --preset HIFI --sort -j 12 -J 4 -m 8G 
                samtools view output/haploid_mapping/{species}_{sample}_remapped.bam | awk -v OFS='\t' '{{print $1,$2,$3,$4,$5,$6}}' | sort -k1,1 > output/haploid_mapping/{species}_{sample}_haploid_locations.tsv
                """

            for contig, length in contig_dict.items(): #removed the -d flag? # minscore chaged to 30 # maxperiod is 50 # added -l 10
                if (int(length) > 10**6):
                    gwf.target(
                        f"haploid_filter_{contig}_{sample}_diploid",
                        inputs=[f"output/haploid_mapping/{species}_{sample}_haploid_locations.tsv",
                                f"output/small_indels/{species}/{sample}/{contig}_small_indels.tsv"],
                        outputs=[f"output/small_indels/{species}/{sample}/{contig}_small_indels_haploid_filtered.tsv"]
                        ) << f"""
                        mkdir -p output/
                        mkdir -p output/small_indels/{species}/
                        mkdir -p output/small_indels/{species}/{sample}/
                        python {scripts}/filtering_haploid_mappping.py output/small_indels/{species}/{sample}/{contig}_small_indels.tsv output/haploid_mapping/{species}_{sample}_haploid_locations.tsv > output/small_indels/{species}/{sample}/{contig}_small_indels_haploid_filtered_tmp.tsv
                        awk 'NR==1 {{print $0 "\t" "species" "\t" "sample"}}; NR>1 {{print $0 "\t" "{species}" "\t" "{sample}"}}' output/small_indels/{species}/{sample}/{contig}_small_indels_haploid_filtered_tmp.tsv > output/small_indels/{species}/{sample}/{contig}_small_indels_haploid_filtered.tsv
                        rm -f output/small_indels/{species}/{sample}/{contig}_small_indels_haploid_filtered_tmp.tsv
                        """
                        
                        #make some platau plots here based on kmerq and mapq of the haploid filtered reads
            gwf.target(
                f"qualtiy_plots_{sample}_diploid",
                inputs=[f"output/small_indels/{species}/{sample}/{contig}_small_indels_haploid_filtered.tsv" for contig, length in contig_dict.items() if (int(length) > 10**6)] +
                       [f"output/opportunities/{species}/{sample}/{contig}_opportunities.tsv" for contig, length in contig_dict.items() if (int(length) > 10**6)],
                outputs=[f"output/final/{sample}_haploid_combined.tsv",
<<<<<<< HEAD:workflow.py
                         f"output/final/{sample}_opportunities_combined.tsv"]
=======
                         f"output/final/{sample}_opportunities_combined.tsv",
                         f"output/plots/{sample}_indel_plots.pdf"],
                cores=8,
                memory="100g",
                walltime='24:00:00',
                account="mutationalscanning"
>>>>>>> development:indel_mutation_call.py
                ) << f"""
                mkdir -p output/
                mkdir -p output/final/
                mkdir -p output/plots/
                awk 'FNR==1 && NR==1 || FNR>1' output/small_indels/{species}/{sample}/*_small_indels_haploid_filtered.tsv > output/final/{sample}_haploid_combined.tsv
                awk 'FNR==1 && NR==1 || FNR>1' output/opportunities/{species}/{sample}/*_opportunities.tsv > output/final/{sample}_opportunities_combined.tsv
                """
#python {scripts}/quality_plots.py output/haploid_mapping/{species}_{sample}_haploid_locations.tsv output/plots/{species}_{sample}_{contig}_quality_plot.pdf

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
                mkdir -p output/
                mkdir -p output/final/
                python {scripts}/filtering_small_indels.py output/final/{sample}_haploid_combined.tsv -o output/final/{sample}_combined_filtered.tsv --haploid-check True --min-quality 93
                python {scripts}/filtering_opportunities.py output/final/{sample}_opportunities_combined.tsv -o output/final/{sample}_opportunities_combined_filtered.tsv  --min-kmerq 93
                """