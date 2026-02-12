import pysam
import sys
import csv
from utilities import indel_opportunities
import signal

# Prevent broken pipe errors when piping to head/less
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

# bam_file = "/faststorage/project/mutationalscanning/DerivedData/bam/HiFi/marmoset/marmoset/kinetics/marmoset_kinetics_diploid.bam"
# fasta_file = "/faststorage/project/mutationalscanning/DerivedData/fasta/HiFi/marmoset/marmoset/diploid_assembly/marmoset_diploid.fa"
# contig = "h1tg000170l"
# cov_file = "/faststorage/project/mutationalscanning/Workspaces/oliver/microsatellites/ratemicrosat/coverage/marmoset_avg_coverage.txt"
# #microsatfile = "microsat_short.txt"
# microsatfile = "/faststorage/project/mutationalscanning/Workspaces/oliver/microsatellites/ratemicrosat/marmoset/marmoset/h1tg000170l_trs_filtered_sorted.tsv"
# tissue = "/faststorage/project/mutationalscanning/Workspaces/oliver/microsatellites/marmosetanalysis/marmoset/marmoset_sorted_tissue_classification.tsv"

if __name__ == "__main__":
    microsatfile = sys.argv[1]
    contig = sys.argv[2]
    bam_file = sys.argv[3]
    fasta_file = sys.argv[4]
    cov_file = sys.argv[5]
    tissue = sys.argv[6]
    indel_opportunities(microsatfile,
                 contig,
                 bam_file,
                 fasta_file,
                 cov_file,
                 tissue)