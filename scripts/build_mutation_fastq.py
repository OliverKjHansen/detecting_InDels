#!/usr/bin/env python3

import sys
import pysam
import csv
from collections import defaultdict


def parse_mutation_file(mutation_tsv):
    """
    Parse the small indels TSV file and extract a mapping from
    genomic position (1-based) to supporting read name.

    Assumes that:
    - position is stored in a column called 'pos' or similar
    - read name is stored in a column called 'read_name'

    Adjust column names here if necessary.
    """
    mutation_dict = defaultdict(list)
    read_check = list() # to avoid duplicates of reads

    with open(mutation_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pos = int(row["position"])
            read_name = row["readID"]

            # Remove suffixes such as '+', if present
            if "+" in read_name:
                read_name = read_name.splist("+")[0]
            
            if read_name in read_check:
                continue
            
            mutation_dict[pos].append(read_name)
            read_check.append(read_name)
    #print(f"Total unique mutation-supporting reads: {len(read_check)}")
    return mutation_dict


def extract_mutation_reads(bam_path, contig, mutation_dict):
    """
    Given a BAM file, a contig name, and a mapping from position
    to mutation-supporting read names, extract FASTQ entries
    for those reads.

    Returns:
        fastq_dict: dict mapping read_name -> FASTQ entry (list of 4 strings)
    """
    fastq_dict = {}
    count = 0

    #reference = pysam.FastaFile(fa_in) dont need fasta for this
    bam = pysam.AlignmentFile(bam_path, "rb")

    for pileupcolumn in bam.pileup(contig=contig, #pileup 
                                       start=None,
                                       stop=None,
                                       region=None,
                                       end=None,
                                       truncate=True,
                                       stepper="nofilter",
                                       min_base_quality=0, #quite light
                                       min_mapping_quality=0 #conservative
                                       ):
        pos_1based = pileupcolumn.pos + 1

        if pos_1based not in mutation_dict:
            continue
        #print("working?1", pos_1based)

        target_reads = mutation_dict[pos_1based]
        #print(pos_1based,target_read_name)
        for pileupread in pileupcolumn.pileups:
            #print("pileupread:", pileupread.alignment.query_name, pos_1based, target_read_name)
            if pileupread.alignment.query_name in target_reads:
                count += 1
                #print("working?2", target_read_name)
                qual = ''.join([chr(x+33) for x in pileupread.alignment.get_forward_qualities()])
                fastq_dict[pileupread.alignment.query_name] = [
                f"@{pileupread.alignment.query_name}",
                pileupread.alignment.get_forward_sequence(),#sequence
                "+",
                qual
            ]
    #print(f"Total extracted mutation-supporting reads from BAM: {count}")
    bam.close()
    return fastq_dict


def write_fastq(fastq_dict, outfile):
    """
    Write FASTQ entries to file.
    """
    with open(outfile, "w") as out:
        for read_name in fastq_dict:
            out.write("\n".join(fastq_dict[read_name]) + "\n")


def main():
    if len(sys.argv) != 5:
        sys.exit(
            "Usage: build_mutation_fastq.py "
            "<small_indels.tsv> <contig> <bam> <output.fastq>"
        )

    mutation_tsv = sys.argv[1]
    bam_path = sys.argv[2]
    contig = sys.argv[3]
    outfile = sys.argv[4]

    mutation_dict = parse_mutation_file(mutation_tsv)
    #print(mutation_dict)
    fastq_dict = extract_mutation_reads(bam_path, contig, mutation_dict)
    write_fastq(fastq_dict, outfile)


if __name__ == "__main__":
    main()
#~/mutationalscanning/DerivedData/bam/HiFi/marmoset/lambiek/kinetics/lambiek_kinetics_diploid.bam 