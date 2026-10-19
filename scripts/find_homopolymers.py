import re
import sys

import pysam

# Usage: python find_homopolymers.py <contig.fa> <min_length> > <homopolymers.bed>
#
# Writes every run of one base at least min_length bp long in the FASTA, as a
# headerless BED: contig, start, end, base. Coordinates are 0-based half-open.
# Case is ignored, so soft-masked bases count; N is not a base, so runs of N
# (assembly gaps) are not homopolymers.


def homopolymer_runs(sequence, min_length):
    """(start, end, base) of every run of one of A, C, G and T at least min_length long."""
    run = re.compile(r"([ACGT])\1{%d,}" % (min_length - 1))
    for match in run.finditer(sequence.upper()):
        yield match.start(), match.end(), match.group(1)


if __name__ == "__main__":
    fasta, min_length = sys.argv[1], int(sys.argv[2])
    with pysam.FastxFile(fasta) as records:
        for record in records:
            for start, end, base in homopolymer_runs(record.sequence, min_length):
                print(record.name, start, end, base, sep="\t")
