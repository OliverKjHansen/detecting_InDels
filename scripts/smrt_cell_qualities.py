"""
Each SMRT Cell's highest base quality, from the reads of a sample's BAM, for
the Q40 rule of final filtering (Q12, Q17): a cell whose reads never go above
Q40 had its base qualities capped there.

    python smrt_cell_qualities.py <alignment.bam> <smrt_cells.tsv> --reads N

Looks at the first N alignments of the BAM, which in a position-sorted BAM of
HiFi reads come from every cell. Writes one row per cell: smrt_cell, the
reads looked at, and max_quality.
"""
import argparse
from collections import defaultdict
from itertools import islice

import pysam

from strata import smrt_cell


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bam")
    parser.add_argument("out")
    parser.add_argument("--reads", type=int, required=True, metavar="N")
    args = parser.parse_args()

    reads, best = defaultdict(int), defaultdict(int)
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for alignment in islice(bam.fetch(until_eof=True), args.reads):
            qualities = alignment.query_qualities
            if qualities is None:
                continue
            cell = smrt_cell(alignment.query_name)
            reads[cell] += 1
            best[cell] = max(best[cell], max(qualities))
    with open(args.out, "w") as out:
        print("smrt_cell", "reads", "max_quality", sep="\t", file=out)
        for cell in sorted(reads):
            print(cell, reads[cell], best[cell], sep="\t", file=out)


if __name__ == "__main__":
    main()
