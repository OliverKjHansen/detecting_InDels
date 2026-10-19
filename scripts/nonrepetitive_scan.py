import argparse
import os
import signal

from utilities import nonrepetitive_scan

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

# Calls 1-3 bp InDels and counts opportunities outside the mask in one pass over
# the window [start, end) of the BAM, and writes the reads carrying a call to a
# FASTQ.
#
# Usage: python nonrepetitive_scan.py <mask_bed> <contig> <bam> <fasta> <avg_coverage> <tissue> <start> <end> \
#            <calls.tsv> <opportunities.tsv> <reads.fastq> \
#            --coverage-band <low> <high> --indel-sizes <smallest> <largest>
#
# The two options are the `scan` parameters in config.yaml.
#
# Every output is written to a .tmp file and renamed only when the whole window
# is done, so a job that dies never leaves a finished-looking output behind.

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Call InDels and count opportunities outside the mask in one pass over a BAM window.")
    for name in ("mask_bed", "contig", "bam", "fasta", "avg_coverage", "tissue"):
        parser.add_argument(name)
    parser.add_argument("start", type=int)
    parser.add_argument("end", type=int)
    for name in ("calls", "opportunities", "fastq"):
        parser.add_argument(name)
    parser.add_argument("--coverage-band", type=float, nargs=2, required=True, metavar=("LOW", "HIGH"))
    parser.add_argument("--indel-sizes", type=int, nargs=2, required=True, metavar=("SMALLEST", "LARGEST"))
    args = parser.parse_args()
    calls, opportunities, fastq = args.calls, args.opportunities, args.fastq

    with open(f"{calls}.tmp", "w") as calls_out, \
         open(f"{opportunities}.tmp", "w") as opportunities_out, \
         open(f"{fastq}.tmp", "w") as fastq_out:
        nonrepetitive_scan(args.mask_bed, args.contig, args.bam, args.fasta, args.avg_coverage, args.tissue,
                           args.start, args.end,
                           calls_out=calls_out, opportunities_out=opportunities_out, fastq_out=fastq_out,
                           coverage_band=args.coverage_band, indel_sizes=args.indel_sizes)

    for path in (calls, opportunities, fastq):
        os.replace(f"{path}.tmp", path)
