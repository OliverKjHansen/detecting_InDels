import argparse
import json
import os
import signal

from text_files import open_text
from utilities import nonrepetitive_scan

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

# Calls 1-3 bp InDels and counts opportunities outside the mask in one pass over
# the window [start, end) of the BAM, and writes the reads carrying a call to a
# FASTQ.
#
# Usage: python nonrepetitive_scan.py <mask_bed> <contig> <bam> <fasta> <avg_coverage> <tissue> <start> <end> \
#            <calls.tsv> <opportunities.tsv> <reads.fastq> \
#            --coverage-band <low> <high> --indel-sizes <smallest> <largest> --flank-width <bp> \
#            --key-parameters <json>
#
# The first three options are the `scan` parameters in config.yaml, and
# --key-parameters its `strata` block, as JSON.
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
    parser.add_argument("--flank-width", type=int, required=True, metavar="BP")
    parser.add_argument("--read-classification", default="NA", metavar="TSV",
                        help="read_classification.sh's table of the sample's reads (default: none)")
    parser.add_argument("--key-parameters", type=json.loads, required=True, metavar="JSON")
    args = parser.parse_args()
    calls, opportunities, fastq = args.calls, args.opportunities, args.fastq

    # An output whose name ends in .gz is written gzipped.
    with open_text(f"{calls}.tmp", "w", gzipped=calls.endswith(".gz")) as calls_out, \
         open_text(f"{opportunities}.tmp", "w", gzipped=opportunities.endswith(".gz")) as opportunities_out, \
         open_text(f"{fastq}.tmp", "w", gzipped=fastq.endswith(".gz")) as fastq_out:
        nonrepetitive_scan(args.mask_bed, args.contig, args.bam, args.fasta, args.avg_coverage, args.tissue,
                           args.start, args.end,
                           calls_out=calls_out, opportunities_out=opportunities_out, fastq_out=fastq_out,
                           coverage_band=args.coverage_band, indel_sizes=args.indel_sizes,
                           key_parameters=args.key_parameters, flank_width=args.flank_width,
                           classification_file=args.read_classification)

    for path in (calls, opportunities, fastq):
        os.replace(f"{path}.tmp", path)
