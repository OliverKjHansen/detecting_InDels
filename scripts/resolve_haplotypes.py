import argparse
from utilities import resolve_haplotypes
import signal

# Prevent broken pipe errors when piping to head/less
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

# Usage: python resolve_haplotypes.py <contig> <start> <end> <bam> <assembly.fa> <avg_coverage.txt> \
#            --min-base-quality <q> --min-reads <n>
# Scans positions [start, end) of the contig (0-based, end exclusive) and writes
# the unresolved-haplotype regions found there to stdout as a headerless
# 3-column BED: contig, start, end. The options are the `precompute`
# parameters in config.yaml.

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unresolved-haplotype regions of one window of a contig.")
    parser.add_argument("contig")
    parser.add_argument("start", type=int)
    parser.add_argument("end", type=int)
    parser.add_argument("bam")
    parser.add_argument("fasta")
    parser.add_argument("avg_coverage")
    parser.add_argument("--min-base-quality", type=int, required=True)
    parser.add_argument("--min-reads", type=int, required=True)
    args = parser.parse_args()
    resolve_haplotypes(args.contig,
                       args.start,
                       args.end,
                       args.bam,
                       args.fasta,
                       args.avg_coverage,
                       min_base_quality=args.min_base_quality,
                       min_reads=args.min_reads)
