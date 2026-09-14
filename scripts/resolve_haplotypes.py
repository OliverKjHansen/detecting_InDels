import sys
from utilities import resolve_haplotypes
import signal

# Prevent broken pipe errors when piping to head/less
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

# Usage: python resolve_haplotypes.py <contig> <start> <end> <bam> <assembly.fa> <avg_coverage.txt>
# Scans positions [start, end) of the contig (0-based, end exclusive) and writes
# the unresolved-haplotype regions found there to stdout as a headerless
# 3-column BED: contig, start, end.

if __name__ == "__main__":
    contig = sys.argv[1]
    window_start = int(sys.argv[2])
    window_end = int(sys.argv[3])
    bam_file = sys.argv[4]
    fasta_file = sys.argv[5]
    cov_file = sys.argv[6]
    resolve_haplotypes(contig,
                       window_start,
                       window_end,
                       bam_file,
                       fasta_file,
                       cov_file)
