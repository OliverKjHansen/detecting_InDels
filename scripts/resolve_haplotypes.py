import sys
from utilities import resolve_haplotypes
import signal

# Prevent broken pipe errors when piping to head/less
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

# Usage: python resolve_haplotypes.py <contig> <bam> <assembly.fa> <avg_coverage.txt>
# Writes the contig's unresolved-haplotype regions to stdout as a headerless
# 3-column BED: contig, start, end.

if __name__ == "__main__":
    contig = sys.argv[1]
    bam_file = sys.argv[2]
    fasta_file = sys.argv[3]
    cov_file = sys.argv[4]
    resolve_haplotypes(contig,
                       bam_file,
                       fasta_file,
                       cov_file)
