import pysam
import sys
import csv
from utilities import small_indels
import signal

# Prevent broken pipe errors when piping to head/less
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

# Usage: small_indels.py <mask_bed> <contig> <bam> <fasta> <avg_coverage> <tissue> [<start> <end>]
#
# Calls 1-3 bp InDels at positions OUTSIDE the mask BED. The optional
# <start> <end> limit the scan to positions [start, end) of the contig
# (0-based, end exclusive); without them the whole contig is scanned.

if __name__ == "__main__":
    if len(sys.argv) not in (7, 9):
        sys.exit(
            "Usage: small_indels.py "
            "<mask_bed> <contig> <bam> <fasta> <avg_coverage> <tissue> [<start> <end>]"
        )
    mask_bed = sys.argv[1]
    contig = sys.argv[2]
    bam_file = sys.argv[3]
    fasta_file = sys.argv[4]
    cov_file = sys.argv[5]
    tissue = sys.argv[6]
    window_start, window_end = (int(sys.argv[7]), int(sys.argv[8])) if len(sys.argv) == 9 else (None, None)
    small_indels(mask_bed,
                 contig,
                 bam_file,
                 fasta_file,
                 cov_file,
                 tissue,
                 window_start,
                 window_end)
