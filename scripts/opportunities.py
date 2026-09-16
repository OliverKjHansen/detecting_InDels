import pysam
import sys
import csv
from utilities import indel_opportunities
import signal

# Prevent broken pipe errors when piping to head/less
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

# Usage: opportunities.py <mask_bed> <contig> <bam> <fasta> <avg_coverage> <tissue> [<start> <end>]
#
# Counts InDel opportunities at positions OUTSIDE the mask BED, the denominator
# for the rate. The optional <start> <end> limit the count to positions
# [start, end) of the contig (0-based, end exclusive).

if __name__ == "__main__":
    if len(sys.argv) not in (7, 9):
        sys.exit(
            "Usage: opportunities.py "
            "<mask_bed> <contig> <bam> <fasta> <avg_coverage> <tissue> [<start> <end>]"
        )
    mask_bed = sys.argv[1]
    contig = sys.argv[2]
    bam_file = sys.argv[3]
    fasta_file = sys.argv[4]
    cov_file = sys.argv[5]
    tissue = sys.argv[6]
    window_start, window_end = (int(sys.argv[7]), int(sys.argv[8])) if len(sys.argv) == 9 else (None, None)
    indel_opportunities(mask_bed,
                        contig,
                        bam_file,
                        fasta_file,
                        cov_file,
                        tissue,
                        window_start,
                        window_end)
