import pysam
import sys
from utilities import tandem_repeat_opportunities
import signal

# Prevent broken pipe errors when piping to head/less
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

# The optional <start> <end> limit the count to positions [start, end) of the
# contig (0-based, end exclusive); without them the whole contig is counted.

if __name__ == "__main__":
    if len(sys.argv) not in (7, 9):
        sys.exit(
            "Usage: tandem_repeat_opportunities.py "
            "<tr_bed> <contig> <bam> <fasta> <avg_coverage> <tissue> [<start> <end>]"
        )
    tr_bedfile           = sys.argv[1]
    contig               = sys.argv[2]
    bam_file             = sys.argv[3]
    fasta_file           = sys.argv[4]
    avg_coverage_file    = sys.argv[5]
    tissue               = sys.argv[6]
    window_start, window_end = (int(sys.argv[7]), int(sys.argv[8])) if len(sys.argv) == 9 else (None, None)
    tandem_repeat_opportunities(tr_bedfile, contig, bam_file, fasta_file,
                                avg_coverage_file, tissue, window_start, window_end)
