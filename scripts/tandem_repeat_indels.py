import pysam
import sys
from utilities import tandem_repeat_indels
import signal

# Prevent broken pipe errors when piping to head/less
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

if __name__ == "__main__":
    if len(sys.argv) != 7:
        sys.exit(
            "Usage: tandem_repeat_indels.py "
            "<tr_bed> <contig> <bam> <fasta> <avg_coverage> <tissue>"
        )
    tr_bedfile           = sys.argv[1]
    contig               = sys.argv[2]
    bam_file             = sys.argv[3]
    fasta_file           = sys.argv[4]
    avg_coverage_file    = sys.argv[5]
    tissue               = sys.argv[6]
    tandem_repeat_indels(tr_bedfile, contig, bam_file, fasta_file,
                         avg_coverage_file, tissue)
