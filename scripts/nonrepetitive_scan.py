import os
import signal
import sys

from utilities import nonrepetitive_scan

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

# Calls 1-3 bp InDels and counts opportunities outside the mask in one pass over
# the window [start, end) of the BAM, and writes the reads carrying a call to a
# FASTQ.
#
# Usage: python nonrepetitive_scan.py <mask_bed> <contig> <bam> <fasta> <avg_coverage> <tissue> <start> <end> \
#            <calls.tsv> <opportunities.tsv> <reads.fastq>
#
# Every output is written to a .tmp file and renamed only when the whole window
# is done, so a job that dies never leaves a finished-looking output behind.

if __name__ == "__main__":
    if len(sys.argv) != 12:
        sys.exit(
            "Usage: nonrepetitive_scan.py <mask_bed> <contig> <bam> <fasta> <avg_coverage> <tissue> "
            "<start> <end> <calls.tsv> <opportunities.tsv> <reads.fastq>"
        )
    mask_bed, contig, bam_file, fasta_file, cov_file, tissue = sys.argv[1:7]
    window_start, window_end = int(sys.argv[7]), int(sys.argv[8])
    calls, opportunities, fastq = sys.argv[9:12]

    with open(f"{calls}.tmp", "w") as calls_out, \
         open(f"{opportunities}.tmp", "w") as opportunities_out, \
         open(f"{fastq}.tmp", "w") as fastq_out:
        nonrepetitive_scan(mask_bed, contig, bam_file, fasta_file, cov_file, tissue,
                           window_start, window_end,
                           calls_out=calls_out, opportunities_out=opportunities_out, fastq_out=fastq_out)

    for path in (calls, opportunities, fastq):
        os.replace(f"{path}.tmp", path)
