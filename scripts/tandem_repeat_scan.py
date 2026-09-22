import argparse
import os
import signal

from utilities import tandem_repeat_scan

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

# Calls InDels and counts opportunities inside the repeat tracts of every model
# given, in one pass over the window [start, end) of the BAM, and writes the
# reads carrying a call in any model to one FASTQ.
#
# Usage: python tandem_repeat_scan.py <contig> <bam> <fasta> <avg_coverage> <tissue> <start> <end> <reads.fastq> \
#            --model M0 <M0.bed> <M0_calls.tsv> <M0_opportunities.tsv> \
#            --model M1 ...
#
# Every output is written to a .tmp file and renamed only when the whole window
# is done, so a job that dies never leaves a finished-looking output behind.


def main():
    parser = argparse.ArgumentParser(
        description="Call InDels and count opportunities for several repeat models in one pass over a BAM window.")
    parser.add_argument("contig")
    parser.add_argument("bam")
    parser.add_argument("fasta")
    parser.add_argument("avg_coverage")
    parser.add_argument("tissue")
    parser.add_argument("start", type=int)
    parser.add_argument("end", type=int)
    parser.add_argument("fastq", help="reads carrying a call in any model")
    parser.add_argument("--model", nargs=4, action="append", required=True,
                        metavar=("NAME", "BED", "CALLS", "OPPORTUNITIES"))
    args = parser.parse_args()

    outputs = [args.fastq] + [path for _, _, calls, opportunities in args.model
                              for path in (calls, opportunities)]
    handles = {path: open(f"{path}.tmp", "w") for path in outputs}

    tandem_repeat_scan(
        {name: bed for name, bed, _, _ in args.model},
        args.contig, args.bam, args.fasta, args.avg_coverage, args.tissue,
        args.start, args.end,
        calls_out={name: handles[calls] for name, _, calls, _ in args.model},
        opportunities_out={name: handles[opportunities] for name, _, _, opportunities in args.model},
        fastq_out=handles[args.fastq],
    )

    for path, handle in handles.items():
        handle.close()
        os.replace(f"{path}.tmp", path)


if __name__ == "__main__":
    main()
