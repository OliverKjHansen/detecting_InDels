import argparse
import json
import os
import signal

from utilities import tandem_repeat_scan

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

# Calls InDels and counts opportunities inside the repeat tracts of every model
# given, in one pass over the window [start, end) of the BAM, and writes the
# reads carrying a call in any model to one FASTQ.
#
# Usage: python tandem_repeat_scan.py <contig> <bam> <fasta> <avg_coverage> <tissue> <start> <end> <reads.fastq> \
#            --coverage-band <low> <high> --flank-width <bp> --key-parameters <json> \
#            --model M0 <M0.bed> <M0_calls.tsv> <M0_opportunities.tsv> \
#            --model M1 ...
#
# --coverage-band and --flank-width are `scan` parameters in config.yaml, and
# --key-parameters its `strata` block, as JSON.
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
    parser.add_argument("--coverage-band", type=float, nargs=2, required=True, metavar=("LOW", "HIGH"))
    parser.add_argument("--flank-width", type=int, required=True, metavar="BP")
    parser.add_argument("--read-classification", default="NA", metavar="TSV",
                        help="read_classification.sh's table of the sample's reads (default: none)")
    parser.add_argument("--key-parameters", type=json.loads, required=True, metavar="JSON")
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
        coverage_band=args.coverage_band,
        key_parameters=args.key_parameters,
        flank_width=args.flank_width,
        classification_file=args.read_classification,
    )

    for path, handle in handles.items():
        handle.close()
        os.replace(f"{path}.tmp", path)


if __name__ == "__main__":
    main()
