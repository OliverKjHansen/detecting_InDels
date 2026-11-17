"""
The measurable genome of a sample: the bases of the contigs analysed that
workflow 3 can test, the positions where its scan can count an opportunity.
That is the contigs' bases outside the mask (every repeat tract, the excluded
group and the unresolved-haplotype regions) whose depth is within the coverage
band. A rate per bp times the measurable genome is the InDels it holds.

    python measurable_genome.py band <avg_coverage.txt> --coverage-band LOW HIGH

prints mosdepth's --quantize thresholds, 0:MIN:MAX+1:, for the depths the scan
takes: from LOW to HIGH times the sample's mean coverage, both included, as
whole depths, and at least 1 (the scan has no pileup column where no read is).

    python measurable_genome.py count <contigs.bed> <bands.quantized.bed.gz> <out.tsv> \\
        --coverage <avg_coverage.txt> --coverage-band LOW HIGH --species S --sample X --masks <mask.bed> ...

reads the depth bands mosdepth writes with those thresholds, labelled below,
in and above (MOSDEPTH_Q0, Q1 and Q2), and writes one row: species, sample,
mean_coverage, min_depth, max_depth, contigs_bp (the contigs analysed),
masked_bp (their bases in the mask), below_band_bp and above_band_bp (their
unmasked bases at a depth below or above the band; a base no band covers is
below), and measurable_bp (the rest: unmasked, in the band). The four parts
add up to contigs_bp. Only the contigs of contigs.bed count, and the masks and
bands are clipped to them; the masks' intervals are joined where they overlap.

The depth must be counted as the scan's pileup counts it (pysam's nsegments,
stepper "nofilter", no mapping or base quality threshold): every read, with
any flag, and a read with a deletion at the position. mosdepth --fast-mode
--flag 0 --mapq 0 does.
"""

import argparse
import math
import os
import sys

import numpy
import pandas

LABELS = {"below", "in", "above"}


def depth_band(coverage_file, low, high):
    """(min_depth, max_depth): the whole depths from low to high times the mean coverage, at least 1."""
    with open(coverage_file) as handle:
        mean = float(handle.read().strip())
    return mean, max(1, math.ceil(mean * low)), math.floor(mean * high)


def read_bed(path, columns):
    """A headerless BED's first columns; an empty file has no rows."""
    if os.path.getsize(path) == 0:
        return pandas.DataFrame({column: pandas.Series(dtype=str if column in ("contig", "label") else "int64")
                                 for column in columns})
    return pandas.read_csv(path, sep="\t", header=None, usecols=range(len(columns)), names=columns,
                           dtype={"contig": str, "start": numpy.int64, "end": numpy.int64, "label": str})


def joined(starts, ends):
    """Intervals sorted, with those that overlap or touch joined, as (starts, ends)."""
    if len(starts) == 0:
        return starts, ends
    order = numpy.argsort(starts, kind="stable")
    starts, ends = starts[order], ends[order]
    # An interval starts a new one if it begins past every end before it.
    first = numpy.ones(len(starts), dtype=bool)
    first[1:] = starts[1:] > numpy.maximum.accumulate(ends)[:-1]
    return starts[first], numpy.maximum.reduceat(ends, numpy.flatnonzero(first))


def clipped(starts, ends, start, end):
    """Intervals cut to [start, end), empty ones left out."""
    starts, ends = numpy.clip(starts, start, end), numpy.clip(ends, start, end)
    keep = ends > starts
    return starts[keep], ends[keep]


def overlap(a, b):
    """The bases two sets of sorted, disjoint intervals share."""
    (a_starts, a_ends), (b_starts, b_ends) = a, b
    if len(a_starts) == 0 or len(b_starts) == 0:
        return 0
    before = numpy.concatenate(([0], numpy.cumsum(a_ends - a_starts)))

    def covered(x):
        # The bases of a before position x: the intervals ending by x whole,
        # and the part of the next one started before x.
        i = numpy.searchsorted(a_ends, x, side="right")
        partial = numpy.where(i < len(a_starts), numpy.maximum(0, x - a_starts[numpy.minimum(i, len(a_starts) - 1)]), 0)
        return before[i] + partial

    return int((covered(b_ends) - covered(b_starts)).sum())


def by_contig(table):
    return {contig: (rows["start"].to_numpy(), rows["end"].to_numpy()) for contig, rows in table.groupby("contig")}


def measurable(contigs, masks, bands):
    """{contigs_bp, masked_bp, below_band_bp, above_band_bp, measurable_bp} over the contigs."""
    unknown = set(bands["label"]) - LABELS
    if unknown:
        sys.exit(f"measurable_genome.py: the depth bands must be labelled below, in and above "
                 f"(MOSDEPTH_Q0, Q1 and Q2), not {', '.join(sorted(map(str, unknown)))}")
    masks = by_contig(masks)
    in_band = by_contig(bands[bands["label"] == "in"])
    above = by_contig(bands[bands["label"] == "above"])
    empty = (numpy.zeros(0, dtype=numpy.int64), numpy.zeros(0, dtype=numpy.int64))
    totals = dict.fromkeys(["contigs_bp", "masked_bp", "below_band_bp", "above_band_bp", "measurable_bp"], 0)
    for contig, start, end in contigs.itertuples(index=False):
        length = end - start
        mask = joined(*clipped(*masks.get(contig, empty), start, end))
        masked = int((mask[1] - mask[0]).sum())
        parts = {}
        for name, bands_of in (("measurable_bp", in_band), ("above_band_bp", above)):
            kept = joined(*clipped(*bands_of.get(contig, empty), start, end))
            parts[name] = int((kept[1] - kept[0]).sum()) - overlap(kept, mask)
        totals["contigs_bp"] += length
        totals["masked_bp"] += masked
        totals["measurable_bp"] += parts["measurable_bp"]
        totals["above_band_bp"] += parts["above_band_bp"]
        totals["below_band_bp"] += length - masked - parts["measurable_bp"] - parts["above_band_bp"]
    return totals


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    band_parser = commands.add_parser("band")
    band_parser.add_argument("coverage")
    band_parser.add_argument("--coverage-band", type=float, nargs=2, required=True, metavar=("LOW", "HIGH"))
    count_parser = commands.add_parser("count")
    count_parser.add_argument("contigs")
    count_parser.add_argument("bands")
    count_parser.add_argument("out")
    count_parser.add_argument("--coverage", required=True)
    count_parser.add_argument("--coverage-band", type=float, nargs=2, required=True, metavar=("LOW", "HIGH"))
    count_parser.add_argument("--species", required=True)
    count_parser.add_argument("--sample", required=True)
    count_parser.add_argument("--masks", nargs="*", default=[])
    args = parser.parse_args()

    mean, min_depth, max_depth = depth_band(args.coverage, *args.coverage_band)
    if args.command == "band":
        print(f"0:{min_depth}:{max_depth + 1}:")
        return

    contigs = read_bed(args.contigs, ["contig", "start", "end"])
    masks = pandas.concat([read_bed(path, ["contig", "start", "end"]) for path in args.masks]
                          or [read_bed(os.devnull, ["contig", "start", "end"])])
    bands = read_bed(args.bands, ["contig", "start", "end", "label"])
    row = {"species": args.species, "sample": args.sample, "mean_coverage": mean,
           "min_depth": min_depth, "max_depth": max_depth, **measurable(contigs, masks, bands)}
    pandas.DataFrame([row]).to_csv(args.out + ".tmp", sep="\t", index=False)
    os.replace(args.out + ".tmp", args.out)


if __name__ == "__main__":
    main()
