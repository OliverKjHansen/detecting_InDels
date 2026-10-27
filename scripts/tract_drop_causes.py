"""
Why workflow 4's step 0 drops tracts under each rule ticket 34 weighed, for
one sample; option 2 was chosen (decisions Q16). Run from the repository root, where
config.yaml is, after workflows 1 and 2:

    python scripts/tract_drop_causes.py <species> <sample> [--edges 1,2,5,10] > tract_drop_causes.tsv

Before ticket 34, step 0 dropped a tract whole when it shared a base with an
unresolved-haplotype region or the excluded group: homopolymer runs,
long-period tracts, and the sub-windows the long-period TRF couldn't finish in
time (most likely satellite arrays). Per model, the table counts the tracts
and their bases:

  all                   every tract of the model
  dropped_unresolved    sharing a base with an unresolved region
  dropped_homopolymer   else, sharing a base with a homopolymer run
  dropped_timed_out     else, sharing a base with a sub-window TRF couldn't finish
  dropped_long_period   else, sharing a base with a long-period tract
  option_1              the rest: what step 0 measured before ticket 34
  option_2              long-period tracts no longer drop a tract; the other
                        three causes still do: what step 0 measures now
  option_3_edge_N       as option 1, but kept when each excluded region it meets
                        reaches in from one of its ends, holding N of its
                        bases or fewer at both ends together
  option_2_3_edge_N     as option 2, but kept when the homopolymer runs and
                        unfinished sub-windows it meets do so

Under every rule a tract sharing a base with an unresolved region is dropped,
and so is a tract an excluded region covers from end to end.
"""
import argparse
import os
import sys
from bisect import bisect_right
from collections import Counter

from merge_regions import merge, overlaps, read_regions

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "workflows"))
from pipeline import TIMED_OUT, Pipeline  # noqa: E402


def tracts(path):
    """(start, end) of each tract in a repeat BED, header skipped."""
    with open(path) as handle:
        handle.readline()
        for line in handle:
            start, end = line.split("\t", 3)[1:3]
            yield int(start), int(end)


def timed_out(paths):
    """(start, end) of the sub-windows marked TIMED_OUT in the long-period window BEDs."""
    regions = []
    for path in paths:
        with open(path) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) > 3 and fields[3] == TIMED_OUT:
                    regions.append((int(fields[1]), int(fields[2])))
    return regions


def bases_at_edges(merged, start, end):
    """
    The bases [start, end) shares with the merged intervals when each of those
    it meets reaches in from one of its ends, or None when one lies inside it,
    away from both ends, or covers it from end to end.
    """
    i = max(bisect_right(merged, (start, float("inf"))) - 1, 0)   # the interval starting at or before `start`
    shared = 0
    while i < len(merged) and merged[i][0] < end:
        first, last = merged[i]
        i += 1
        if last <= start:
            continue
        if (first <= start) == (last >= end):
            return None
        shared += min(last, end) - max(first, start)
    return shared


def measures(regions, start, end, edges):
    """The measures the tract [start, end) counts towards, by the contig's merged regions."""
    if overlaps(regions["unresolved"], start, end):
        return ["all", "dropped_unresolved"]
    if overlaps(regions["homopolymers"], start, end):
        found = ["all", "dropped_homopolymer"]
    elif overlaps(regions["timed_out"], start, end):
        found = ["all", "dropped_timed_out"]
    elif overlaps(regions["excluded"], start, end):
        found = ["all", "dropped_long_period", "option_2"]
    else:
        found = ["all", "option_1", "option_2"]
    excluded = bases_at_edges(regions["excluded"], start, end)
    option_2_excluded = bases_at_edges(regions["option_2_excluded"], start, end)
    found += [f"option_3_edge_{edge}" for edge in edges if excluded is not None and excluded <= edge]
    found += [f"option_2_3_edge_{edge}" for edge in edges
              if option_2_excluded is not None and option_2_excluded <= edge]
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("species")
    parser.add_argument("sample")
    parser.add_argument("--edges", default="1,2,5,10",
                        help="the N of option_3_edge_N and option_2_3_edge_N, comma-separated")
    args = parser.parse_args()
    edges = [int(edge) for edge in args.edges.split(",")]

    pipeline = Pipeline()
    models = pipeline.config["repeats"]["models"]
    tract_counts, bp_counts = Counter(), Counter()
    for contig, length in pipeline.long_contigs(args.species, args.sample).items():
        path = lambda name, **fields: pipeline.path(name, args.species, args.sample, contig=contig, **fields)  # noqa: E731
        regions = {name: merge(read_regions([path(name)])[contig])
                   for name in ("unresolved", "homopolymers", "excluded")}
        regions["timed_out"] = merge(timed_out(pipeline.long_period_window_beds(args.species, args.sample,
                                                                                contig, length)))
        regions["option_2_excluded"] = merge(regions["homopolymers"] + regions["timed_out"])
        for model in models:
            for start, end in tracts(path("tandem_repeats", model=model)):
                for measure in measures(regions, start, end, edges):
                    tract_counts[model, measure] += 1
                    bp_counts[model, measure] += end - start

    print("model", "measure", "tracts", "bp", "species", "sample", sep="\t")
    for model in models:
        for measure in ["all", "dropped_unresolved", "dropped_homopolymer", "dropped_timed_out",
                        "dropped_long_period", "option_1", "option_2",
                        *(f"option_3_edge_{edge}" for edge in edges),
                        *(f"option_2_3_edge_{edge}" for edge in edges)]:
            print(model, measure, tract_counts[model, measure], bp_counts[model, measure],
                  args.species, args.sample, sep="\t")


if __name__ == "__main__":
    main()
