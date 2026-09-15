import sys
from bisect import bisect_right
from collections import defaultdict

# Usage: python resolve_repeat_overlaps.py <repeats.bed> > <resolved.bed>
#
# Removes overlapping tandem repeat annotations so every base has at most one
# repeat state. When tracts overlap, the one with the smallest k is kept and the
# others are dropped whole. Among tracts with the same k, the longest is kept
# (ties: the one listed first).
#
# Input and output are the unified 11-column repeat BED with a header line:
#   contig, start, end, canonical_motif, k, copies,
#   trf_purity, model_purity, model, species, sample
# Coordinates are 0-based half-open, so tracts that only touch do not overlap.
#
# How it stays fast on millions of tracts: the periods are handled one at a
# time, smallest first. The tracts kept so far never overlap each other, so kept
# starts and ends are both sorted, and whether a new tract hits one of them is
# a binary search. Within one period, only tracts that overlap each other (a
# cluster) have to be compared, and clusters are small.


def overlaps_kept(starts, ends, start, end):
    """True if [start, end) overlaps a kept tract (starts/ends sorted, disjoint)."""
    i = bisect_right(starts, start)
    if i > 0 and ends[i - 1] > start:   # kept tract starting at or before `start`
        return True
    return i < len(starts) and starts[i] < end   # next kept tract starts inside


def longest_first(tracts):
    """
    Resolve overlaps among tracts of one period: longest kept, ties by input
    order. tracts are (start, end, index); returns the kept ones.
    """
    kept = []
    tracts.sort()
    cluster, cluster_end = [], -1
    for tract in tracts + [(float("inf"), float("inf"), -1)]:
        if tract[0] >= cluster_end and cluster:
            if len(cluster) == 1:
                kept.extend(cluster)
            else:
                chosen = []
                for c in sorted(cluster, key=lambda t: (-(t[1] - t[0]), t[2])):
                    if all(c[1] <= o[0] or o[1] <= c[0] for o in chosen):
                        chosen.append(c)
                kept.extend(chosen)
            cluster, cluster_end = [], -1
        cluster.append(tract)
        cluster_end = max(cluster_end, tract[1])
    return kept


def resolve_contig(tracts_by_k):
    """tracts_by_k: {k: [(start, end, index)]} for one contig -> kept tracts."""
    kept = []
    for k in sorted(tracts_by_k):
        starts = [t[0] for t in kept]
        ends = [t[1] for t in kept]
        free = [t for t in tracts_by_k[k] if not overlaps_kept(starts, ends, t[0], t[1])]
        kept = sorted(kept + longest_first(free))
    return kept


if __name__ == "__main__":
    bed = sys.argv[1]

    with open(bed) as handle:
        header = handle.readline()
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]

    tracts = defaultdict(lambda: defaultdict(list))   # contig -> k -> tracts
    for index, row in enumerate(rows):
        start, end = int(row[1]), int(row[2])
        if end > start:
            tracts[row[0]][int(row[4])].append((start, end, index))

    kept_rows = []
    for contig in sorted(tracts):
        kept_rows.extend(rows[t[2]] for t in resolve_contig(tracts[contig]))

    sys.stdout.write(header)
    for row in kept_rows:
        print("\t".join(row))
