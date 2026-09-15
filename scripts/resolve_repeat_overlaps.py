import sys
from intervaltree import IntervalTree

# Usage: python resolve_repeat_overlaps.py <repeats.bed> > <resolved.bed>
#
# Removes overlapping tandem repeat annotations so every base has at most one
# repeat state. When tracts overlap, the one with the smallest k is kept and the
# others are dropped whole. Among tracts with the same k, the longest is kept.
#
# Input and output are the unified 11-column repeat BED with a header line:
#   contig, start, end, canonical_motif, k, copies,
#   trf_purity, model_purity, model, species, sample
# Coordinates are 0-based half-open, so tracts that only touch do not overlap.

if __name__ == "__main__":
    bed = sys.argv[1]

    with open(bed) as handle:
        header = handle.readline()
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]

    # Smallest k first; within a k, longest tract first.
    rows.sort(key=lambda row: (int(row[4]), -(int(row[2]) - int(row[1]))))

    kept_trees = {}
    kept_rows = []
    for row in rows:
        contig, start, end = row[0], int(row[1]), int(row[2])
        tree = kept_trees.setdefault(contig, IntervalTree())
        if tree.overlaps(start, end):
            continue
        tree.addi(start, end)
        kept_rows.append(row)

    kept_rows.sort(key=lambda row: (row[0], int(row[1])))

    sys.stdout.write(header)
    for row in kept_rows:
        print("\t".join(row))
