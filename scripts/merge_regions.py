import sys
from bisect import bisect_right
from collections import defaultdict

# Usage: python merge_regions.py <bed> [<bed> ...] > <merged.bed>
#
# Writes the bases covered by any interval of the BEDs as a headerless BED of
# merged intervals, contig, start, end, sorted by contig and start, as
# `sort -k1,1 -k2,2n | bedtools merge` would. Intervals that overlap or touch
# become one. Only the first three columns are read, and a line whose start or
# end is not a number, such as the header of a repeat BED, is skipped, so the
# repeat BEDs, the unresolved-haplotype regions and the excluded regions can
# be given as they are.


def read_regions(paths):
    """{contig: [(start, end)]} from the BEDs, header lines and empty intervals left out."""
    regions = defaultdict(list)
    for path in paths:
        with open(path) as handle:
            for line in handle:
                fields = line.split("\t", 3)
                if len(fields) < 3 or not (fields[1].isdigit() and fields[2].strip().isdigit()):
                    continue
                start, end = int(fields[1]), int(fields[2])
                if end > start:
                    regions[fields[0]].append((start, end))
    return regions


def merge(intervals):
    """The intervals sorted, with those that overlap or touch joined."""
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(interval) for interval in merged]


def overlaps(merged, start, end):
    """True if [start, end) shares a base with one of the merged intervals."""
    i = bisect_right(merged, (start, float("inf")))
    if i > 0 and merged[i - 1][1] > start:   # the interval starting at or before `start`
        return True
    return i < len(merged) and merged[i][0] < end   # the next one starts inside


if __name__ == "__main__":
    regions = read_regions(sys.argv[1:])
    for contig in sorted(regions):
        for start, end in merge(regions[contig]):
            print(contig, start, end, sep="\t")
