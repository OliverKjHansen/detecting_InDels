import sys

from merge_regions import merge, overlaps, read_regions

# Usage: python drop_excluded_tracts.py <tracts.bed> <regions.bed> [<regions.bed> ...] > <measured.bed>
#
# The tracts workflow 4 measures: the rows of a repeat BED (header kept) that
# share no base with the regions, the homopolymer runs and the sub-windows the
# long-period TRF couldn't finish from repetitive_dna.py, and the
# unresolved-haplotype regions from precompute.py; not the long-period tracts,
# which a shorter period wins over (ticket 34). A tract that overlaps a region
# is dropped whole, so a tract is measured whole or not at all. Its bases stay
# in workflow 3's mask, so they are in neither workflow. Tracts that only touch
# a region do not overlap it.

if __name__ == "__main__":
    tracts, region_beds = sys.argv[1], sys.argv[2:]
    regions = {contig: merge(intervals) for contig, intervals in read_regions(region_beds).items()}

    with open(tracts) as handle:
        sys.stdout.write(handle.readline())
        for line in handle:
            contig, start, end = line.rstrip("\n").split("\t", 3)[:3]
            if not overlaps(regions.get(contig, []), int(start), int(end)):
                sys.stdout.write(line)
