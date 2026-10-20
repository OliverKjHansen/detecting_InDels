"""
Read shapes that real data has and the main fixture lacks, each run through the
real scans: reverse-strand reads, base qualities below 93, reads that end
inside a tract or cross a window edge, soft clips, mapQ below 60, and
overlapping tracts. A bug fix builds its failing case from the same pieces.

Each test writes its own small BAM: the planted reads plus N_CLEAN clean reads
spanning the contig, so every position stays inside the coverage window.
"""

import pysam

import make_fixtures as mf
from scans import nonrepetitive_scan, rows_of, summed
from read_tags import decode_runs

REFERENCE = mf.build_reference()
N_CLEAN = 10


def dataset(tmp_path, reads, mask=()):
    """Write the reads plus N_CLEAN clean ones, with the given mask intervals."""
    clean = [mf.make_read(f"clean{i:02d}", REFERENCE) for i in range(N_CLEAN)]
    return mf.write_dataset(str(tmp_path), REFERENCE, reads + clean, mask=mask)


def calls_by_read(paths, window_start=None, window_end=None):
    calls, _, _ = nonrepetitive_scan(paths, window_start, window_end)
    return {call["readID"]: call for call in rows_of(calls)}


# -- reverse strand -----------------------------------------------------------


def test_calls_on_reverse_strand_reads_are_found_at_their_positions(tmp_path):
    reads = [mf.make_read("rev_del", REFERENCE, "del", 500, 1, reverse=True),
             mf.make_read("rev_ins", REFERENCE, "ins", 1003, 1, reverse=True)]
    calls = calls_by_read(dataset(tmp_path, reads))
    # `position` is the 1-based anchor: the base before the gap.
    assert (calls["rev_del"]["position"], calls["rev_del"]["diploid_direction"]) == ("500", "-")
    assert (calls["rev_ins"]["position"], calls["rev_ins"]["diploid_direction"]) == ("1004", "-")


def test_subread_tags_of_a_reverse_strand_read_run_in_sequencing_order(tmp_path):
    """
    The precondition of S1: pbmm2 leaves the per-base subread tags in the order
    the read was sequenced, so on a reverse-strand read the tags for aligned
    base i sit at index query_length - 1 - i.
    """
    read = mf.make_read("rev", REFERENCE, reverse=True, subreads={1000: (4, 6, 12)})
    paths = dataset(tmp_path, [read])
    with pysam.AlignmentFile(paths["bam"]) as bam:
        stored = next(r for r in bam.fetch(mf.CONTIG) if r.query_name == "rev")
    aligned = stored.get_aligned_pairs(matches_only=True)
    i = next(query for query, reference in aligned if reference == 1000)
    own = stored.query_length - 1 - i
    assert (stored.get_tag("sm")[own], stored.get_tag("sx")[own],
            decode_runs(stored.get_tag("sa"))[own]) == (4, 6, 12)
    assert (stored.get_tag("sm")[i], stored.get_tag("sx")[i],
            decode_runs(stored.get_tag("sa"))[i]) == (mf.SUBREADS, 0, mf.SUBREADS)


# -- base qualities below 93 --------------------------------------------------


def test_a_calls_kmerq_is_the_lowest_quality_beside_its_gap(tmp_path):
    reads = [
        # 1 bp deletion of base 500: the two bases on each side are 498, 499
        # and 501, 502. Base 400 is far away.
        mf.make_read("del_q30", REFERENCE, "del", 500, 1, qualities={502: 30, 400: 20}),
        # 1 bp insertion after base 1003: the two bases on each side are
        # 1002, 1003 and 1004, 1005, not the inserted base.
        mf.make_read("ins_q25", REFERENCE, "ins", 1003, 1, qualities={1005: 25, "inserted": 20}),
    ]
    calls = calls_by_read(dataset(tmp_path, reads))
    assert calls["del_q30"]["kmerQ"] == "30"
    assert calls["ins_q25"]["kmerQ"] == "25"


# -- soft clips ---------------------------------------------------------------


def test_soft_clipped_bases_count_in_read_coordinates_but_not_in_reference_ones(tmp_path):
    read = mf.make_read("clipped", REFERENCE, "del", 500, 1, soft_clip=(50, 30))
    call = calls_by_read(dataset(tmp_path, [read]))["clipped"]
    # Anchor 499 (0-based) is read base 50 + 499; `start` is 1-based.
    assert (call["position"], call["start"]) == ("500", "550")
    # The whole read: 6000 reference bases, one deleted, 80 clipped.
    assert call["read_length"] == "6079"


# -- mapQ below 60 ------------------------------------------------------------


def test_a_reads_opportunities_are_counted_under_its_mapq(tmp_path):
    read = mf.make_read("mapq5", REFERENCE, mapq=5)
    _, opportunities, _ = nonrepetitive_scan(dataset(tmp_path, [read]))
    by_mapq = {}
    for row in rows_of(opportunities):
        by_mapq[row["mapQ_bin"]] = by_mapq.get(row["mapQ_bin"], 0) + int(row["count"])
    # Every read is a clean copy of the contig, so each one counts the same
    # positions: the mapQ 5 read a tenth as many as the N_CLEAN at mapQ 60.
    assert set(by_mapq) == {"1-20", ">=60"}
    assert by_mapq["1-20"] * N_CLEAN == by_mapq[">=60"]


# -- reads that end early -----------------------------------------------------


def test_a_calls_distance_to_the_read_end_is_measured_on_its_own_read(tmp_path):
    read = mf.make_read("short", REFERENCE, "del", 500, 1, start=400, end=3000)
    call = calls_by_read(dataset(tmp_path, [read]))["short"]
    # The anchor, base 499, is 99 bases from the read's first aligned base.
    assert call["ToreadEnd"] == "99"


def test_reads_crossing_a_window_edge_add_up_to_the_whole_contig(tmp_path):
    reads = [mf.make_read("left", REFERENCE, "del", 1050, 1, start=0, end=1500),
             mf.make_read("right", REFERENCE, "ins", 1300, 1, start=1000, end=4000),
             mf.make_read("inside", REFERENCE, start=1090, end=1110)]
    paths = dataset(tmp_path, reads)
    whole = nonrepetitive_scan(paths)
    parts = [nonrepetitive_scan(paths, start, end) for start, end in [(0, 1100), (1100, 6000)]]

    assert [row for calls, _, _ in parts for row in rows_of(calls)] == rows_of(whole[0])
    assert {row["readID"] for row in rows_of(whole[0])} == {"left", "right"}
    assert summed(opportunities for _, opportunities, _ in parts) == summed([whole[1]])


# -- tandem repeat tracts -----------------------------------------------------
#
# Partial reads and overlapping tracts in workflow 4 are in test_tract_unit.py.
