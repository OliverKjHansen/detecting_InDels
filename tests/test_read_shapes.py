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
from scans import nonrepetitive_scan, rows_of, tandem_repeat_scan
from utilities import subread_array

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
            subread_array(stored.get_tag("sa"), own)) == (4, 6, 12)
    assert (stored.get_tag("sm")[i], stored.get_tag("sx")[i],
            subread_array(stored.get_tag("sa"), i)) == (mf.SUBREADS, 0, mf.SUBREADS)


# -- base qualities below 93 --------------------------------------------------


def test_a_calls_kmerq_is_the_lowest_quality_beside_its_gap(tmp_path):
    reads = [
        # 1 bp deletion of base 500: the two bases on each side are 498, 499
        # and 501, 502. Base 400 is far away.
        mf.make_read("del_q30", REFERENCE, "del", 500, 1, qualities={502: 30, 400: 20}),
        # 1 bp insertion after base 1003: the window is 1003, the inserted
        # base and 1004.
        mf.make_read("ins_q25", REFERENCE, "ins", 1003, 1, qualities={"inserted": 25}),
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
        by_mapq[row["mapQ"]] = by_mapq.get(row["mapQ"], 0) + int(row["count"])
    # Every read is a clean copy of the contig, so each one counts the same
    # positions: the mapQ 5 read a tenth as many as the N_CLEAN at mapQ 60.
    assert set(by_mapq) == {"5", "60"}
    assert by_mapq["5"] * N_CLEAN == by_mapq["60"]


# -- reads that end early -----------------------------------------------------


def test_a_calls_distance_to_the_read_end_is_measured_on_its_own_read(tmp_path):
    read = mf.make_read("short", REFERENCE, "del", 500, 1, start=400, end=3000)
    call = calls_by_read(dataset(tmp_path, [read]))["short"]
    # The anchor, base 499, is 99 bases from the read's first aligned base.
    assert call["ToreadEnd"] == "99"


def summed_opportunities(texts, key_columns):
    totals = {}
    for text in texts:
        for row in rows_of(text):
            key = tuple(row[column] for column in key_columns)
            totals[key] = totals.get(key, 0) + int(row["count"])
    return {key: count for key, count in totals.items() if count}


def test_reads_crossing_a_window_edge_add_up_to_the_whole_contig(tmp_path):
    reads = [mf.make_read("left", REFERENCE, "del", 1050, 1, start=0, end=1500),
             mf.make_read("right", REFERENCE, "ins", 1300, 1, start=1000, end=4000),
             mf.make_read("inside", REFERENCE, start=1090, end=1110)]
    paths = dataset(tmp_path, reads)
    whole = nonrepetitive_scan(paths)
    parts = [nonrepetitive_scan(paths, start, end) for start, end in [(0, 1100), (1100, 6000)]]

    key = [column for column in rows_of(whole[1])[0] if column != "count"]
    assert [row for calls, _, _ in parts for row in rows_of(calls)] == rows_of(whole[0])
    assert {row["readID"] for row in rows_of(whole[0])} == {"left", "right"}
    assert summed_opportunities([opportunities for _, opportunities, _ in parts], key) \
        == summed_opportunities([whole[1]], key)


# -- tandem repeat tracts -----------------------------------------------------


def tract_opportunities(paths, bed):
    """Opportunity count per period k over the whole contig, for one model."""
    _, opportunities, _ = tandem_repeat_scan({"M1": bed}, paths)
    by_k = {}
    for row in rows_of(opportunities["M1"]):
        by_k[row["k"]] = by_k.get(row["k"], 0) + int(row["count"])
    return {k: count for k, count in by_k.items() if count}


def test_overlapping_tracts_give_calls_the_largest_period_and_opportunities_the_smallest(tmp_path):
    bed = str(tmp_path / "tracts.bed")
    mf.write_tr_bed(bed, [(2100, 2300, "AC", 2), (2190, 2250, "ACAG", 4, 90.0)])
    # A 4 bp deletion inside both tracts fits either period.
    read = mf.make_read("del4", REFERENCE, "del", 2200, 4)
    paths = dataset(tmp_path, [read])
    calls, _, _ = tandem_repeat_scan({"M1": bed}, paths)
    assert [(call["readID"], call["k"]) for call in rows_of(calls["M1"])] == [("del4", "4")]
    assert set(tract_opportunities(paths, bed)) == {"2"}


def test_a_read_ending_inside_a_tract_counts_the_tract_bases_it_covers(tmp_path):
    """
    Today's unit is one opportunity per read base in a tract. Ticket 22
    replaces it with one per tract per read spanning the tract and its flanks.
    """
    bed = str(tmp_path / "tracts.bed")
    mf.write_tr_bed(bed, [(2100, 2300, "AC", 2)])
    without = tract_opportunities(dataset(tmp_path / "without", []), bed)
    partial = mf.make_read("partial", REFERENCE, start=0, end=2150)
    with_partial = tract_opportunities(dataset(tmp_path / "with", [partial]), bed)
    # Tract bases 2100-2148: all it covers except its last base, which has no
    # base after it.
    assert with_partial["2"] - without["2"] == 49
