"""
Denominators that match the calls (ticket 16: D1, D3 and D4 of
documents/nonrepetitive_indel_review.md), outside repeats.

Calls and opportunities alike are tested on two read bases on each side of
their event, never on the event's own bases: a deletion can't test the bases
it removes, so neither does its opportunity, and an insertion's opportunity
has no inserted bases to test. An event without two read bases on a side is
neither called nor counted. Insertions have their own opportunities, at every
gap between two read bases for each size, keyed by the two bases around the
gap.
"""

import pytest

import make_fixtures as mf
from scans import nonrepetitive_scan, rows_of, summed

REFERENCE = mf.build_reference()


def dataset(tmp_path, reads, clean=0):
    reads = reads + [mf.make_read(f"clean{i:02d}", REFERENCE) for i in range(clean)]
    return mf.write_dataset(str(tmp_path), REFERENCE, reads)


def opportunities(paths, **match):
    """{(indel_type, size): count} over the rows with these values."""
    _, table, _ = nonrepetitive_scan(paths)
    totals = {}
    for row in rows_of(table):
        if all(row[column] == value for column, value in match.items()):
            key = (row["indel_type"], int(row["size"]))
            totals[key] = totals.get(key, 0) + int(row["count"])
    return totals


def test_each_size_is_counted_where_its_event_has_two_read_bases_on_each_side(tmp_path):
    # One read over [0, 1000), the ACGT region, where no k-mer is a homopolymer.
    # Deleting the k-mer at read base p keeps p-2, p-1, p+k and p+k+1, so p runs
    # from 2 to 997-k (D3: none is cut short at the read's end). Inserting
    # between bases g-1 and g keeps g-2, g-1, g and g+1, so g runs from 2 to
    # 998, whatever the size.
    paths = dataset(tmp_path, [mf.make_read("read", REFERENCE, end=1000)])
    assert opportunities(paths) == {("deletion", 1): 996, ("deletion", 2): 995, ("deletion", 3): 994,
                                    ("insertion", 1): 997, ("insertion", 2): 997, ("insertion", 3): 997}


def test_a_deletions_context_is_its_k_mer(tmp_path):
    _, table, _ = nonrepetitive_scan(dataset(tmp_path, [], clean=1))
    assert all(len(row["context"]) == int(row["size"]) for row in rows_of(table) if row["indel_type"] == "deletion")


def test_an_insertions_context_is_the_pair_of_bases_around_its_gap(tmp_path):
    _, table, _ = nonrepetitive_scan(dataset(tmp_path, [], clean=1))
    contexts = {row["context"] for row in rows_of(table) if row["indel_type"] == "insertion"}
    # All 16 pairs occur in the contig; each is grouped with the same gap read
    # from the other strand.
    assert contexts == {"AA", "AC", "AG", "AT", "CA", "CC", "CG", "GA", "GC", "TA"}


def test_a_low_quality_base_fails_the_four_opportunities_that_test_it(tmp_path):
    # Base 500 at Q30: four deletions and four gaps of each size test it.
    paths = dataset(tmp_path, [mf.make_read("low", REFERENCE, qualities={500: 30})])
    assert opportunities(paths, kmerQ_pass="fail") == {("deletion", 1): 4, ("deletion", 2): 4, ("deletion", 3): 4,
                                                       ("insertion", 1): 4, ("insertion", 2): 4, ("insertion", 3): 4}


def kmerq_pass(tmp_path, read):
    calls, _, _ = nonrepetitive_scan(dataset(tmp_path, [read], clean=4))
    return [row["kmerQ_pass"] for row in rows_of(calls) if row["readID"] == read.query_name]


@pytest.mark.parametrize("low, tested", [(497, False), (498, True), (499, True), (501, True), (502, True),
                                         (503, False)])
def test_a_deletion_call_tests_two_read_bases_on_each_side_of_its_gap(tmp_path, low, tested):
    # Base 500 deleted: the read runs ... 498 499 | 501 502 ...
    read = mf.make_read("del", REFERENCE, "del", 500, 1, qualities={low: 30})
    assert kmerq_pass(tmp_path, read) == ["fail" if tested else "pass"]


@pytest.mark.parametrize("low, tested", [(1001, False), (1002, True), (1003, True), ("inserted", False),
                                         (1004, True), (1005, True), (1006, False)])
def test_an_insertion_call_tests_two_read_bases_on_each_side_of_its_gap(tmp_path, low, tested):
    # A base inserted after 1003: the read runs ... 1002 1003 | 1004 1005 ...
    read = mf.make_read("ins", REFERENCE, "ins", 1003, 1, qualities={low: 30})
    assert kmerq_pass(tmp_path, read) == ["fail" if tested else "pass"]


def test_a_call_without_two_read_bases_on_a_side_is_not_called(tmp_path):
    # The read starts at the deletion's anchor: one read base before the gap.
    assert kmerq_pass(tmp_path, mf.make_read("edge", REFERENCE, "del", 500, 1, start=499)) == []
    assert kmerq_pass(tmp_path / "inside", mf.make_read("inside", REFERENCE, "del", 500, 1, start=498)) == ["pass"]


# -- coverage at the anchor (ticket 31) --------------------------------------
#
# Ten reads up to base 999, four of them on to the end of the contig: depth 10
# up to base 999 and 4 after it, against a mean of 10, so the coverage band is
# 5-20. Only bases 996-1003 are unmasked.

STEP_MASK = [(0, 996), (1004, mf.CONTIG_LENGTH)]


def step_reads(*extra, full=4):
    """Ten reads: `full` over the whole contig, the extra ones, and the rest up to base 999."""
    return ([mf.make_read(f"full/{i}/ccs", REFERENCE) for i in range(full)]
            + [mf.make_read(f"short/{i}/ccs", REFERENCE, end=1000) for i in range(10 - full - len(extra))]
            + list(extra))


def test_an_opportunity_is_counted_where_its_anchors_depth_is_in_the_band(tmp_path):
    # An event of k bases anchored on 996 to 1003-k touches only unmasked
    # bases; of those, the ones anchored up to 999 have their anchor in the
    # band: four per size, on each of the four full reads.
    paths = mf.write_dataset(str(tmp_path), REFERENCE, step_reads(), mask=STEP_MASK)
    assert opportunities(paths, smrt_cell="full") == {(kind, k): 16 for kind in ("deletion", "insertion")
                                                      for k in (1, 2, 3)}


def test_an_opportunity_is_not_counted_where_only_its_k_mers_depth_is_in_the_band(tmp_path):
    # The other way round: four full reads, and six from base 1000 on, so
    # depth 4 up to base 999 and 10 after it. Anchors from 1000 to 1003-k are
    # in the band: three, two and one per full read.
    reads = ([mf.make_read(f"full/{i}/ccs", REFERENCE) for i in range(4)]
             + [mf.make_read(f"late/{i}/ccs", REFERENCE, start=1000) for i in range(6)])
    paths = mf.write_dataset(str(tmp_path), REFERENCE, reads, mask=STEP_MASK)
    assert opportunities(paths, smrt_cell="full") == {(kind, k): 4 * (4 - k) for kind in ("deletion", "insertion")
                                                      for k in (1, 2, 3)}


def test_a_call_is_counted_where_its_anchors_depth_is_in_the_band(tmp_path):
    # Two full reads and the two with calls go past base 999: depth 4 there.
    reads = step_reads(mf.make_read("in_band/1/ccs", REFERENCE, "del", 1000, 1, end=1100),
                       mf.make_read("out_of_band/1/ccs", REFERENCE, "del", 1001, 1, start=990), full=2)
    calls, _, _ = nonrepetitive_scan(mf.write_dataset(str(tmp_path), REFERENCE, reads, mask=STEP_MASK))
    assert [row["readID"] for row in rows_of(calls)] == ["in_band/1/ccs"]


@pytest.mark.parametrize("split", [999, 1000, 1001])
def test_windows_still_add_up_where_the_depth_steps(tmp_path, split):
    paths = mf.write_dataset(str(tmp_path), REFERENCE, step_reads(), mask=STEP_MASK)
    tables = [nonrepetitive_scan(paths, *window)[1] for window in [(None, None), (0, split), (split, mf.CONTIG_LENGTH)]]
    assert summed(tables[1:]) == summed(tables[:1])
