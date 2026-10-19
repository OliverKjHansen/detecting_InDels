"""
InDels at the edge of a tract (ticket 17: L1 of the workflow 3 review, R1 and
R2 of the workflow 4 review). pbmm2 places an InDel that can slide at its
leftmost position, so slippage in a pure tract is anchored on the base before
the tract, or before a partial unit in front of it.

An event of k bases anchored at base a is workflow 3's if bases a to a+k are
all unmasked, for calls and opportunities alike. Otherwise it is workflow 4's,
in the tract of the first base after the anchor that is in one (looking at
most the longest period on), or else the anchor's. So each event is counted
once, by one workflow.
"""

import pytest

import make_fixtures as mf
from scans import nonrepetitive_scan, rows_of, tandem_repeat_scan

FLANK = ("ACGT" * 250)[:1000]
# (CA)30 at [1001, 1061); the A before it is a partial unit.
REFERENCE = FLANK + "A" + "CA" * 30 + FLANK
TRACT = (1001, 1061, "AC", 2)


def callers(tmp_path, reads, reference=REFERENCE, tracts=(TRACT,)):
    """{read: the workflows that call it} for reads near the tract(s), workflow 3 masking them."""
    paths = mf.write_dataset(str(tmp_path), reference, reads, mask=[tract[:2] for tract in tracts])
    bed = str(tmp_path / "tracts.bed")
    mf.write_tr_bed(bed, list(tracts))
    outside, _, _ = nonrepetitive_scan(paths)
    inside, _, _ = tandem_repeat_scan({"M1": bed}, paths)
    called = {}
    for workflow, table in (("workflow 3", outside), ("workflow 4", inside["M1"])):
        for row in rows_of(table):
            called.setdefault(row["readID"], []).append(workflow)
    return called, rows_of(inside["M1"])


def near_tract(name, kind, position, length=2, reference=REFERENCE):
    return mf.make_read(name, reference, kind, position, length, start=900, end=1160)


def test_slippage_anchored_before_a_tract_is_the_tracts(tmp_path):
    # A CA unit deleted, placed on the tract's first unit (anchor 1000) or on
    # the partial unit before it (anchor 999).
    reads = [near_tract("first_unit", "del", 1001), near_tract("partial_unit", "del", 1000),
             near_tract("insertion", "ins", 1000)]
    called, _ = callers(tmp_path, reads)
    assert called == {"first_unit": ["workflow 4"], "partial_unit": ["workflow 4"], "insertion": ["workflow 4"]}


def test_every_event_near_a_tract_is_counted_by_one_workflow(tmp_path):
    anchors = list(range(990, 1012)) + list(range(1050, 1072))
    reads = ([near_tract(f"del{a}", "del", a + 1) for a in anchors]
             + [near_tract(f"ins{a}", "ins", a) for a in anchors])
    called, _ = callers(tmp_path, reads)
    # Bases a to a+2 touch the tract for anchors 999 to 1060.
    expected = {f"{kind}{a}": ["workflow 4" if 999 <= a <= 1060 else "workflow 3"]
                for kind in ("del", "ins") for a in anchors}
    assert called == expected


def test_slippage_of_the_next_tract_is_that_tracts(tmp_path):
    # (AC)15 then (AT)15: an AT unit deleted, placed leftmost, is anchored on
    # the last base of the AC tract.
    reference = FLANK + "AC" * 15 + "AT" * 15 + FLANK
    tracts = [(1000, 1030, "AC", 2), (1030, 1060, "AT", 2)]
    read = mf.make_read("at", reference, "del", 1030, 2, start=900, end=1160)
    called, calls = callers(tmp_path, [read], reference, tracts)
    assert called == {"at": ["workflow 4"]}
    assert [(call["canonical_motif"], call["indel_string"]) for call in calls] == [("AT", "AT")]


def test_opportunities_test_the_same_bases_against_the_mask(tmp_path):
    # One read over the whole reference. An opportunity of size k for the
    # k-mer or gap at read base p has its anchor at p-1, so it is counted if
    # p-1 to p+k-1 are unmasked: p up to 1001-k before the tract, from 1062
    # after it. With two read bases on each side of the event, p runs from 2,
    # and to 2059-k for a deletion or 2059 for an insertion.
    paths = mf.write_dataset(str(tmp_path), REFERENCE, [mf.make_read("read", REFERENCE)], mask=[TRACT[:2]])
    _, table, _ = nonrepetitive_scan(paths)
    totals = {}
    for row in rows_of(table):
        key = (row["indel_type"], int(row["size"]))
        totals[key] = totals.get(key, 0) + int(row["count"])
    assert totals == {**{("deletion", k): (1000 - k) + (998 - k) for k in (1, 2, 3)},
                      **{("insertion", k): (1000 - k) + 998 for k in (1, 2, 3)}}


@pytest.mark.parametrize("window_end", [1000, 1001, 1002])
def test_windows_still_add_up_at_a_tract_edge(tmp_path, window_end):
    reads = [near_tract("partial_unit", "del", 1000), near_tract("clean", None, None, 0)]
    paths = mf.write_dataset(str(tmp_path), REFERENCE, reads, mask=[TRACT[:2]])
    bed = str(tmp_path / "tracts.bed")
    mf.write_tr_bed(bed, [TRACT])
    whole = tandem_repeat_scan({"M1": bed}, paths)[0]["M1"]
    split = [tandem_repeat_scan({"M1": bed}, paths, start, end)[0]["M1"]
             for start, end in ((0, window_end), (window_end, len(REFERENCE)))]
    assert sorted(rows_of(whole), key=str) == sorted(rows_of(split[0]) + rows_of(split[1]), key=str)
    assert [row["readID"] for row in rows_of(whole)] == ["partial_unit"]
