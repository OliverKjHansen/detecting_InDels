"""
Workflow 4's unit (ticket 22: R6-R8 of the workflow 4 review). An opportunity
is one tract in one read that spans the tract and scan.min_tract_flank bases
on each side; its calls come from the same reads only. A read ending inside
the tract or its flank can't show a change in the tract's length, so it counts
for neither. Each record carries the tract's copies and purity and the read's
flank on each side, binned, so rates by tract length and purity are exact.
Calls and opportunities of a tract are tested on the same bases: the read's
bases from two before the tract to two after it.
"""

import pytest

import make_fixtures as mf
from scans import MIN_TRACT_FLANK, rows_of, tandem_repeat_scan

FLANK = ("ACGT" * 250)[:1000]
REFERENCE = FLANK + "A" + "CA" * 30 + FLANK
TRACT = (1001, 1061, "AC", 2)     # (CA)30, 60 bp


def scan(tmp_path, reads, tracts=(TRACT,), window=(None, None)):
    """(calls, opportunities) of model M1 for these reads, as row dicts."""
    paths = mf.write_dataset(str(tmp_path), REFERENCE, reads)
    bed = str(tmp_path / "tracts.bed")
    mf.write_tr_bed(bed, list(tracts))
    calls, opportunities, _ = tandem_repeat_scan({"M1": bed}, paths, *window)
    return rows_of(calls["M1"]), rows_of(opportunities["M1"])


def count(rows, **match):
    return sum(int(row["count"]) for row in rows if all(row[c] == v for c, v in match.items()))


def read(name, start=900, end=1160, **kwargs):
    """A read over [start, end); by default it spans the tract with about 100 bp on each side."""
    return mf.make_read(name, REFERENCE, start=start, end=end, **kwargs)


def deletion(name, start=900, end=1160, **kwargs):
    """A read with one CA unit deleted inside the tract (anchored on base 1030)."""
    return mf.make_read(name, REFERENCE, "del", 1031, 2, start=start, end=end, **kwargs)


def test_a_read_spanning_the_tract_and_its_flanks_is_one_opportunity(tmp_path):
    calls, opportunities = scan(tmp_path, [read("plain"), deletion("del")])
    assert count(opportunities) == 2
    assert [call["readID"] for call in calls] == ["del"]


@pytest.mark.parametrize("start, end", [(900, 1040), (1020, 1160)])
def test_a_read_ending_inside_the_tract_counts_for_neither(tmp_path, start, end):
    reads = [read("plain", start, end), deletion("del", start, end)]
    calls, opportunities = scan(tmp_path, reads)
    assert (calls, count(opportunities)) == ([], 0)


def test_a_read_ending_inside_a_flank_counts_for_neither(tmp_path):
    short = TRACT[1] + MIN_TRACT_FLANK - 1      # one base short of the flank after the tract
    late = TRACT[0] - MIN_TRACT_FLANK + 1       # one base short before it
    reads = [read("short", end=short), deletion("short_del", end=short),
             read("late", start=late), deletion("late_del", start=late),
             read("enough", TRACT[0] - MIN_TRACT_FLANK, TRACT[1] + MIN_TRACT_FLANK)]
    calls, opportunities = scan(tmp_path, reads)
    assert (calls, count(opportunities)) == ([], 1)


def test_a_tract_is_one_opportunity_per_read_however_long(tmp_path):
    reference_tracts = [(1001, 1061, "AC", 2), (1500, 1520, "ACGT", 4)]
    _, opportunities = scan(tmp_path, [read("whole", 0, 2061)], reference_tracts)
    assert {(row["canonical_motif"], row["count"]) for row in opportunities} == {("AC", "1"), ("ACGT", "1")}


def test_the_tract_and_the_reads_flanks_are_keys_of_both_tables(tmp_path):
    # Flanks of 101 and 99 bp; 30 copies, 100% pure.
    calls, opportunities = scan(tmp_path, [deletion("del")])
    keys = {"tract_copies_bin": "20-50", "tract_purity_bin": ">=100",
            "left_flank_bin": "100-300", "right_flank_bin": "25-100", "read_end_distance_bin": "25-100"}
    assert {column: calls[0][column] for column in keys} == keys
    assert {column: opportunities[0][column] for column in keys} == keys


@pytest.mark.parametrize("low, tested", [(998, False), (999, True), (1030, True), (1062, True), (1063, False)])
def test_a_tracts_calls_and_opportunities_test_the_same_bases(tmp_path, low, tested):
    # The tract's bases and two on each side, 999 to 1062, whichever base the call is at.
    expected = "fail" if tested else "pass"
    calls, opportunities = scan(tmp_path, [read("plain", qualities={low: 30}),
                                           deletion("del", qualities={low: 30})])
    assert [call["kmerQ_pass"] for call in calls] == [expected]
    assert {row["kmerQ_pass"] for row in opportunities} == {expected}


def test_a_tracts_reads_carry_their_nearest_event_outside_the_tract(tmp_path):
    # A mismatch 10 bases after the tract; the call itself isn't another event.
    other = "G" if REFERENCE[1070] != "G" else "T"
    calls, opportunities = scan(tmp_path, [deletion("del", substitutions={1070: other}),
                                           read("plain", substitutions={1070: other}), deletion("alone")])
    assert {call["readID"]: call["nearby_event_bin"] for call in calls} == {"del": "5-25", "alone": "none"}
    assert count(opportunities, nearby_event_bin="5-25") == 2 and count(opportunities, nearby_event_bin="none") == 1


def test_tracts_of_one_model_must_not_overlap(tmp_path):
    with pytest.raises(ValueError, match="overlap"):
        scan(tmp_path, [read("plain")], [(1001, 1061, "AC", 2), (1050, 1080, "ACGT", 4)])


@pytest.mark.parametrize("split", [950, 1001, 1002, 1030, 1061, 1100])
def test_windows_add_up_with_each_tract_in_the_window_of_its_first_base(tmp_path, split):
    reads = [read("plain"), deletion("del"), mf.make_read("before", REFERENCE, "del", 1001, 2, start=900, end=1160)]
    whole = scan(tmp_path / "whole", reads)
    parts = [scan(tmp_path / f"part{i}", reads, window=window)
             for i, window in enumerate([(0, split), (split, len(REFERENCE))])]
    assert sorted(map(str, whole[0])) == sorted(map(str, parts[0][0] + parts[1][0]))
    assert count(whole[1]) == count(parts[0][1]) + count(parts[1][1]) == 3


# -- classes (ticket 23) ------------------------------------------------------


@pytest.mark.parametrize("kind, position, length, size, indel_class", [
    ("del", 1031, 2, "2", "slippage"),       # one CA unit
    ("del", 1032, 4, "4", "slippage"),       # two units, starting on the A: ACAC
    ("ins", 1030, 2, "2", "other"),          # a multiple of 2, but two copies of one base
    ("del", 1031, 1, "1", "other"),          # not a multiple of the period
    ("del", 1031, 3, "3", "other"),
])
def test_every_indel_counted_in_a_tract_is_written_with_its_class_and_size(
        tmp_path, kind, position, length, size, indel_class):
    calls, _ = scan(tmp_path, [mf.make_read("read", REFERENCE, kind, position, length, start=900, end=1160)])
    assert [(call["size"], call["indel_class"]) for call in calls] == [(size, indel_class)]


def test_slippage_is_the_motif_on_either_strand_whatever_the_case(tmp_path):
    # A lower-case (soft-masked) tract of (TG)n is the same tract as (CA)n.
    reference = FLANK + "A" + "tg" * 30 + FLANK
    paths = mf.write_dataset(str(tmp_path), reference,
                             [mf.make_read("read", reference.upper(), "del", 1031, 2, start=900, end=1160)])
    bed = str(tmp_path / "tracts.bed")
    mf.write_tr_bed(bed, [TRACT])
    calls, _, _ = tandem_repeat_scan({"M1": bed}, paths)
    assert [call["indel_class"] for call in rows_of(calls["M1"])] == ["slippage"]
