"""
Subread support (S1): the per-base tags sm, sx and sa are stored in the order
the read was sequenced, so on a reverse-strand read the scans must read them
from the other end. Calls carry the corrected tags; calls and opportunities
both carry the binned support of the bases they test.
"""

import make_fixtures as mf
from scans import nonrepetitive_scan, rows_of, tandem_repeat_scan

REFERENCE = mf.build_reference()


def dataset(tmp_path, reads):
    clean = [mf.make_read(f"clean{i:02d}", REFERENCE) for i in range(10)]
    return mf.write_dataset(str(tmp_path), REFERENCE, reads + clean)


def call_of(tmp_path, read):
    calls, _, _ = nonrepetitive_scan(dataset(tmp_path, [read]))
    return next(row for row in rows_of(calls) if row["readID"] == read.query_name)


# A 1 bp insertion after base 1003: the anchor is base 1003, the base after it
# the inserted one.
PLANTED = {1003: (9, 1, 11), "inserted": (4, 6, 12)}


def test_a_reverse_strand_call_reports_the_subread_tags_of_its_own_bases(tmp_path):
    read = mf.make_read("rev", REFERENCE, "ins", 1003, 1, reverse=True, subreads=PLANTED)
    call = call_of(tmp_path, read)
    # sm and sx at the base after the anchor; sa, sm_pre and sx_pre at the anchor.
    assert (call["sm"], call["sx"], call["sa"], call["sm_pre"], call["sx_pre"]) == ("4", "6", "11", "9", "1")


def test_a_forward_strand_call_reports_the_same(tmp_path):
    read = mf.make_read("fwd", REFERENCE, "ins", 1003, 1, subreads=PLANTED)
    call = call_of(tmp_path, read)
    assert (call["sm"], call["sx"], call["sa"], call["sm_pre"], call["sx_pre"]) == ("4", "6", "11", "9", "1")


def test_a_call_carries_the_subread_support_of_the_bases_it_tests(tmp_path):
    # The anchor's 9 of 11 is the lowest of the bases tested; the inserted
    # base's 4 of 12 isn't tested (D4).
    read = mf.make_read("rev", REFERENCE, "ins", 1003, 1, reverse=True, subreads=PLANTED)
    assert call_of(tmp_path, read)["subread_support"] == "0.8-0.9"
    clean_call = mf.make_read("plain", REFERENCE, "ins", 1003, 1, reverse=True)
    assert call_of(tmp_path / "plain", clean_call)["subread_support"] == ">=1"


def opportunities_by_support(tmp_path, read):
    _, opportunities, _ = nonrepetitive_scan(dataset(tmp_path, [read]))
    totals = {}
    for row in rows_of(opportunities):
        totals[row["subread_support"]] = totals.get(row["subread_support"], 0) + int(row["count"])
    return totals


def test_opportunities_testing_a_poorly_supported_base_fall_in_its_bin_on_either_strand(tmp_path):
    # Half the subreads match base 1500 of one read. Each opportunity whose
    # tested bases include it, and only those, falls below 0.8.
    for reverse in (False, True):
        read = mf.make_read("low", REFERENCE, reverse=reverse, subreads={1500: (5, 3, 10)})
        totals = opportunities_by_support(tmp_path / str(reverse), read)
        assert totals["<0.8"] == 24, (reverse, totals)   # four deletions and four gaps per size
        assert set(totals) == {"<0.8", ">=1"}, (reverse, totals)


def test_reads_without_subread_tags_have_no_support_rather_than_none(tmp_path):
    read = mf.make_read("sequel", REFERENCE, "ins", 1003, 1, subread_tags=False)
    call = call_of(tmp_path, read)
    assert (call["sm"], call["sa"], call["subread_support"]) == ("NA", "NA", "NA")
    assert opportunities_by_support(tmp_path / "opportunities", read)["NA"] > 0


def test_the_repeat_scan_reads_the_tags_the_same_way(tmp_path):
    bed = str(tmp_path / "tracts.bed")
    mf.write_tr_bed(bed, [(2100, 2300, "AC", 2)])
    # A 2 bp deletion inside the (CA)n tract: the anchor is base 2203, the
    # base after the gap 2206.
    read = mf.make_read("rev", REFERENCE, "del", 2204, 2, reverse=True,
                        subreads={2203: (9, 1, 11), 2206: (4, 6, 12)})
    calls, _, _ = tandem_repeat_scan({"M1": bed}, dataset(tmp_path, [read]))
    call = next(row for row in rows_of(calls["M1"]) if row["readID"] == "rev")
    assert (call["sm"], call["sx"], call["sa"], call["sm_pre"], call["sx_pre"]) == ("4", "6", "11", "9", "1")
