"""
Read quality and context (ticket 14): the read's mapping quality to the diploid
assembly (N1), its passes, its predicted accuracy, its full passes on the
weaker strand of the molecule, and its nearest other InDel or mismatch are keys
on calls and opportunities alike. Calls also carry the assembly sequence
around their anchor.
"""

import make_fixtures as mf
from scans import FLANK_WIDTH, KEY_PARAMETERS, nonrepetitive_scan, rows_of, tandem_repeat_scan

REFERENCE = mf.build_reference()


def dataset(tmp_path, reads):
    clean = [mf.make_read(f"clean{i:02d}", REFERENCE) for i in range(10)]
    return mf.write_dataset(str(tmp_path), REFERENCE, reads + clean)


def call_of(tmp_path, read):
    calls, _, _ = nonrepetitive_scan(dataset(tmp_path, [read]))
    return next(row for row in rows_of(calls) if row["readID"] == read.query_name)


def opportunity_totals(paths, column):
    _, opportunities, _ = nonrepetitive_scan(paths)
    totals = {}
    for row in rows_of(opportunities):
        totals[row[column]] = totals.get(row[column], 0) + int(row["count"])
    return totals


def test_the_bins_come_from_the_key_parameters():
    assert KEY_PARAMETERS["mapq_bins"] == [1, 20, 60]
    assert KEY_PARAMETERS["passes_bins"] == [5, 10, 15, 20, 30]
    assert KEY_PARAMETERS["accuracy_bins"] == [0.99, 0.999, 0.9999]
    assert KEY_PARAMETERS["strand_passes_bins"] == [2, 5, 10]
    assert KEY_PARAMETERS["nearby_event_bins"] == [5, 25, 100]


def test_a_call_carries_its_reads_mapq_passes_and_accuracy(tmp_path):
    read = mf.make_read("low", REFERENCE, "del", 500, 1, mapq=5, tags={"np": 12, "rq": 0.995, "fn": 7, "rn": 1})
    call = call_of(tmp_path, read)
    assert (call["mapQ"], call["mapQ_bin"]) == ("5", "1-20")
    assert (call["np"], call["np_bin"], call["rq_bin"]) == ("12", "10-15", "0.99-0.999")
    assert (call["fn"], call["rn"], call["strand_passes_bin"]) == ("7", "1", "<2")


def test_opportunities_split_by_mapq_without_changing_their_total(tmp_path):
    paths = dataset(tmp_path, [mf.make_read("mapq5", REFERENCE, mapq=5)])
    by_mapq = opportunity_totals(paths, "mapQ_bin")
    # The mapQ 5 read counts the same positions as each of the ten at mapQ 60.
    assert set(by_mapq) == {"1-20", ">=60"}
    assert by_mapq["1-20"] * 10 == by_mapq[">=60"]
    assert sum(by_mapq.values()) == sum(opportunity_totals(paths, "np_bin").values())


def test_reads_without_a_tag_have_no_value_rather_than_none(tmp_path):
    read = mf.make_read("bare", REFERENCE, "del", 500, 1, tags={"np": None, "rq": None, "fn": None, "rn": None},
                        subread_tags=False)
    call = call_of(tmp_path, read)
    assert (call["np"], call["np_bin"], call["rq_bin"], call["strand_passes_bin"]) == ("NA", "NA", "NA", "NA")
    assert opportunity_totals(dataset(tmp_path / "opportunities", [read]), "np_bin")["NA"] > 0


def test_a_call_carries_the_distance_to_its_reads_nearest_other_event(tmp_path):
    # A mismatch 3 bases after the deletion's anchor (base 499).
    near = mf.make_read("near", REFERENCE, "del", 500, 1, substitutions={502: "T" if REFERENCE[502] != "T" else "G"})
    assert call_of(tmp_path, near)["nearby_event_bin"] == "<5"
    # Its own InDel is not another event.
    alone = mf.make_read("alone", REFERENCE, "del", 500, 1)
    assert call_of(tmp_path / "alone", alone)["nearby_event_bin"] == "none"


def test_opportunities_near_another_event_in_the_read_fall_in_its_bin(tmp_path):
    read = mf.make_read("mismatch", REFERENCE, substitutions={1500: "T" if REFERENCE[1500] != "T" else "G"})
    nearby = opportunity_totals(dataset(tmp_path, [read]), "nearby_event_bin")
    # Clean reads have no events; the one read's anchors 1-4 bases from the
    # mismatch, on either side, fall below 5, about 8 anchors for each size
    # of deletion and of insertion.
    assert set(nearby) == {"none", "<5", "5-25", "25-100", ">=100"}
    assert 0 < nearby["<5"] <= 6 * 8


def test_a_call_carries_the_assembly_around_its_anchor(tmp_path):
    call = call_of(tmp_path, mf.make_read("del", REFERENCE, "del", 500, 1))
    assert call["flank_sequence"] == REFERENCE[499 - FLANK_WIDTH:499 + FLANK_WIDTH + 1]


def test_the_repeat_scan_carries_the_same_keys(tmp_path):
    bed = str(tmp_path / "tracts.bed")
    mf.write_tr_bed(bed, [(2100, 2300, "AC", 2)])
    calls, opportunities, _ = tandem_repeat_scan({"M1": bed}, dataset(tmp_path, [
        mf.make_read("del", REFERENCE, "del", 2204, 2, mapq=5)]))
    for row in rows_of(calls["M1"]) + rows_of(opportunities["M1"]):
        assert {"mapQ_bin", "np_bin", "rq_bin", "strand_passes_bin", "nearby_event_bin"} <= set(row)
    assert rows_of(calls["M1"])[0]["mapQ"] == "5"
