"""
Read geometry (ticket 13): read length, distance to the read's end, the strand
the read is aligned to, and its SMRT Cell are keys on calls and opportunities
alike, and base quality is pass/fail at the configured threshold.
"""

import make_fixtures as mf
from scans import KEY_PARAMETERS, nonrepetitive_scan, rows_of, tandem_repeat_scan

REFERENCE = mf.build_reference()
MOVIE = "m84108_250930_153107_s2"


def dataset(tmp_path, reads):
    clean = [mf.make_read(f"{MOVIE}/{i}/ccs", REFERENCE) for i in range(10)]
    return mf.write_dataset(str(tmp_path), REFERENCE, reads + clean)


def calls_by_read(paths):
    calls, _, _ = nonrepetitive_scan(paths)
    return {row["readID"]: row for row in rows_of(calls)}


def opportunity_totals(paths, column):
    _, opportunities, _ = nonrepetitive_scan(paths)
    totals = {}
    for row in rows_of(opportunities):
        totals[row[column]] = totals.get(row[column], 0) + int(row["count"])
    return totals


def test_the_bins_come_from_the_key_parameters():
    assert KEY_PARAMETERS["read_length_bins"] == [4000, 8000, 12000, 16000, 20000]
    assert KEY_PARAMETERS["read_end_distance_bins"] == [25, 100, 300, 1000]
    assert KEY_PARAMETERS["min_kmerq"] == 93


def test_a_call_carries_its_reads_length_strand_and_smrt_cell(tmp_path):
    short = mf.make_read("m64101_220701_134454/7/ccs", REFERENCE, "del", 500, 1, start=0, end=3000, reverse=True)
    call = calls_by_read(dataset(tmp_path, [short]))["m64101_220701_134454/7/ccs"]
    assert (call["read_length_bin"], call["aligned_strand"], call["smrt_cell"]) == ("<4000", "-", "m64101_220701_134454")


def test_a_call_carries_its_distance_to_the_read_end_and_its_raw_and_pass_fail_quality(tmp_path):
    # The anchor, base 499, is 99 bases from the read's first aligned base.
    near_end = mf.make_read("near", REFERENCE, "del", 500, 1, start=400, end=3000, qualities={501: 40})
    call = calls_by_read(dataset(tmp_path, [near_end]))["near"]
    assert (call["read_end_distance_bin"], call["kmerQ"], call["kmerQ_pass"]) == ("25-100", "40", "fail")


def test_opportunities_split_by_read_length_strand_and_smrt_cell_without_changing_their_total(tmp_path):
    reads = [mf.make_read("m64101_220701_134454/7/ccs", REFERENCE, start=0, end=3000, reverse=True)]
    paths = dataset(tmp_path, reads)
    lengths = opportunity_totals(paths, "read_length_bin")
    strands = opportunity_totals(paths, "aligned_strand")
    cells = opportunity_totals(paths, "smrt_cell")
    assert set(lengths) == {"<4000", "4000-8000"}
    assert sum(lengths.values()) == sum(strands.values()) == sum(cells.values())
    assert lengths["<4000"] == strands["-"] == cells["m64101_220701_134454"]


def test_opportunities_near_a_read_end_fall_in_the_shortest_distance_bins(tmp_path):
    distances = opportunity_totals(dataset(tmp_path, []), "read_end_distance_bin")
    # Ten reads spanning the contig: on each, the anchors 0-24 bases from
    # either end are in the first bin, and there are about 2 x 25 of them per
    # size, so a few hundred opportunities.
    assert 0 < distances["<25"] < distances["25-100"] < distances[">=1000"]


def test_opportunities_are_pass_or_fail_at_the_threshold(tmp_path):
    low = mf.make_read("low", REFERENCE, qualities={1500: 40})
    passes = opportunity_totals(dataset(tmp_path, [low]), "kmerQ_pass")
    assert set(passes) == {"pass", "fail"}
    assert passes["fail"] == 24   # the windows that include base 1500: four deletions and four gaps per size


def test_the_repeat_scan_carries_the_same_keys(tmp_path):
    bed = str(tmp_path / "tracts.bed")
    mf.write_tr_bed(bed, [(2100, 2300, "AC", 2)])
    calls, opportunities, _ = tandem_repeat_scan({"M1": bed}, dataset(tmp_path, [
        mf.make_read("del", REFERENCE, "del", 2204, 2)]))
    for row in rows_of(calls["M1"]) + rows_of(opportunities["M1"]):
        assert {"read_length_bin", "aligned_strand", "smrt_cell", "read_end_distance_bin", "kmerQ_pass"} <= set(row)
