"""
End-to-end tests of the pileup-scanning callers against tests/fixtures/.

The fixture reads carry InDels planted at known positions (see
tests/make_fixtures.py), so these assert an exact call set rather than merely
that output was produced. None of this needs a cluster, and none of it needs the
external binaries - it is the part of the pipeline that can genuinely be
developed on a laptop.
"""

import os
import subprocess
import sys

import pysam
import pytest

import utilities
from scans import nonrepetitive_scan, rows_of
from scans import tandem_repeat_scan as scan

CONTIG = "testcontig1"


@pytest.fixture(scope="module")
def small_indel_calls(fixture_paths):
    calls, _, _ = nonrepetitive_scan(fixture_paths)
    return rows_of(calls)


@pytest.fixture(scope="module")
def tr_indel_calls(fixture_paths):
    calls, _, _ = scan({"M1": fixture_paths["tr_bed"]}, fixture_paths)
    return rows_of(calls["M1"])


# -- workflow B: InDels outside repeats -------------------------------------


def test_small_indels_recovers_exactly_the_expected_reads(small_indel_calls):
    assert {call["readID"] for call in small_indel_calls} == {
        "read_del1",
        "read_ins1",
        "read_del2",
    }


def test_small_indels_reports_the_planted_positions(small_indel_calls):
    by_read = {call["readID"]: call for call in small_indel_calls}
    # `position` is 1-based, so it is the planted 0-based index + 1 for an
    # insertion, and the index of the base before the gap for a deletion.
    assert by_read["read_del1"]["position"] == "500"
    assert by_read["read_ins1"]["position"] == "1004"
    assert by_read["read_del2"]["position"] == "1500"


def test_small_indels_reports_the_right_types_and_sizes(small_indel_calls):
    by_read = {call["readID"]: call for call in small_indel_calls}
    assert by_read["read_del1"]["indel_type"] == "deletion"
    assert by_read["read_ins1"]["indel_type"] == "insertion"
    assert by_read["read_del2"]["indel_type"] == "deletion"
    assert len(by_read["read_del2"]["indel_string"]) == 2


def test_homopolymer_indels_are_dropped(small_indel_calls):
    """read_del_homo sits in the poly-A tract and must not be called."""
    assert "read_del_homo" not in {call["readID"] for call in small_indel_calls}


def test_masked_regions_are_skipped(small_indel_calls):
    """
    read_tr_del1 and read_tr_del2 sit inside the microsatellite interval listed
    in the filter BED, so workflow B must not see them at all.
    """
    called = {call["readID"] for call in small_indel_calls}
    assert "read_tr_del1" not in called
    assert "read_tr_del2" not in called


def test_single_read_indels_are_marked_unique(small_indel_calls):
    assert all(call["unique"] == "unique" for call in small_indel_calls)


def test_tissue_is_joined_onto_every_call(small_indel_calls):
    assert all(call["tissue"] in {"liver", "brain"} for call in small_indel_calls)


def test_pacbio_tags_are_carried_through(small_indel_calls):
    """np/ec/sa come from the BAM tags and feed the later quality modelling."""
    for call in small_indel_calls:
        assert call["np"] == "20"
        assert call["sa"] == "10"
        assert call["kmerQ"] == "93"


# -- workflow C: InDels inside repeats --------------------------------------


def test_tandem_repeat_caller_keeps_only_period_multiples(tr_indel_calls):
    """
    Inside a k=2 repeat, the 2 bp deletion is a whole-motif contraction and is
    kept; the 1 bp deletion is not a multiple of the period and is dropped.
    """
    assert {call["readID"] for call in tr_indel_calls} == {"read_tr_del2"}


def test_tandem_repeat_calls_carry_the_repeat_annotation(tr_indel_calls):
    call = tr_indel_calls[0]
    assert call["canonical_motif"] == "AC"
    assert call["k"] == "2"
    assert call["model"] == "M1"
    assert float(call["model_purity"]) == 100.0


def test_tandem_repeat_caller_ignores_positions_outside_the_bed(tr_indel_calls):
    """The InDels in the complex region are outside every TR interval."""
    called = {call["readID"] for call in tr_indel_calls}
    assert called.isdisjoint({"read_del1", "read_ins1", "read_del2"})


def test_homopolymer_filter_is_not_applied_inside_repeats(tr_indel_calls):
    """
    k=1 slippage is the signal workflow C measures, so it must not reuse
    workflow B's homopolymer filter. The (CA)n call surviving with a fully
    repetitive context demonstrates the filter is off.
    """
    assert tr_indel_calls
    assert tr_indel_calls[0]["indel_string"] in {"CA", "AC"}


# -- opportunities ----------------------------------------------------------


def test_opportunities_are_counted_outside_the_mask(fixture_paths):
    _, opportunities, _ = nonrepetitive_scan(fixture_paths)
    rows = rows_of(opportunities)
    assert rows
    assert {"context", "tissue", "mapQ", "kmerQ_pass", "count"} <= set(rows[0])
    assert sum(int(row["count"]) for row in rows) > 0


def test_tandem_repeat_opportunities_are_binned_by_motif(fixture_paths):
    _, opportunities, _ = scan({"M1": fixture_paths["tr_bed"]}, fixture_paths)
    rows = rows_of(opportunities["M1"])
    counted = [row for row in rows if int(row["count"]) > 0]
    assert counted
    assert all(row["canonical_motif"] == "AC" and row["k"] == "2" for row in counted)


# -- one pass over the BAM for several repeat models ------------------------


@pytest.fixture(scope="module")
def three_model_beds(tmp_path_factory):
    """
    Three repeat models over the fixture that disagree: M0 covers only part of
    the (CA)n tract, M1 all of it, and M2 all of it plus an overlapping k=4
    interval, so the largest-k (calls) and smallest-k (opportunities) choices
    come into play.
    """
    directory = tmp_path_factory.mktemp("tr_beds")
    header = "contig\tstart\tend\tcanonical_motif\tk\tcopies\ttrf_purity\tmodel_purity\tmodel\n"
    rows = {
        "M0": [f"{CONTIG}\t2150\t2300\tAC\t2\t75\t100.00\t100.00\tM0\n"],
        "M1": [f"{CONTIG}\t2100\t2300\tAC\t2\t100\t100.00\t100.00\tM1\n"],
        "M2": [f"{CONTIG}\t2100\t2300\tAC\t2\t100\t100.00\t100.00\tM2\n",
               f"{CONTIG}\t2190\t2250\tACAG\t4\t15\t90.00\t90.00\tM2\n"],
    }
    beds = {}
    for model, lines in rows.items():
        beds[model] = str(directory / f"{CONTIG}_{model}.bed")
        with open(beds[model], "w") as bed:
            bed.write(header + "".join(lines))
    return beds


def test_scanning_models_together_matches_scanning_each_alone(three_model_beds, fixture_paths):
    calls, opportunities, _ = scan(three_model_beds, fixture_paths)
    for model, bed in three_model_beds.items():
        alone_calls, alone_opportunities, _ = scan({model: bed}, fixture_paths)
        assert calls[model] == alone_calls[model]
        assert opportunities[model] == alone_opportunities[model]
    assert all(rows_of(calls[model]) for model in three_model_beds)


def test_scan_writes_each_read_carrying_a_call_once(three_model_beds, fixture_paths):
    calls, _, fastq = scan(three_model_beds, fixture_paths)
    called = {row["readID"] for model in calls for row in rows_of(calls[model])}
    lines = fastq.strip().split("\n")
    names = [line[1:] for line in lines[0::4]]
    assert sorted(names) == sorted(called)

    with pysam.AlignmentFile(fixture_paths["bam"]) as bam:
        reads = {read.query_name: read for read in bam.fetch(CONTIG) if read.query_name in called}
    for i, name in enumerate(names):
        read = reads[name]
        assert lines[4 * i + 1] == read.get_forward_sequence()
        assert lines[4 * i + 3] == "".join(chr(q + 33) for q in read.get_forward_qualities())


def test_scan_windows_add_up_to_the_whole_contig(three_model_beds, fixture_paths):
    whole_calls, whole_opportunities, whole_fastq = scan(three_model_beds, fixture_paths)
    windows = [(0, 2150), (2150, 2204), (2204, 2205), (2205, 6000)]
    parts = [scan(three_model_beds, fixture_paths, start, end) for start, end in windows]

    def summed(texts):
        totals = {}
        for text in texts:
            for row in rows_of(text):
                key = tuple(value for column, value in row.items() if column != "count")
                totals[key] = totals.get(key, 0) + int(row["count"])
        return {key: count for key, count in totals.items() if count}

    for model in three_model_beds:
        assert [row for calls, _, _ in parts for row in rows_of(calls[model])] == rows_of(whole_calls[model])
        assert summed(opportunities[model] for _, opportunities, _ in parts) == summed([whole_opportunities[model]])
    assert "".join(fastq for _, _, fastq in parts) == whole_fastq


def test_haploid_filter_checks_several_call_tables_against_one_remapping(tmp_path):
    """Several tables in one run give the same rows as one run per table, with one header."""
    header = "#contig\tstart\tend\tindel_type\tindel_string\treadID\tdiploid_direction\tread_length\n"
    tables = {
        "a.tsv": header + "c1\t11\t11\tdeletion\tAC\tr1\t+\t100\nc1\t41\t41\tdeletion\tAC\tr3\t+\t100\n",
        "b.tsv": header + "c2\t21\t23\tinsertion\tTT\tr2\t-\t100\n",
    }
    for name, text in tables.items():
        (tmp_path / name).write_text(text)
    (tmp_path / "locations.tsv").write_text("r1\t0\th1\t5\t60\t10M2D90M\nr2\t16\th1\t5\t60\t78M2I20M\n")

    script = os.path.join(os.path.dirname(utilities.__file__), "filtering_haploid_mappping.py")

    def run(*tables):
        return subprocess.run([sys.executable, script, *[str(tmp_path / t) for t in tables],
                               str(tmp_path / "locations.tsv")],
                              capture_output=True, text=True, check=True).stdout

    separate = [run(name).split("\n", 1) for name in tables]
    assert run(*tables) == separate[0][0] + "\n" + "".join(body for _, body in separate)


# -- one pass over the BAM outside repeats ----------------------------------


def test_nonrepetitive_scan_writes_each_read_carrying_a_call_once(fixture_paths):
    calls, _, fastq = nonrepetitive_scan(fixture_paths)
    called = {row["readID"] for row in rows_of(calls)}
    assert called == {"read_del1", "read_ins1", "read_del2"}
    assert sorted(line[1:] for line in fastq.strip().split("\n")[0::4]) == sorted(called)


def test_nonrepetitive_scan_windows_add_up_to_the_whole_contig(fixture_paths):
    whole_calls, whole_opportunities, whole_fastq = nonrepetitive_scan(fixture_paths)
    parts = [nonrepetitive_scan(fixture_paths, start, end)
             for start, end in [(0, 499), (499, 1003), (1003, 1004), (1004, 6000)]]

    def summed(texts):
        totals = {}
        for text in texts:
            for row in rows_of(text):
                key = tuple(value for column, value in row.items() if column != "count")
                totals[key] = totals.get(key, 0) + int(row["count"])
        return {key: count for key, count in totals.items() if count}

    assert [row for calls, _, _ in parts for row in rows_of(calls)] == rows_of(whole_calls)
    assert summed(opportunities for _, opportunities, _ in parts) == summed([whole_opportunities])
    assert "".join(fastq for _, _, fastq in parts) == whole_fastq


# -- parameters from config.yaml ---------------------------------------------


def test_positions_outside_the_coverage_band_are_skipped(fixture_paths):
    # Every fixture position has exactly the mean depth, so a band above the
    # mean leaves nothing to call or count.
    calls, opportunities, _ = nonrepetitive_scan(fixture_paths, coverage_band=(1.5, 3.0))
    assert rows_of(calls) == [] and rows_of(opportunities) == []
    calls, opportunities, _ = scan({"M1": fixture_paths["tr_bed"]}, fixture_paths, coverage_band=(1.5, 3.0))
    assert rows_of(calls["M1"]) == [] and rows_of(opportunities["M1"]) == []


def test_only_indels_of_the_configured_sizes_are_called_and_counted(fixture_paths):
    calls, opportunities, _ = nonrepetitive_scan(fixture_paths, indel_sizes=(2, 3))
    assert {call["readID"] for call in rows_of(calls)} == {"read_del2"}
    assert {len(row["context"]) for row in rows_of(opportunities)} == {2, 3}
