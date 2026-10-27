"""
Why workflow 4's step 0 drops tracts, and what each rule of ticket 34 would
measure instead: tract_drop_causes.py on a sample with one contig, its repeat
BEDs, excluded group and unresolved regions written where workflows 1 and 2
leave them.
"""

import contextlib
import io
import os
import runpy
import sys

import targets
from targets import CONTIG, SAMPLE, SPECIES
from workflow_config import REPO_ROOT

MODELS = ["M0", "M1", "M2"]   # the models config.yaml lists
# Six 40 bp M1 tracts, each with one way of meeting the regions below.
TRACTS = [(100, 140, "AC", 2),   # meets nothing
          (200, 240, "AC", 2),   # an unresolved region from base 230
          (300, 340, "AC", 2),   # a homopolymer run over its last 2 bases and beyond
          (400, 440, "AC", 2),   # inside a long-period tract
          (700, 740, "AC", 2),   # a homopolymer run inside it, away from its ends
          (800, 840, "AC", 2)]   # a long-period tract over its last 5 bases and beyond
HOMOPOLYMERS = [(338, 345), (715, 722)]
LONG_PERIODS = [(380, 600), (835, 900)]
UNRESOLVED = [(230, 300)]


def bed(regions, *columns):
    """BED lines of the contig's (start, end) regions, with the extra columns on each."""
    return ["\t".join([CONTIG, str(start), str(end), *columns]) for start, end in regions]


def sample(tmp_path, monkeypatch, tracts=TRACTS, homopolymers=HOMOPOLYMERS, long_periods=LONG_PERIODS,
           unresolved=UNRESOLVED, timed_out=()):
    """
    The sample's M1 tracts and regions, as workflows 1 and 2 leave them: a
    sub-window TRF couldn't finish (`timed_out`) is a row of its window's BED
    marked timed_out, and the excluded group holds it, the homopolymer runs and
    the long-period tracts.
    """
    s = targets.sample(tmp_path, monkeypatch, "ACGT" * 250)
    for model in MODELS:
        s.write_tracts("tandem_repeats", tracts if model == "M1" else [], model=model, contig=CONTIG)
    s.write("homopolymers", bed(homopolymers), contig=CONTIG)
    s.write("long_periods_window", bed(timed_out, "timed_out"), contig=CONTIG, start=0)
    s.write("excluded", bed(sorted([*homopolymers, *long_periods, *timed_out])), contig=CONTIG)
    s.write("unresolved", bed(unresolved), contig=CONTIG)
    return s


def causes(*args):
    """tract_drop_causes.py's table, run from the sample's checkout: {(model, measure): (tracts, bp)}."""
    out = io.StringIO()
    argv, path, sys.argv = sys.argv, sys.path[:], ["tract_drop_causes.py", SPECIES, SAMPLE, *args]
    try:
        with contextlib.redirect_stdout(out):
            runpy.run_path(os.path.join(REPO_ROOT, "scripts", "tract_drop_causes.py"), run_name="__main__")
    finally:
        sys.argv, sys.path[:] = argv, path
    header, *rows = [line.split("\t") for line in out.getvalue().splitlines()]
    assert header == ["model", "measure", "tracts", "bp", "species", "sample"]
    assert all(row[4:] == [SPECIES, SAMPLE] for row in rows)
    return {(row[0], row[1]): (int(row[2]), int(row[3])) for row in rows}


def test_each_dropped_tract_is_counted_once_under_its_first_cause(tmp_path, monkeypatch):
    sample(tmp_path, monkeypatch)
    table = causes()
    assert [table["M1", measure] for measure in
            ("all", "dropped_unresolved", "dropped_homopolymer", "dropped_long_period", "option_1")] == \
        [(6, 240), (1, 40), (2, 80), (2, 80), (1, 40)]
    assert table["M0", "all"] == table["M2", "all"] == (0, 0)


def test_option_2_drops_tracts_only_for_unresolved_regions_and_homopolymer_runs(tmp_path, monkeypatch):
    sample(tmp_path, monkeypatch)
    assert causes()["M1", "option_2"] == (3, 120)   # the first, fourth and last


def test_option_3_keeps_tracts_whose_excluded_bases_reach_in_from_an_end_up_to_the_edge(tmp_path, monkeypatch):
    sample(tmp_path, monkeypatch)
    table = causes("--edges", "1,2,5")
    # 2 bases of a homopolymer run at the third tract's end, 5 of a long-period
    # tract at the last one's; the run inside the fifth is never an edge.
    assert [table["M1", f"option_3_edge_{edge}"] for edge in (1, 2, 5)] == [(1, 40), (2, 80), (3, 120)]
    # Options 2 and 3 together: long periods ignored, homopolymer runs at the edge.
    assert [table["M1", f"option_2_3_edge_{edge}"] for edge in (1, 2, 5)] == [(3, 120), (4, 160), (4, 160)]


def test_option_3_counts_the_excluded_bases_at_both_ends_together(tmp_path, monkeypatch):
    sample(tmp_path, monkeypatch, tracts=[(100, 140, "AC", 2)], homopolymers=[(90, 101), (139, 150)], long_periods=[])
    table = causes("--edges", "1,2")
    assert [table["M1", "option_3_edge_1"], table["M1", "option_3_edge_2"]] == [(0, 0), (1, 40)]


def test_a_tract_inside_an_excluded_region_is_never_kept_by_option_3(tmp_path, monkeypatch):
    # 8 bp, fewer than the edge, but covered from end to end.
    sample(tmp_path, monkeypatch, tracts=[(400, 408, "AC", 2)], homopolymers=[], long_periods=[(380, 600)])
    table = causes("--edges", "10")
    assert [table["M1", "option_3_edge_10"], table["M1", "option_2_3_edge_10"]] == [(0, 0), (1, 8)]


def test_a_sub_window_trf_couldnt_finish_still_drops_its_tracts_under_option_2(tmp_path, monkeypatch):
    # Not a long-period tract TRF found, but most likely a satellite array.
    sample(tmp_path, monkeypatch, tracts=[(100, 140, "AC", 2), (500, 540, "AC", 2)], homopolymers=[],
           long_periods=[(90, 150)], unresolved=[], timed_out=[(400, 600)])
    table = causes()
    assert [table["M1", measure] for measure in ("dropped_timed_out", "dropped_long_period", "option_2")] == \
        [(1, 40), (1, 40), (1, 40)]


def test_an_empty_repeat_bed_has_no_tracts(tmp_path, monkeypatch):
    s = sample(tmp_path, monkeypatch)
    open(s.path("tandem_repeats", contig=CONTIG, model="M0"), "w").close()
    assert causes()["M0", "all"] == (0, 0)


def test_option_1_is_what_workflow_4_measures_today(tmp_path, monkeypatch):
    s = sample(tmp_path, monkeypatch)
    s.run("repetitive_indel_call", f"TRTracts_{SPECIES}_{SAMPLE}")
    measured = {row["model"]: (int(row["tracts_measured"]), int(row["bp_measured"]))
                for row in s.table("repetitive_tracts_summary")}
    table = causes()
    assert {model: table[model, "option_1"] for model in MODELS} == measured
