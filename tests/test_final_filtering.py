"""
The final filtering stage (ticket 24, Q7): one definition of a final call,
applied by scripts/final_filtering.py to a sample's calls and opportunities
alike, so a filter never lowers the rate by what it removes from the calls
alone. SMRT Cells whose base qualities top out at Q40 are dropped (Q12, Q17); a
filter on a tag the sample lacks is skipped and recorded (Q15); every row
carries its platform.
"""

import os
import subprocess
import sys

from scans import rows_of

SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts", "final_filtering.py")


def table(path, header, rows):
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(map(str, row)) + "\n" for row in rows))
    return path


def calls(tmp_path, rows):
    """Calls: (readID, smrt_cell, kmerQ, kmerQ_pass, depth_bin, subread_support, haploid_check)."""
    return table(tmp_path / "calls.tsv", ["#contig", "readID", "smrt_cell", "kmerQ", "kmerQ_pass", "depth_bin",
                                          "subread_support", "haploid_check"], [("c1", *row) for row in rows])


def opportunities(tmp_path, rows):
    """Opportunities: (smrt_cell, kmerQ_pass, depth_bin, subread_support, count)."""
    return table(tmp_path / "opportunities.tsv", ["contig", "smrt_cell", "kmerQ_pass", "depth_bin",
                                                  "subread_support", "count"], [("c1", *row) for row in rows])


def smrt_cells(tmp_path, best):
    """The sample's SMRT Cells, as precompute writes them: {cell: its reads' highest base quality}."""
    return table(tmp_path / "smrt_cells.tsv", ["smrt_cell", "reads", "max_quality"],
                 [(cell, 100, quality) for cell, quality in best.items()])


def final(tmp_path, call_table, opportunity_table, *options, platform="Revio", cells=None):
    """(final calls, final opportunities, report), as row dicts."""
    out = [tmp_path / name for name in ("final_calls.tsv", "final_opportunities.tsv", "report.tsv")]
    cells = smrt_cells(tmp_path, cells or {"m1": 93})
    subprocess.run([sys.executable, SCRIPT, str(call_table), str(opportunity_table), *map(str, out),
                    "--platform", platform, "--smrt-cells", str(cells), "--capped-quality", "40", *options],
                   check=True, capture_output=True)
    return [rows_of(path.read_text()) for path in out]


def rate(final_calls, final_opportunities):
    return len(final_calls) / sum(int(row["count"]) for row in final_opportunities)


def test_a_filter_on_a_feature_spread_alike_leaves_the_rate_unchanged(tmp_path):
    # Calls and opportunities split 3:1 between two values of a feature.
    call_rows = [(f"r{i}", "m1", 93, "pass", "low" if i % 4 else "high", ">=1", "True") for i in range(40)]
    opportunity_rows = [("m1", "pass", "low", ">=1", 3000), ("m1", "pass", "high", ">=1", 1000)]
    before = final(tmp_path, calls(tmp_path, call_rows), opportunities(tmp_path, opportunity_rows))
    after = final(tmp_path, calls(tmp_path, call_rows), opportunities(tmp_path, opportunity_rows),
                  "--keep", "depth_bin=low")
    assert len(after[0]) == 30
    assert rate(*after[:2]) == rate(*before[:2]) == 0.01


def test_a_filter_keeps_the_rows_with_any_of_its_values(tmp_path):
    call_rows = [("a", "m1", 93, "pass", "low", ">=1", "True"), ("b", "m1", 93, "pass", "mid", ">=1", "True"),
                 ("c", "m1", 93, "pass", "high", ">=1", "True")]
    opportunity_rows = [("m1", "pass", value, ">=1", 10) for value in ("low", "mid", "high")]
    final_calls, final_opportunities, _ = final(tmp_path, calls(tmp_path, call_rows),
                                                opportunities(tmp_path, opportunity_rows), "--keep", "depth_bin=low,high")
    assert [row["readID"] for row in final_calls] == ["a", "c"]
    assert sorted(row["depth_bin"] for row in final_opportunities) == ["high", "low"]


def test_smrt_cells_capped_at_q40_are_dropped_from_both_tables_and_listed(tmp_path):
    # Capped or not is the cell's reads' highest base quality, the same for
    # every table of the sample, whatever its calls show: m_silent has none.
    call_rows = [("a", "m_capped", 93, "pass", "low", ">=1", "True"), ("b", "m_capped", 31, "fail", "low", ">=1", "NA"),
                 ("c", "m_full", 93, "pass", "low", ">=1", "True"), ("d", "m_full", 40, "fail", "low", ">=1", "True")]
    opportunity_rows = [("m_capped", "fail", "low", ">=1", 500), ("m_full", "pass", "low", ">=1", 400),
                        ("m_full", "fail", "low", ">=1", 100), ("m_silent", "fail", "low", ">=1", 50)]
    final_calls, final_opportunities, report = final(tmp_path, calls(tmp_path, call_rows),
                                                     opportunities(tmp_path, opportunity_rows),
                                                     cells={"m_capped": 40, "m_full": 93, "m_silent": 39})
    assert [row["readID"] for row in final_calls] == ["c", "d"]
    assert {row["smrt_cell"] for row in final_opportunities} == {"m_full"}
    step = next(row for row in report if row["filter"] == "capped_smrt_cells")
    assert (step["calls_before"], step["calls_after"], step["opportunities_before"], step["opportunities_after"],
            step["note"]) == ("4", "2", "1050", "500", "m_capped,m_silent")


def test_the_smrt_cells_highest_base_quality_comes_from_its_reads(tmp_path):
    import make_fixtures as mf
    reference = mf.build_reference()
    capped = {p: 40 for p in range(1000)}
    reads = [mf.make_read("m_capped/1/ccs", reference, end=1000, qualities=capped),
             mf.make_read("m_capped/2/ccs", reference, start=10, end=1000, qualities={**capped, 500: 30}),
             mf.make_read("m_full/1/ccs", reference, end=1000, qualities={500: 30})]
    paths = mf.write_dataset(str(tmp_path), reference, reads)
    out = tmp_path / "smrt_cells.tsv"
    subprocess.run([sys.executable, SCRIPT.replace("final_filtering.py", "smrt_cell_qualities.py"),
                    paths["bam"], str(out), "--reads", "10"], check=True, capture_output=True)
    assert {row["smrt_cell"]: (row["reads"], row["max_quality"]) for row in rows_of(out.read_text())} == {
        "m_capped": ("2", "40"), "m_full": ("1", "93")}


def test_a_filter_on_a_tag_the_sample_lacks_is_skipped_and_recorded(tmp_path):
    # Sequel II: no subread tags, so subread_support is NA throughout.
    call_rows = [("a", "m1", 93, "pass", "low", "NA", "True"), ("b", "m1", 93, "pass", "low", "NA", "True")]
    opportunity_rows = [("m1", "pass", "low", "NA", 200)]
    final_calls, final_opportunities, report = final(
        tmp_path, calls(tmp_path, call_rows), opportunities(tmp_path, opportunity_rows),
        "--keep", "subread_support=>=1", "--keep", "kmerQ_pass=pass", platform="Sequel II")
    assert len(final_calls) == 2 and sum(int(row["count"]) for row in final_opportunities) == 200
    step = next(row for row in report if row["filter"] == "subread_support=>=1")
    assert (step["calls_after"], step["note"]) == ("2", "skipped: NA throughout the sample")


def test_every_row_carries_its_platform(tmp_path):
    call_rows = [("a", "m1", 93, "pass", "low", "NA", "True")]
    final_calls, final_opportunities, _ = final(tmp_path, calls(tmp_path, call_rows),
                                                opportunities(tmp_path, [("m1", "pass", "low", "NA", 5)]),
                                                platform="Sequel II")
    assert {row["platform"] for row in final_calls + final_opportunities} == {"Sequel II"}


def test_the_report_lists_every_step_with_rows_before_and_after(tmp_path):
    call_rows = [("a", "m1", 93, "pass", "low", ">=1", "True"), ("b", "m1", 93, "pass", "high", ">=1", "False"),
                 ("c", "m1", 93, "pass", "high", ">=1", "NA"), ("d", "m1", 10, "fail", "high", ">=1", "NA")]
    opportunity_rows = [("m1", "pass", "low", ">=1", 30), ("m1", "pass", "high", ">=1", 60),
                        ("m1", "fail", "high", ">=1", 10)]
    _, _, report = final(tmp_path, calls(tmp_path, call_rows), opportunities(tmp_path, opportunity_rows),
                         "--keep", "kmerQ_pass=pass", "--keep", "depth_bin=high")
    columns = ("filter", "calls_before", "calls_after", "opportunities_before", "opportunities_after",
               "haploid_check_missing")
    assert [tuple(row[column] for column in columns) for row in report] == [
        ("input", "4", "4", "100", "100", "2"),
        ("capped_smrt_cells", "4", "4", "100", "100", "2"),
        ("kmerQ_pass=pass", "4", "3", "100", "90", "1"),
        ("depth_bin=high", "3", "2", "90", "60", "1"),
    ]


# -- big tables (ticket 32) ---------------------------------------------------


def test_the_result_is_the_same_whatever_the_chunk_size(tmp_path):
    # Every step falls across chunk edges: capped cells, a skipped filter and
    # two kept filters, on tables read a row or a few rows at a time.
    call_rows = [(f"r{i}", "m_capped" if i % 5 == 0 else "m1", 93 if i % 3 else 40, "pass" if i % 3 else "fail",
                  "low" if i % 4 else "high", "NA", "NA" if i % 7 == 0 else "True") for i in range(23)]
    opportunity_rows = [("m_capped" if i % 5 == 0 else "m1", "pass" if i % 3 else "fail", "low" if i % 4 else "high",
                         "NA", 10 + i) for i in range(17)]
    options = ["--keep", "subread_support=>=1", "--keep", "kmerQ_pass=pass", "--keep", "depth_bin=low"]
    cells = {"m_capped": 40, "m1": 93}
    whole = final(tmp_path, calls(tmp_path, call_rows), opportunities(tmp_path, opportunity_rows), *options, cells=cells)
    for rows in (1, 2, 5):
        chunked = final(tmp_path, calls(tmp_path, call_rows), opportunities(tmp_path, opportunity_rows), *options,
                        "--chunk-rows", str(rows), cells=cells)
        assert chunked == whole, rows


def test_gzipped_tables_are_read_and_written(tmp_path):
    import gzip
    call_rows = [("a", "m1", 93, "pass", "low", ">=1", "True"), ("b", "m1", 93, "pass", "high", ">=1", "True")]
    opportunity_rows = [("m1", "pass", "low", ">=1", 30), ("m1", "pass", "high", ">=1", 10)]
    plain = opportunities(tmp_path, opportunity_rows)
    zipped = tmp_path / "opportunities.tsv.gz"
    zipped.write_bytes(gzip.compress(plain.read_bytes()))
    out = [tmp_path / name for name in ("final_calls.tsv", "final_opportunities.tsv.gz", "report.tsv")]
    subprocess.run([sys.executable, SCRIPT, str(calls(tmp_path, call_rows)), str(zipped), *map(str, out),
                    "--platform", "Revio", "--smrt-cells", str(smrt_cells(tmp_path, {"m1": 93})),
                    "--capped-quality", "40", "--keep", "depth_bin=low"], check=True, capture_output=True)
    final_opportunities = rows_of(gzip.decompress(out[1].read_bytes()).decode())
    assert [(row["depth_bin"], row["count"], row["platform"]) for row in final_opportunities] == [("low", "30", "Revio")]
