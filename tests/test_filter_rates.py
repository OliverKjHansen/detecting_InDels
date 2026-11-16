"""
Rates by filter metric, and a funnel, for one sample (ticket 56): from the
candidate tables, before final filtering, the unique calls, the non-unique
calls and the opportunities of each value of each metric, and after each
step of final filtering and of the funnel's own steps, with the same filter
on both tables.
"""

import os
import subprocess
import sys

import duckdb
import pytest

from scans import rows_of
from test_final_filtering import final, smrt_cells, table

SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts", "filter_rates.py")
CALL_COLUMNS = ["#contig", "indel_type", "unique", "indel_string", "readID", "smrt_cell", "rq_bin", "haploid_check"]
OPPORTUNITY_COLUMNS = ["contig", "indel_type", "size", "context", "smrt_cell", "rq_bin", "count"]


def calls(tmp_path, rows):
    """Calls: (indel_type, unique, indel_string, smrt_cell, rq_bin, haploid_check)."""
    return table(tmp_path / "calls.tsv", CALL_COLUMNS,
                 [("c1", kind, seen, bases, f"r{i}", *rest) for i, (kind, seen, bases, *rest) in enumerate(rows)])


def opportunities(tmp_path, rows):
    """Opportunities: (indel_type, size, smrt_cell, rq_bin, count)."""
    return table(tmp_path / "opportunities.tsv", OPPORTUNITY_COLUMNS,
                 [("c1", kind, size, "A" * size, *rest) for kind, size, *rest in rows])


def rates(tmp_path, call_table, opportunity_table, *options, cells=None):
    """(rows by filter, funnel rows), as row dicts."""
    out = [tmp_path / "by_filter.tsv", tmp_path / "funnel.tsv"]
    cells = smrt_cells(tmp_path, cells or {"m1": 93})
    subprocess.run([sys.executable, SCRIPT, str(call_table), str(opportunity_table), *map(str, out),
                    "--species", "ape", "--sample", "s1", "--platform", "Revio",
                    "--smrt-cells", str(cells), "--capped-quality", "40", *options],
                   check=True, capture_output=True)
    return [rows_of(path.read_text()) for path in out]


def counts(row):
    return int(row["unique_calls"]), int(row["non_unique_calls"]), int(row["opportunities"]), float(row["rate"])


def given(unique_calls, non_unique_calls, opportunities_):
    """The counts of a row, and its rate as written (6 significant digits)."""
    return pytest.approx((unique_calls, non_unique_calls, opportunities_, unique_calls / opportunities_), rel=1e-5)


def test_each_value_of_a_metric_gets_its_calls_its_opportunities_and_their_rate(tmp_path):
    call_rows = ([("deletion", "unique", "A", "m1", ">=0.9999", "True")] * 3
                 + [("deletion", "unique", "A", "m1", "0.99-0.999", "True")]
                 + [("deletion", "non-unique", "A", "m1", ">=0.9999", "True")] * 2)
    opportunity_rows = [("deletion", 1, "m1", ">=0.9999", 1000), ("deletion", 1, "m1", "0.99-0.999", 100)]
    by_filter, _ = rates(tmp_path, calls(tmp_path, call_rows), opportunities(tmp_path, opportunity_rows))
    rq = {row["value"]: row for row in by_filter if row["metric"] == "rq_bin"}
    assert set(rq) == {">=0.9999", "0.99-0.999"}
    assert counts(rq[">=0.9999"]) == given(3, 2, 1000)
    assert counts(rq["0.99-0.999"]) == given(1, 0, 100)
    assert {(row["indel_type"], row["size"], row["species"], row["sample"]) for row in rq.values()} == {
        ("deletion", "1", "ape", "s1")}


def test_a_calls_only_check_divides_its_calls_by_all_the_opportunities(tmp_path):
    # The opportunities have no haploid check: each of its values is a share
    # of the calls over the same opportunities.
    call_rows = ([("deletion", "unique", "A", "m1", ">=0.9999", "True")] * 2
                 + [("deletion", "unique", "A", "m1", ">=0.9999", "False")])
    opportunity_rows = [("deletion", 1, "m1", ">=0.9999", 1000), ("deletion", 1, "m1", "0.99-0.999", 100)]
    by_filter, _ = rates(tmp_path, calls(tmp_path, call_rows), opportunities(tmp_path, opportunity_rows),
                         "--calls-only", "haploid_check=True,NA")
    check = {row["value"]: row for row in by_filter if row["metric"] == "haploid_check"}
    assert counts(check["True"]) == given(2, 0, 1100)
    assert counts(check["False"]) == given(1, 0, 1100)
    assert {row["applies_to"] for row in check.values()} == {"calls"}
    assert {row["applies_to"] for row in by_filter if row["metric"] == "rq_bin"} == {"both"}


def test_the_funnel_counts_after_each_step_with_final_filterings_steps_first(tmp_path):
    # m2's base qualities top out at Q40: final filtering drops its reads
    # from both tables first; then the haploid check, on the calls only; then
    # the funnel's own step, on both tables.
    call_rows = ([("deletion", "unique", "A", "m1", ">=0.9999", "True")] * 2
                 + [("deletion", "unique", "A", "m1", ">=0.9999", "False"),
                    ("deletion", "unique", "A", "m1", "0.99-0.999", "True"),
                    ("deletion", "unique", "A", "m2", ">=0.9999", "True"),
                    ("deletion", "non-unique", "A", "m1", ">=0.9999", "True")])
    opportunity_rows = [("deletion", 1, "m1", ">=0.9999", 1000), ("deletion", 1, "m1", "0.99-0.999", 100),
                        ("deletion", 1, "m2", ">=0.9999", 500)]
    _, funnel = rates(tmp_path, calls(tmp_path, call_rows), opportunities(tmp_path, opportunity_rows),
                      "--calls-only", "haploid_check=True,NA", "--step", "rq_bin=>=0.9999",
                      cells={"m1": 93, "m2": 40})
    assert [(row["step"], row["filter"], row["applies_to"]) for row in funnel] == [
        ("1", "input", "both"), ("2", "capped_smrt_cells", "both"), ("3", "haploid_check=True,NA", "calls"),
        ("4", "rq_bin=>=0.9999", "both")]
    assert [counts(row) for row in funnel] == [given(5, 1, 1600), given(4, 1, 1100), given(3, 1, 1100),
                                              given(2, 1, 1000)]


def test_each_indel_type_and_size_gets_rows_of_its_own(tmp_path):
    # A call's size is its sequence's length.
    call_rows = [("deletion", "unique", "A", "m1", ">=0.9999", "True"),
                 ("deletion", "unique", "AC", "m1", ">=0.9999", "True"),
                 ("deletion", "unique", "AC", "m1", ">=0.9999", "True"),
                 ("insertion", "unique", "ACG", "m1", ">=0.9999", "True")]
    opportunity_rows = [("deletion", 1, "m1", ">=0.9999", 100), ("deletion", 2, "m1", ">=0.9999", 200),
                        ("insertion", 3, "m1", ">=0.9999", 300), ("insertion", 1, "m1", ">=0.9999", 400)]
    by_filter, funnel = rates(tmp_path, calls(tmp_path, call_rows), opportunities(tmp_path, opportunity_rows))
    expected = {("deletion", "1"): given(1, 0, 100), ("deletion", "2"): given(2, 0, 200),
                ("insertion", "1"): (0, 0, 400, 0.0), ("insertion", "3"): given(1, 0, 300)}
    for rows in ([row for row in by_filter if row["metric"] == "rq_bin"],
                 [row for row in funnel if row["filter"] == "input"]):
        assert {(row["indel_type"], row["size"]): counts(row) for row in rows} == expected


def test_the_funnel_after_final_filterings_steps_counts_what_final_filtering_keeps(tmp_path):
    call_rows = ([("deletion", "unique", "A", "m1", ">=0.9999", "True")] * 3
                 + [("deletion", "unique", "A", "m1", ">=0.9999", "False"),
                    ("deletion", "unique", "A", "m1", "0.99-0.999", "NA"),
                    ("deletion", "unique", "A", "m2", ">=0.9999", "True"),
                    ("deletion", "non-unique", "A", "m1", ">=0.9999", "True"),
                    ("insertion", "unique", "AC", "m1", ">=0.9999", "NA")])
    opportunity_rows = [("deletion", 1, "m1", ">=0.9999", 1000), ("deletion", 1, "m1", "0.99-0.999", 100),
                        ("deletion", 1, "m2", ">=0.9999", 500), ("insertion", 2, "m1", ">=0.9999", 700)]
    options = ["--keep", "rq_bin=>=0.9999", "--calls-only", "haploid_check=True,NA"]
    cells = {"m1": 93, "m2": 40}
    call_table, opportunity_table = calls(tmp_path, call_rows), opportunities(tmp_path, opportunity_rows)
    _, funnel = rates(tmp_path, call_table, opportunity_table, *options, "--step", "smrt_cell=m1", cells=cells)
    final_calls, final_opportunities, _ = final(tmp_path, call_table, opportunity_table, *options, cells=cells)

    def kept(kind, size):
        return (sum(row["indel_type"] == kind and len(row["indel_string"]) == size and row["unique"] == "unique"
                    for row in final_calls),
                sum(row["indel_type"] == kind and len(row["indel_string"]) == size and row["unique"] == "non-unique"
                    for row in final_calls),
                sum(int(row["count"]) for row in final_opportunities
                    if row["indel_type"] == kind and int(row["size"]) == size))

    last_final_step = [row for row in funnel if row["filter"] == "haploid_check=True,NA"]
    assert {(row["indel_type"], int(row["size"])): counts(row)[:3] for row in last_final_step} == {
        ("deletion", 1): kept("deletion", 1), ("insertion", 2): kept("insertion", 2)}
    assert kept("deletion", 1) == (3, 1, 1000)


def test_a_funnel_step_on_a_tag_the_sample_lacks_is_skipped_and_noted(tmp_path):
    # A Sequel II sample has no rq: the step would remove everything.
    call_rows = [("deletion", "unique", "A", "m1", "NA", "True")] * 2
    opportunity_rows = [("deletion", 1, "m1", "NA", 1000)]
    _, funnel = rates(tmp_path, calls(tmp_path, call_rows), opportunities(tmp_path, opportunity_rows),
                      "--step", "rq_bin=>=0.9999")
    step = [row for row in funnel if row["filter"] == "rq_bin=>=0.9999"]
    assert [counts(row) for row in step] == [given(2, 0, 1000)]
    assert step[0]["note"] == "skipped: NA throughout the sample"


def test_a_parquet_cube_is_read_as_workflow_3_writes_it(tmp_path):
    # Typed columns: size an integer, count a BIGINT.
    call_rows = [("deletion", "unique", "A", "m1", ">=0.9999", "True"),
                 ("insertion", "unique", "AC", "m1", ">=0.9999", "True")]
    tsv = opportunities(tmp_path, [("deletion", 1, "m1", ">=0.9999", 100), ("insertion", 2, "m1", ">=0.9999", 50)])
    cube = tmp_path / "opportunities.parquet"
    duckdb.execute(f"COPY (SELECT * REPLACE (CAST(size AS INTEGER) AS size, CAST(count AS BIGINT) AS count) "
                   f"FROM read_csv('{tsv}', delim='\\t', header=true, all_varchar=true)) TO '{cube}' (FORMAT parquet)")
    by_filter, funnel = rates(tmp_path, calls(tmp_path, call_rows), cube)
    expected = {("deletion", "1"): given(1, 0, 100), ("insertion", "2"): given(1, 0, 50)}
    assert {(row["indel_type"], row["size"]): counts(row) for row in by_filter if row["metric"] == "rq_bin"} == expected
    assert {(row["indel_type"], row["size"]): counts(row) for row in funnel if row["filter"] == "input"} == expected


def test_the_metrics_are_the_keys_both_tables_carry(tmp_path):
    # Not the contig, the context, the type or the size; not the calls' own
    # columns (readID, haploid_check) unless a calls-only check names one.
    call_rows = [("deletion", "unique", "A", "m1", ">=0.9999", "True")]
    by_filter, _ = rates(tmp_path, calls(tmp_path, call_rows),
                         opportunities(tmp_path, [("deletion", 1, "m1", ">=0.9999", 100)]))
    assert {row["metric"] for row in by_filter} == {"smrt_cell", "rq_bin"}


def test_a_sample_without_calls_gets_its_opportunities_and_full_headers(tmp_path):
    by_filter, funnel = rates(tmp_path, calls(tmp_path, []),
                              opportunities(tmp_path, [("deletion", 1, "m1", ">=0.9999", 100)]))
    assert [counts(row) for row in funnel if row["filter"] == "input"] == [(0, 0, 100, 0.0)]
    by_filter_path, funnel_path = tmp_path / "by_filter.tsv", tmp_path / "funnel.tsv"
    assert by_filter_path.read_text().split("\n")[0].split("\t") == [
        "species", "sample", "platform", "metric", "value", "applies_to", "indel_type", "size",
        "unique_calls", "non_unique_calls", "opportunities", "rate"]
    assert funnel_path.read_text().split("\n")[0].split("\t") == [
        "species", "sample", "platform", "step", "filter", "applies_to", "note", "indel_type", "size",
        "unique_calls", "non_unique_calls", "opportunities", "rate"]
    assert not list(tmp_path.glob("*.tmp"))


def test_two_empty_tables_give_headers_only(tmp_path):
    by_filter, funnel = rates(tmp_path, calls(tmp_path, []), opportunities(tmp_path, []))
    assert by_filter == funnel == []
    assert (tmp_path / "funnel.tsv").read_text().startswith("species\tsample\tplatform\tstep\t")


def test_the_funnel_is_split_by_the_by_columns_and_its_strata_add_up_to_the_sample(tmp_path):
    # Per tissue and methylation bin, so the rate of each can be plotted after
    # any step; every stratum gets the same steps.
    call_table = table(tmp_path / "calls.tsv",
                       ["#contig", "indel_type", "unique", "indel_string", "smrt_cell", "rq_bin", "tissue", "methylation_bin"],
                       [("c1", "deletion", "unique", "A", "m1", ">=0.9999", "testis", "-12"),
                        ("c1", "deletion", "unique", "A", "m1", ">=0.9999", "testis", "-12"),
                        ("c1", "deletion", "unique", "A", "m1", "0.99-0.999", "testis", "-12"),
                        ("c1", "deletion", "unique", "A", "m1", ">=0.9999", "blood", "4")])
    opportunity_table = table(tmp_path / "opportunities.tsv",
                              ["contig", "indel_type", "size", "context", "smrt_cell", "rq_bin", "tissue", "methylation_bin",
                               "count"],
                              [("c1", "deletion", 1, "A", "m1", ">=0.9999", "testis", "-12", 1000),
                               ("c1", "deletion", 1, "A", "m1", "0.99-0.999", "testis", "-12", 1000),
                               ("c1", "deletion", 1, "A", "m1", ">=0.9999", "blood", "4", 500)])
    _, funnel = rates(tmp_path, call_table, opportunity_table, "--step", "rq_bin=>=0.9999",
                      "--by", "tissue", "--by", "methylation_bin")
    strata = {(row["filter"], row["tissue"], row["methylation_bin"]): counts(row) for row in funnel}
    assert strata == {("input", "testis", "-12"): given(3, 0, 2000), ("input", "blood", "4"): given(1, 0, 500),
                      ("capped_smrt_cells", "testis", "-12"): given(3, 0, 2000),
                      ("capped_smrt_cells", "blood", "4"): given(1, 0, 500),
                      ("rq_bin=>=0.9999", "testis", "-12"): given(2, 0, 1000),
                      ("rq_bin=>=0.9999", "blood", "4"): given(1, 0, 500)}
    header = (tmp_path / "funnel.tsv").read_text().split("\n")[0].split("\t")
    assert header[header.index("note") + 1:header.index("indel_type")] == ["tissue", "methylation_bin"]


def test_a_by_column_a_table_lacks_stops_the_run_with_its_name(tmp_path):
    call_table = calls(tmp_path, [("deletion", "unique", "A", "m1", ">=0.9999", "True")])
    opportunity_table = opportunities(tmp_path, [("deletion", 1, "m1", ">=0.9999", 100)])
    with pytest.raises(subprocess.CalledProcessError) as failure:
        rates(tmp_path, call_table, opportunity_table, "--by", "tissue")
    assert b"--by tissue: the calls have no column tissue" in failure.value.stderr
