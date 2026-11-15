"""
Rates by filter metric, and a funnel, for one sample (ticket 56): small
tables of workflow 3's candidate calls and opportunities, before final
filtering, from which a rate can be plotted without the opportunity cube.

    python filter_rates.py <calls.tsv> <opportunities.parquet> <by_filter.tsv> <funnel.tsv> \\
        --species S --sample X --platform <Revio | "Sequel II"> --smrt-cells <smrt_cells.tsv> --capped-quality <Q> \\
        [--keep COLUMN=VALUE[,VALUE...] ...] [--calls-only COLUMN=VALUE[,VALUE...] ...] [--step COLUMN=VALUE[,VALUE...] ...]

Every row counts, for one InDel type and size, the unique calls, the
non-unique calls and the opportunities (summed counts), and gives the rate:
unique calls per opportunity. An InDel's size is its sequence's length.

by_filter.tsv: one row per metric, value, type and size. A metric is a key
both tables carry, a value of it kept in both, so the rows of one metric add
up to the sample's; the rate of a threshold is that of its values added up.
The --calls-only checks are metrics of the calls alone (applies_to "calls"):
the opportunities can't take them, so each of their values divides its calls
by all the opportunities of the type and size.

funnel.tsv: one row per step, type and size: the counts after the step and
every step before it. The steps are final filtering's, in its order and with
its rules (final_filtering.py: the capped SMRT Cells, each --keep, each
--calls-only), then each --step, on both tables like a --keep. A step on a
tag the sample lacks is skipped, as there, and its note says so. The step
after the last --calls-only is final filtering's output; the last step is the
most filtered rate.
"""
import argparse
import os

import pandas

from duckdb_tables import connect, identifier
from final_filtering import Table, capped_smrt_cells, final_steps, keep_step, kept

# Columns both tables have that aren't filter metrics.
NOT_METRICS = {"contig", "indel_type", "size", "context", "count", "species", "sample", "platform"}
GROUPS = ["indel_type", "size"]
SAMPLE = ["species", "sample", "platform"]
COUNTS = ["indel_type", "size", "unique_calls", "non_unique_calls", "opportunities", "rate"]
BY_FILTER_COLUMNS = SAMPLE + ["metric", "value", "applies_to"] + COUNTS
FUNNEL_COLUMNS = SAMPLE + ["step", "filter", "applies_to", "note"] + COUNTS


def with_size(table):
    """
    A table's rows as an SQL table, with each InDel's size as an integer: the
    opportunities' own column (text in a TSV), the length of a call's sequence.
    """
    if table.count_column is None:
        return f"(SELECT *, CAST(length({identifier('indel_string')}) AS INTEGER) AS size FROM {table.relation})"
    return f"(SELECT * REPLACE (CAST(size AS INTEGER) AS size) FROM {table.relation})"


def by_metric(con, relation, metrics, size, groups):
    """
    {(metric, value, *groups): size} over the SQL table `relation`, in one
    pass: each metric's values, with the `groups` columns; and the whole
    table's, under (None, None, *groups).
    """
    which = " ".join(f"WHEN grouping({identifier(m)}) = 0 THEN '{m}'" for m in metrics)
    value = " ".join(f"WHEN grouping({identifier(m)}) = 0 THEN CAST({identifier(m)} AS VARCHAR)" for m in metrics)
    keys = ", ".join(map(identifier, groups))
    sets = ", ".join([f"({keys})"] + [f"({keys}, {identifier(m)})" for m in metrics])
    which, value = (f"CASE {which} END", f"CASE {value} END") if metrics else ("NULL", "NULL")
    rows = con.execute(f"SELECT {which}, {value}, {keys}, {size} FROM {relation} "
                       f"GROUP BY GROUPING SETS ({sets})").fetchall()
    return {tuple(row[:-1]): row[-1] or 0 for row in rows}


def rate(unique_calls, opportunities):
    return unique_calls / opportunities if opportunities else float("nan")


def by_filter_rows(con, calls, opportunities, metrics, calls_only):
    """The rows of by_filter.tsv, as dicts: the metrics of both tables, then the calls-only ones."""
    called = by_metric(con, with_size(calls), metrics + calls_only, "count(*)", GROUPS + ["unique"])
    counted = by_metric(con, with_size(opportunities), metrics, opportunities.size(), GROUPS)
    keys = sorted({key[:-1] for key in called if key[0] is not None}
                  | {key for key in counted if key[0] is not None}, key=str)
    rows = []
    for metric, value, indel_type, size in keys:
        alone = metric in calls_only
        unique_calls = called.get((metric, value, indel_type, size, "unique"), 0)
        opportunities_ = counted.get((None, None, indel_type, size) if alone else (metric, value, indel_type, size), 0)
        rows.append({"metric": metric, "value": value, "applies_to": "calls" if alone else "both",
                     "indel_type": indel_type, "size": size, "unique_calls": unique_calls,
                     "non_unique_calls": called.get((metric, value, indel_type, size, "non-unique"), 0),
                     "opportunities": opportunities_, "rate": rate(unique_calls, opportunities_)})
    return rows


def after_each_step(con, table, relation, steps, groups):
    """{groups: [size after each step]} of `table`, read as the SQL table `relation`, in one pass."""
    sums = [f"{table.size()} FILTER (WHERE {kept(steps[:i + 1], table.name)})" for i in range(len(steps))]
    keys = ", ".join(map(identifier, groups))
    rows = con.execute(f"SELECT {keys}, {', '.join(sums)} FROM {relation} GROUP BY {keys}").fetchall()
    return {tuple(row[:len(groups)]): [count or 0 for count in row[len(groups):]] for row in rows}


def funnel_rows(con, calls, opportunities, steps):
    """The rows of funnel.tsv, as dicts."""
    called = after_each_step(con, calls, with_size(calls), steps, GROUPS + ["unique"])
    counted = after_each_step(con, opportunities, with_size(opportunities), steps, GROUPS)
    kinds = sorted(set(counted) | {key[:-1] for key in called}, key=str)
    none = [0] * len(steps)
    rows = []
    for i, step in enumerate(steps):
        for indel_type, size in kinds:
            unique_calls = called.get((indel_type, size, "unique"), none)[i]
            opportunities_ = counted.get((indel_type, size), none)[i]
            rows.append({"step": i + 1, "filter": step.name,
                         "applies_to": "both" if len(step.tables) == 2 else ",".join(step.tables), "note": step.note,
                         "indel_type": indel_type, "size": size, "unique_calls": unique_calls,
                         "non_unique_calls": called.get((indel_type, size, "non-unique"), none)[i],
                         "opportunities": opportunities_, "rate": rate(unique_calls, opportunities_)})
    return rows


def write(rows, columns, args, path):
    """Writes the rows, with the sample's columns, to a temporary file first, renamed once complete."""
    sample = {"species": args.species, "sample": args.sample, "platform": args.platform}
    table = pandas.DataFrame([{**sample, **row} for row in rows], columns=columns)
    table.to_csv(path + ".tmp", sep="\t", index=False, float_format="%.6g", na_rep="NA")
    os.replace(path + ".tmp", path)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("calls")
    parser.add_argument("opportunities")
    parser.add_argument("by_filter")
    parser.add_argument("funnel")
    parser.add_argument("--species", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--smrt-cells", required=True, metavar="TSV")
    parser.add_argument("--capped-quality", type=int, required=True, metavar="Q")
    parser.add_argument("--keep", action="append", default=[], metavar="COLUMN=VALUE[,VALUE...]")
    parser.add_argument("--calls-only", action="append", default=[], metavar="COLUMN=VALUE[,VALUE...]")
    parser.add_argument("--step", action="append", default=[], metavar="COLUMN=VALUE[,VALUE...]")
    args = parser.parse_args()

    con = connect()
    calls = Table(con, "calls", args.calls, None)
    opportunities = Table(con, "opportunities", args.opportunities, None, count_column="count")
    metrics = [column for column in calls.columns if column in opportunities.columns and column not in NOT_METRICS]
    calls_only = [column for column in (keep.split("=", 1)[0] for keep in args.calls_only)
                  if column in calls.columns and column not in opportunities.columns]
    write(by_filter_rows(con, calls, opportunities, metrics, calls_only), BY_FILTER_COLUMNS, args, args.by_filter)

    tables = [calls, opportunities]
    steps = final_steps(con, tables, capped_smrt_cells(args.smrt_cells, args.capped_quality), args.keep,
                        args.calls_only, lambda message: parser.error(f"--keep {message}"))
    for step in args.step:
        steps.append(keep_step(con, tables, steps, step, lambda message: parser.error(f"--step {message}")))
    write(funnel_rows(con, calls, opportunities, steps), FUNNEL_COLUMNS, args, args.funnel)


if __name__ == "__main__":
    main()
