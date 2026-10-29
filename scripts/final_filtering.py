"""
The final filtering stage (Q7): the one definition of a final call, applied to
a sample's calls and its opportunities alike, so a filter never lowers the rate
by what it removes from the calls alone. Writes the tables the notebooks read.

    python final_filtering.py <calls.tsv> <opportunities.parquet> <final_calls.tsv> <final_opportunities.parquet> \\
        <filter_report.tsv> --platform <Revio | "Sequel II"> --smrt-cells <smrt_cells.tsv> --capped-quality <Q> \\
        [--keep COLUMN=VALUE[,VALUE...] ...]

A table is read and written as its name says: Parquet (the opportunity cube,
ticket 40), or TSV, gzipped if the name ends in .gz, its values as written.
DuckDB filters the tables with the job's cores, spilling to disk past its
memory (scripts/duckdb_tables.py).

The steps, in order, each on both tables:

  capped_smrt_cells  drop the SMRT Cells whose reads never have a base quality
                     above --capped-quality, by smrt_cell_qualities.py's table
                     of the sample (--smrt-cells): their base qualities top
                     out there (Q12, Q17). The same cells for every table of
                     the sample. The report's note lists them.
  COLUMN=VALUES      keep the rows whose COLUMN, a key of both tables, has one
                     of the values. Skipped, and noted, when every row that
                     reaches it is NA in the column: a tag the sample's
                     platform lacks (Q15).

Every row of both outputs gets a `platform` column. The report has one row per
step, `input` first: the calls and the opportunities (summed counts) before and
after it, the calls left without a haploid check, and a note.
"""
import argparse

import pandas

from duckdb_tables import columns, connect, identifier, literal, source, write


class Table:
    """
    One of the two tables a step filters: where it is read from and written
    to, and how it is counted, by rows, or by the sum of a count column.
    """

    def __init__(self, con, table_name, path, target, count_column=None):
        self.name, self.target, self.count_column = table_name, target, count_column
        self.header = columns(con, source(path))   # as written: "#contig" stays
        self.columns = [column.lstrip("#") for column in self.header]
        # The rows as an SQL table, its columns named without the "#".
        renamed = ", ".join(f"{identifier(as_written)} AS {identifier(column)}"
                            for as_written, column in zip(self.header, self.columns))
        self.relation = f"(SELECT {renamed} FROM {source(path)})"

    def size(self):
        """The SQL that counts the rows: their number, or the sum of the count column."""
        return "count(*)" if self.count_column is None else f"sum(CAST({identifier(self.count_column)} AS BIGINT))"

    def output_columns(self, platform):
        """The SQL of the columns written, named as read: the platform goes before the count, or last."""
        written = [f"{identifier(column)} AS {identifier(as_written)}"
                   for as_written, column in zip(self.header, self.columns)]
        platform_column = f"{literal(platform)} AS platform"
        if self.count_column is None:
            return written + [platform_column]
        at = self.header.index(self.count_column)
        return written[:at] + [platform_column] + written[at:]


class Step:
    """
    One step: the rows it keeps, as an SQL condition, and the rows of each
    table before and after it.
    """

    def __init__(self, step_name, keeps="true", note=""):
        self.name, self.keeps, self.note = step_name, keeps, note
        self.before, self.after = {}, {}
        self.haploid_check_missing = 0

    def row(self, has_haploid_check):
        return {"filter": self.name, "calls_before": self.before.get("calls", 0),
                "calls_after": self.after.get("calls", 0),
                "opportunities_before": self.before.get("opportunities", 0),
                "opportunities_after": self.after.get("opportunities", 0),
                "haploid_check_missing": self.haploid_check_missing if has_haploid_check else "NA",
                "note": self.note}


def kept(steps):
    """The SQL condition of the rows every one of the steps keeps."""
    return " AND ".join(f"({step.keeps})" for step in steps) or "true"


def has_values(con, tables, steps, column):
    """
    Whether a row that the steps leave, in either table, has a value other
    than NA in `column`: whether the sample has that tag (Q15).
    """
    return any(con.execute(f"SELECT bool_or(coalesce({identifier(column)}, '') <> 'NA') FROM {table.relation} "
                           f"WHERE {kept(steps)}").fetchone()[0] for table in tables)


def record_rows(con, table, steps):
    """Records each step's rows of `table` before and after it, in one pass over the table."""
    after = [f"{table.size()} FILTER (WHERE {kept(steps[:i + 1])})" for i in range(len(steps))]
    missing = []
    if "haploid_check" in table.columns:
        missing = [f"count(*) FILTER (WHERE {kept(steps[:i + 1])} AND {identifier('haploid_check')} = 'NA')"
                   for i in range(len(steps))]
    total, *counts = con.execute(f"SELECT {', '.join([table.size()] + after + missing)} FROM {table.relation}").fetchone()
    before = total or 0
    for i, step in enumerate(steps):
        step.before[table.name], step.after[table.name] = before, counts[i] or 0
        before = step.after[table.name]
        if missing:
            step.haploid_check_missing += counts[len(steps) + i]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("calls")
    parser.add_argument("opportunities")
    parser.add_argument("final_calls")
    parser.add_argument("final_opportunities")
    parser.add_argument("report")
    parser.add_argument("--platform", required=True)
    parser.add_argument("--smrt-cells", required=True, metavar="TSV")
    parser.add_argument("--capped-quality", type=int, required=True, metavar="Q")
    parser.add_argument("--keep", action="append", default=[], metavar="COLUMN=VALUE[,VALUE...]")
    args = parser.parse_args()

    con = connect()
    tables = [Table(con, "calls", args.calls, args.final_calls),
              Table(con, "opportunities", args.opportunities, args.final_opportunities, count_column="count")]

    cells = pandas.read_csv(args.smrt_cells, sep="\t", dtype={"smrt_cell": str})
    capped = sorted(cells.smrt_cell[cells.max_quality <= args.capped_quality])
    steps = [Step("input"),
             Step("capped_smrt_cells",
                  f"coalesce({identifier('smrt_cell')}, '') NOT IN ({', '.join(map(literal, capped))})" if capped else "true",
                  note=",".join(capped))]
    for keep in args.keep:
        column, values = keep.split("=", 1)
        for table in tables:
            if column not in table.columns:
                parser.error(f"--keep {keep}: the {table.name} have no column {column}")
        # A tag the sample lacks is NA in every row that reaches the step: skip it, rather than remove everything.
        if has_values(con, tables, steps, column):
            steps.append(Step(keep, f"{identifier(column)} IN ({', '.join(map(literal, values.split(',')))})"))
        else:
            steps.append(Step(keep, note="skipped: NA throughout the sample"))

    for table in tables:
        record_rows(con, table, steps)
        write(con, f"SELECT {', '.join(table.output_columns(args.platform))} FROM {table.relation} WHERE {kept(steps)}",
              table.target)

    has_haploid_check = "haploid_check" in tables[0].columns
    pandas.DataFrame([step.row(has_haploid_check) for step in steps]).to_csv(args.report, sep="\t", index=False)


if __name__ == "__main__":
    main()
