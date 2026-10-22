"""
The final filtering stage (Q7): the one definition of a final call, applied to
a sample's calls and its opportunities alike, so a filter never lowers the rate
by what it removes from the calls alone. Writes the tables the notebooks read.

    python final_filtering.py <calls.tsv> <opportunities.tsv> <final_calls.tsv> <final_opportunities.tsv> \\
        <filter_report.tsv> --platform <Revio | "Sequel II"> --smrt-cells <smrt_cells.tsv> --capped-quality <Q> \\
        [--keep COLUMN=VALUE[,VALUE...] ...] [--chunk-rows N]

A table whose name ends in .gz is read or written gzipped. The tables are read
--chunk-rows rows at a time, so memory doesn't grow with them; the result is
the same whatever the chunk size.

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

from text_files import open_text

# Rows read at a time: a few GB of memory at most.
CHUNK_ROWS = 1_000_000


class Table:
    """
    One of the two tables a step filters: where it is read from and written
    to, and how it is counted, by rows, or by the sum of a count column.
    """

    def __init__(self, name, source, target, count_column=None):
        self.name, self.source, self.target, self.count_column = name, source, target, count_column
        with open_text(source) as handle:
            self.header = handle.readline().rstrip("\n").split("\t")   # as written: "#contig" stays
        self.columns = [name.lstrip("#") for name in self.header]

    def chunks(self, rows, columns=None):
        """The table `rows` rows at a time, as text, so NA and the numbers are written back as read."""
        for chunk in pandas.read_csv(self.source, sep="\t", dtype=str, keep_default_na=False, chunksize=rows,
                                     usecols=None if columns is None else sorted(columns)):
            chunk.columns = [name.lstrip("#") for name in chunk.columns]
            yield chunk

    def size(self, rows):
        return len(rows) if self.count_column is None else int(rows[self.count_column].astype(int).sum())

    def output_header(self):
        """The header written: the platform goes before the count, or last."""
        if self.count_column is None:
            return self.header + ["platform"]
        at = self.header.index(self.count_column)
        return self.header[:at] + ["platform"] + self.header[at:]


class Step:
    """
    One step: the rows it keeps (`column`, the key it filters on, or None),
    and the rows of each table before and after it, summed over the chunks.
    """

    def __init__(self, name, keeps, column=None, note=""):
        self.name, self.keeps, self.column, self.note = name, keeps, column, note
        self.before, self.after = {}, {}
        self.haploid_check_missing = 0

    def apply(self, chunk, table):
        """The rows of a chunk of `table` the step keeps, counted before and after."""
        kept = chunk[self.keeps(chunk)]
        self.before[table.name] = self.before.get(table.name, 0) + table.size(chunk)
        self.after[table.name] = self.after.get(table.name, 0) + table.size(kept)
        if "haploid_check" in kept:
            self.haploid_check_missing += int((kept["haploid_check"] == "NA").sum())
        return kept

    def row(self, has_haploid_check):
        return {"filter": self.name, "calls_before": self.before.get("calls", 0),
                "calls_after": self.after.get("calls", 0),
                "opportunities_before": self.before.get("opportunities", 0),
                "opportunities_after": self.after.get("opportunities", 0),
                "haploid_check_missing": self.haploid_check_missing if has_haploid_check else "NA",
                "note": self.note}


def everything(table):
    return pandas.Series(True, index=table.index)


def has_values(tables, steps, column, rows):
    """
    Whether a row that the steps leave, in either table, has a value other
    than NA in `column`: whether the sample has that tag (Q15). Reads only
    the columns the steps and this one need.
    """
    needed = {"smrt_cell", column} | {step.column for step in steps if step.column}
    for table in tables:
        for chunk in table.chunks(rows, needed):
            for step in steps:
                chunk = chunk[step.keeps(chunk)]
            if (chunk[column] != "NA").any():
                return True
    return False


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
    parser.add_argument("--chunk-rows", type=int, default=CHUNK_ROWS, metavar="N")
    args = parser.parse_args()

    tables = [Table("calls", args.calls, args.final_calls),
              Table("opportunities", args.opportunities, args.final_opportunities, count_column="count")]

    cells = pandas.read_csv(args.smrt_cells, sep="\t", dtype={"smrt_cell": str})
    capped = sorted(cells.smrt_cell[cells.max_quality <= args.capped_quality])
    steps = [Step("input", everything),
             Step("capped_smrt_cells", lambda chunk: ~chunk["smrt_cell"].isin(capped), note=",".join(capped))]
    for keep in args.keep:
        column, values = keep.split("=", 1)
        for table in tables:
            if column not in table.columns:
                parser.error(f"--keep {keep}: the {table.name} have no column {column}")
        # A tag the sample lacks is NA in every row that reaches the step: skip it, rather than remove everything.
        if has_values(tables, steps, column, args.chunk_rows):
            steps.append(Step(keep, lambda chunk, column=column, values=values.split(","): chunk[column].isin(values),
                              column=column))
        else:
            steps.append(Step(keep, everything, note="skipped: NA throughout the sample"))

    for table in tables:
        out_header = table.output_header()
        out_columns = [name.lstrip("#") for name in out_header]
        with open_text(table.target, "w") as out:
            out.write("\t".join(out_header) + "\n")
            for chunk in table.chunks(args.chunk_rows):
                for step in steps:
                    chunk = step.apply(chunk, table)
                chunk.assign(platform=args.platform)[out_columns].to_csv(out, sep="\t", index=False, header=False)

    has_haploid_check = "haploid_check" in tables[0].columns
    pandas.DataFrame([step.row(has_haploid_check) for step in steps]).to_csv(args.report, sep="\t", index=False)


if __name__ == "__main__":
    main()
