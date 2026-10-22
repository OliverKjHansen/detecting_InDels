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
                     of the values. Skipped, and noted, for a sample whose
                     tables are NA in the column throughout: a tag its
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


def read_chunks(path, rows):
    """A table, `rows` rows at a time, as text: NA and the numbers are written back as they were read."""
    return pandas.read_csv(path, sep="\t", dtype=str, keep_default_na=False, chunksize=rows)


def header(path):
    with open_text(path) as handle:
        return handle.readline().rstrip("\n").split("\t")


class Step:
    """One step of the report: what it keeps, and the rows before and after it, summed over the chunks."""

    def __init__(self, name, keeps, note=""):
        self.name, self.keeps, self.note = name, keeps, note
        self.calls = [0, 0]
        self.opportunities = [0, 0]
        self.haploid_check_missing = 0

    def apply(self, table, kind):
        """The rows of a chunk the step keeps, counted before and after."""
        kept = table[self.keeps(table)]
        totals = self.calls if kind == "calls" else self.opportunities
        size = (lambda rows: len(rows)) if kind == "calls" else (lambda rows: int(rows["count"].astype(int).sum()))
        totals[0] += size(table)
        totals[1] += size(kept)
        if kind == "calls" and "haploid_check" in kept:
            self.haploid_check_missing += int((kept["haploid_check"] == "NA").sum())
        return kept

    def row(self, has_haploid_check):
        return {"filter": self.name, "calls_before": self.calls[0], "calls_after": self.calls[1],
                "opportunities_before": self.opportunities[0], "opportunities_after": self.opportunities[1],
                "haploid_check_missing": self.haploid_check_missing if has_haploid_check else "NA",
                "note": self.note}


def everything(table):
    return pandas.Series(True, index=table.index)


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

    columns = {"calls": [name.lstrip("#") for name in header(args.calls)], "opportunities": header(args.opportunities)}
    keeps = []
    for keep in args.keep:
        column, values = keep.split("=", 1)
        for kind in columns:
            if column not in columns[kind]:
                parser.error(f"--keep {keep}: the {kind} have no column {column}")
        keeps.append((keep, column, values.split(",")))

    # A first pass for the keys the sample has: a tag its platform lacks is NA throughout.
    present = {column: False for _, column, _ in keeps}
    for path in (args.calls, args.opportunities):
        for table in read_chunks(path, args.chunk_rows):
            for column in present:
                if not present[column]:
                    present[column] = bool((table[column] != "NA").any())

    cells = pandas.read_csv(args.smrt_cells, sep="\t", dtype={"smrt_cell": str})
    capped = sorted(cells.smrt_cell[cells.max_quality <= args.capped_quality])
    steps = [Step("input", everything),
             Step("capped_smrt_cells", lambda table: ~table["smrt_cell"].isin(capped), ",".join(capped))]
    for keep, column, values in keeps:
        if present[column]:
            steps.append(Step(keep, lambda table, column=column, values=values: table[column].isin(values)))
        else:
            steps.append(Step(keep, everything, "skipped: NA throughout the sample"))

    # Every row gets the platform: after the calls' columns, before the opportunities' count.
    for kind, source, target in (("calls", args.calls, args.final_calls),
                                 ("opportunities", args.opportunities, args.final_opportunities)):
        names = header(source)
        if kind == "calls":
            out_names = names + ["platform"]
        else:
            out_names = names[:names.index("count")] + ["platform"] + names[names.index("count"):]
        with open_text(target, "w") as out:
            out.write("\t".join(out_names) + "\n")
            for table in read_chunks(source, args.chunk_rows):
                table.columns = [name.lstrip("#") for name in table.columns]
                for step in steps:
                    table = step.apply(table, kind)
                table = table.assign(platform=args.platform)
                table[[name.lstrip("#") for name in out_names]].to_csv(out, sep="\t", index=False, header=False)

    has_haploid_check = "haploid_check" in columns["calls"]
    pandas.DataFrame([step.row(has_haploid_check) for step in steps]).to_csv(args.report, sep="\t", index=False)


if __name__ == "__main__":
    main()
