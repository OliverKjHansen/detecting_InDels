"""
The final filtering stage (Q7): the one definition of a final call, applied to
a sample's calls and its opportunities alike, so a filter never lowers the rate
by what it removes from the calls alone. Writes the tables the notebooks read.

    python final_filtering.py <calls.tsv> <opportunities.tsv> <final_calls.tsv> <final_opportunities.tsv> \\
        <filter_report.tsv> --platform <Revio | "Sequel II"> --smrt-cells <smrt_cells.tsv> --capped-quality <Q> \\
        [--keep COLUMN=VALUE[,VALUE...] ...]

The steps, in order, each on both tables:

  capped_smrt_cells  drop the SMRT Cells whose reads never have a base quality
                     above --capped-quality, by smrt_cell_qualities.py's table
                     of the sample (--smrt-cells): their base qualities top
                     out there (Q12, Q17). The same cells for every table of
                     the sample. The report's note lists them.
  COLUMN=VALUES      keep the rows whose COLUMN, a key of both tables, has one
                     of the values. Skipped, and noted, for a sample whose rows
                     are NA in the column throughout: a tag its platform lacks
                     (Q15).

Every row of both outputs gets a `platform` column. The report has one row per
step, `input` first: the calls and the opportunities (summed counts) before and
after it, the calls left without a haploid check, and a note.
"""
import argparse

import pandas


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

    # As text, so NA and the numbers are written back as they were read.
    calls = pandas.read_csv(args.calls, sep="\t", dtype=str, keep_default_na=False)
    opportunities = pandas.read_csv(args.opportunities, sep="\t", dtype=str, keep_default_na=False)
    opportunities["count"] = opportunities["count"].astype(int)
    report = []

    def step(name, kept_calls, kept_opportunities, note=""):
        nonlocal calls, opportunities
        row = {"filter": name, "calls_before": len(calls), "calls_after": int(kept_calls.sum()),
               "opportunities_before": int(opportunities["count"].sum())}
        calls, opportunities = calls[kept_calls], opportunities[kept_opportunities]
        row["opportunities_after"] = int(opportunities["count"].sum())
        row["haploid_check_missing"] = (int((calls["haploid_check"] == "NA").sum())
                                        if "haploid_check" in calls else "NA")
        row["note"] = note
        report.append(row)

    def everything(table):
        return pandas.Series(True, index=table.index)

    step("input", everything(calls), everything(opportunities))

    cells = pandas.read_csv(args.smrt_cells, sep="\t", dtype={"smrt_cell": str})
    capped = sorted(cells.smrt_cell[cells.max_quality <= args.capped_quality])
    step("capped_smrt_cells", ~calls["smrt_cell"].isin(capped), ~opportunities["smrt_cell"].isin(capped),
         ",".join(capped))

    for keep in args.keep:
        column, values = keep.split("=", 1)
        values = values.split(",")
        for name, table in (("calls", calls), ("opportunities", opportunities)):
            if column not in table:
                parser.error(f"--keep {keep}: the {name} have no column {column}")
        if not ((calls[column] != "NA").any() or (opportunities[column] != "NA").any()):
            step(keep, everything(calls), everything(opportunities), "skipped: NA throughout the sample")
        else:
            step(keep, calls[column].isin(values), opportunities[column].isin(values))

    calls = calls.assign(platform=args.platform)
    opportunities.insert(opportunities.columns.get_loc("count"), "platform", args.platform)
    calls.to_csv(args.final_calls, sep="\t", index=False)
    opportunities.to_csv(args.final_opportunities, sep="\t", index=False)
    pandas.DataFrame(report).to_csv(args.report, sep="\t", index=False)


if __name__ == "__main__":
    main()
