"""
Compare a table from a local replay of the pipeline with the cluster's.

    python scripts/compare_tables.py <local.tsv> <cluster.tsv> [--label NAME]

Prints "same" or what differs, and exits 1 when they differ. The tables are
compared by content: as multisets of rows, so order doesn't matter, and with
numbers compared by value, so 20, 20.0 and 2e1 are the same: whole numbers
exactly, others to 9 significant digits. The species and sample columns are
ignored when only the cluster's table has them, because the workflows append
them after the scripts have run. Any other column on one side only is a
difference. In a table with a count column, rows counting 0 carry nothing and
are ignored: the scans once wrote one for every kmerQ level.
"""

import argparse
import csv
import re
import sys
from collections import Counter, namedtuple

Comparison = namedtuple("Comparison", "same only_local only_cluster local_columns cluster_columns")

# Appended by the workflows to the tables after the haploid check and the
# final filters.
APPENDED = ("species", "sample")

NUMBER = re.compile(r"[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?")


def normalise(cell):
    """A number in one fixed format; anything else as it is."""
    if not NUMBER.fullmatch(cell):
        return cell
    value = float(cell)
    if value.is_integer() and abs(value) < 2 ** 53:   # exact as a float
        return str(int(value))
    return f"{value:.9g}"


def read_table(path):
    """(column names, rows) of a TSV with one header line, which may start with '#'."""
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = [name.lstrip("#") for name in next(reader)]
        rows = [row for row in reader if row]
    if "count" in header:
        count = header.index("count")
        rows = [row for row in rows if normalise(row[count]) != "0"]
    return header, rows


def rows_on(columns, header, rows):
    """The rows cut down to `columns`, in that order, as a multiset."""
    indices = [header.index(column) for column in columns]
    return Counter(tuple(normalise(row[i]) for i in indices) for row in rows)


def compare(local_path, cluster_path):
    """
    Comparison(same, only_local, only_cluster, local_columns, cluster_columns):
    the rows on one side only, as multisets over the shared columns, and the
    columns on one side only, leaving out the ones the workflows append.
    """
    local_header, local_rows = read_table(local_path)
    cluster_header, cluster_rows = read_table(cluster_path)
    shared = [column for column in local_header if column in cluster_header]
    local_rows = rows_on(shared, local_header, local_rows)
    cluster_rows = rows_on(shared, cluster_header, cluster_rows)
    local_columns = [column for column in local_header if column not in cluster_header]
    cluster_columns = [column for column in cluster_header
                       if column not in local_header and column not in APPENDED]
    return Comparison(same=local_rows == cluster_rows and not local_columns and not cluster_columns,
                      only_local=local_rows - cluster_rows,
                      only_cluster=cluster_rows - local_rows,
                      local_columns=local_columns,
                      cluster_columns=cluster_columns)


def report(label, comparison, shown=5):
    """What the replay prints for one table: "same", or what differs."""
    if comparison.same:
        return f"same:      {label}"
    lines = [f"DIFFERENT: {label}"]
    for side, columns in (("local", comparison.local_columns), ("cluster", comparison.cluster_columns)):
        if columns:
            lines.append(f"  columns only the {side} table has: {', '.join(columns)}")
    for side, rows in (("local", comparison.only_local), ("cluster", comparison.only_cluster)):
        if rows:
            lines.append(f"  {sum(rows.values()):,} rows only in the {side} table, for example:")
            lines += ["    " + "\t".join(row) for row in list(rows.elements())[:shown]]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split("\n")[0])
    parser.add_argument("local")
    parser.add_argument("cluster")
    parser.add_argument("--label", help="what to call the table in the report (default: the local path)")
    args = parser.parse_args()
    comparison = compare(args.local, args.cluster)
    print(report(args.label or args.local, comparison))
    sys.exit(0 if comparison.same else 1)


if __name__ == "__main__":
    main()
