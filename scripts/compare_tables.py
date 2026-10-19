"""
Compare a table from a local replay of the pipeline with the cluster's.

    python scripts/compare_tables.py <local.tsv> <cluster.tsv> [--label NAME]

Prints "same" or what differs, and exits 1 when they differ. The tables are
compared by content: as multisets of rows, so order doesn't matter, and with
numbers compared by value to 9 significant digits, so 20, 20.0 and 2e1 are the
same. Columns only the cluster's table has are ignored: the workflows append
species and sample after the scripts have run. A column only the local table
has is a difference.
"""

import argparse
import csv
import re
import sys
from collections import Counter, namedtuple

Comparison = namedtuple("Comparison", "same only_local only_cluster local_columns ignored_columns")

NUMBER = re.compile(r"[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?")


def normalise(cell):
    """A number in one fixed format; anything else as it is."""
    return f"{float(cell):.9g}" if NUMBER.fullmatch(cell) else cell


def read_table(path):
    """(column names, rows) of a TSV with one header line, which may start with '#'."""
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = [name.lstrip("#") for name in next(reader)]
        return header, [row for row in reader if row]


def rows_on(columns, header, rows):
    """The rows cut down to `columns`, in that order, as a multiset."""
    indices = [header.index(column) for column in columns]
    return Counter(tuple(normalise(row[i]) for i in indices) for row in rows)


def compare(local_path, cluster_path):
    """
    Comparison(same, only_local, only_cluster, local_columns, ignored_columns):
    the rows on one side only, as multisets over the shared columns; the
    columns only the local table has; those only the cluster's has.
    """
    local_header, local_rows = read_table(local_path)
    cluster_header, cluster_rows = read_table(cluster_path)
    shared = [column for column in local_header if column in cluster_header]
    local_rows = rows_on(shared, local_header, local_rows)
    cluster_rows = rows_on(shared, cluster_header, cluster_rows)
    local_columns = [column for column in local_header if column not in cluster_header]
    return Comparison(same=local_rows == cluster_rows and not local_columns,
                      only_local=local_rows - cluster_rows,
                      only_cluster=cluster_rows - local_rows,
                      local_columns=local_columns,
                      ignored_columns=[column for column in cluster_header if column not in local_header])


def report(label, comparison, shown=5):
    """What the replay prints for one table: "same", or what differs."""
    if comparison.same:
        return f"same:      {label}"
    lines = [f"DIFFERENT: {label}"]
    if comparison.local_columns:
        lines.append(f"  columns only the local table has: {', '.join(comparison.local_columns)}")
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
