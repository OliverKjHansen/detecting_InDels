"""
The opportunity cube (ticket 40): a sample's opportunities, the read positions
that could carry an InDel, counted for every combination of the keys the scan
wrote (contig, InDel type and size or motif and period, context, the strata),
and kept in Parquet, for DuckDB to filter and sum by any of them.

    python opportunity_tables.py sum <contig.parquet> <window.tsv.gz> ...
    python opportunity_tables.py join <merged.parquet> --species S --sample S <contig.parquet> ...

sum   a contig's window tables: the count summed per key, sorted by key. The
      keys are text, as the scan wrote them, but `size` and `k`, numbers.
join  a sample's contigs, one after another, with its species and sample.

The output is written whole or not at all. The windows are unpacked into the
job's $TMPDIR first: DuckDB 2.0 can fail to read a gzipped table when short of
memory, as it rereads parts of the file ("Trying to seek a compressed CSV
File").
"""
import argparse
import gzip
import os
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor

from duckdb_tables import columns, connect, cores, identifier, literal, source, write

# Keys the scans write as whole numbers, kept as numbers in the cube.
INTEGER_KEYS = ("size", "k")


def unpacked(windows, folder):
    """The window tables as plain text: the gzipped ones unpacked into `folder`, the job's cores at a time."""
    def unpack(numbered):
        number, window = numbered
        if not str(window).endswith(".gz"):
            return str(window)
        plain = os.path.join(folder, f"{number}.tsv")
        with gzip.open(window, "rb") as packed, open(plain, "wb") as out:
            shutil.copyfileobj(packed, out, 1 << 24)
        return plain
    with ThreadPoolExecutor(cores()) as pool:
        return list(pool.map(unpack, enumerate(windows)))


def summed(con, windows):
    """The query that sums the windows' counts per key."""
    table = source(*windows)
    keys = [f"CAST({identifier(column)} AS INTEGER) AS {identifier(column)}" if column in INTEGER_KEYS
            else identifier(column) for column in columns(con, table) if column != "count"]
    count = identifier("count")
    # A sum of integers is a 128-bit one, which Parquet would store as a double.
    return (f"SELECT {', '.join(keys)}, CAST(sum(CAST({count} AS BIGINT)) AS BIGINT) AS {count} FROM {table} "
            "GROUP BY ALL ORDER BY ALL")


def joined(contigs, species, sample):
    """The query that puts the contigs' tables one after another, with the species and sample."""
    return f"SELECT *, {literal(species)} AS species, {literal(sample)} AS sample FROM {source(*contigs)}"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    sum_parser = commands.add_parser("sum")
    sum_parser.add_argument("output")
    sum_parser.add_argument("windows", nargs="+")
    join_parser = commands.add_parser("join")
    join_parser.add_argument("output")
    join_parser.add_argument("--species", required=True)
    join_parser.add_argument("--sample", required=True)
    join_parser.add_argument("contigs", nargs="+")
    args = parser.parse_args()

    con = connect()
    if args.command == "sum":
        with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as folder:
            write(con, summed(con, unpacked(args.windows, folder)), args.output)
    else:
        write(con, joined(args.contigs, args.species, args.sample), args.output)


if __name__ == "__main__":
    main()
