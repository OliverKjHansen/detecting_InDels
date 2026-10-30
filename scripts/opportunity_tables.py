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
File"). In the pilot it once misread one instead, into a row without a count
(ticket 41), so a contig's cube is kept only if it holds the windows' total,
counted apart from DuckDB.
"""
import argparse
import gzip
import multiprocessing
import os
import shutil
import tempfile
from concurrent.futures import ProcessPoolExecutor

from duckdb_tables import columns, connect, cores, identifier, literal, source, write

# Keys the scans write as whole numbers, kept as numbers in the cube.
INTEGER_KEYS = ("size", "k")


def unpacked(windows, folder):
    """
    The window tables as plain text, the gzipped ones unpacked into `folder`,
    and the sum of their counts. One process per core: counting is Python's
    work, which threads would take turns at. The processes are started fresh
    (spawn), not forked from one that may hold DuckDB's threads.
    """
    plains = [os.path.join(folder, f"{number}.tsv") for number in range(len(windows))]
    with ProcessPoolExecutor(cores(), mp_context=multiprocessing.get_context("spawn")) as pool:
        done = list(pool.map(unpacked_window, map(str, windows), plains))
    return [plain for plain, _ in done], sum(total for _, total in done)


def unpacked_window(window, plain):
    """A window table as plain text, at `plain` if it was gzipped, and the sum of its counts."""
    if not window.endswith(".gz"):
        return window, counted(window, window)
    with gzip.open(window, "rb") as packed, open(plain, "wb") as out:
        shutil.copyfileobj(packed, out, 1 << 24)
    return plain, counted(plain, window)


def counted(table, name):
    """The sum of a window table's counts, its last column, each a whole number."""
    with open(table, "rb") as lines:
        if not next(lines, b"").rstrip(b"\r\n").endswith(b"\tcount"):
            raise ValueError(f"{name}: the last column isn't count")
        total = 0
        for line in lines:
            count = line[line.rfind(b"\t") + 1:].rstrip(b"\r\n")
            if not count.isdigit():
                raise ValueError(f"{name}: a count that isn't a whole number: {line!r}")
            total += int(count)
    return total


def holding(con, total):
    """A check that a contig's cube holds its windows' total, with no row without a count and no key with a NUL byte."""
    def check(cube):
        table = f"read_parquet({literal(cube)})"
        count = identifier("count")
        keys = ", ".join(f"CAST({identifier(column)} AS VARCHAR)" for column in columns(con, table) if column != "count")
        found = con.execute(f"SELECT coalesce(sum({count}), 0), count(*) FILTER (WHERE {count} IS NULL), "
                            f"count(*) FILTER (WHERE contains(concat_ws('', {keys}), chr(0))) FROM {table}").fetchone()
        if tuple(found) != (total, 0, 0):
            raise ValueError(f"{cube}: the cube holds {found[0]:,} opportunities, its windows {total:,}; "
                             f"{found[1]:,} rows have no count, {found[2]:,} a key with a NUL byte")
    return check


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
            windows, total = unpacked(args.windows, folder)
            write(con, summed(con, windows), args.output, check=holding(con, total))
    else:
        write(con, joined(args.contigs, args.species, args.sample), args.output)


if __name__ == "__main__":
    main()
