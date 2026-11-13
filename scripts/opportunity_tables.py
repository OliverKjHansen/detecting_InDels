"""
The opportunity cube (ticket 40): a sample's opportunities, the read positions
that could carry an InDel, counted for every combination of the keys the scan
wrote (contig, InDel type and size or motif and period, context, the strata),
and kept in Parquet, for DuckDB to filter and sum by any of them.

    python opportunity_tables.py sum <contig.parquet> <window.tsv.gz> ...
    python opportunity_tables.py join <merged.parquet> --species S --sample S <contig.parquet> ...

sum   a contig's window tables: the count summed per key, sorted by key as
      text, byte by byte. The keys are text, as the scan wrote them, but `size`
      and `k`, numbers.
join  a sample's contigs, one after another, with its species and sample.

The output is written whole or not at all. The sum never holds the keys in
memory (ticket 50): a contig can have hundreds of millions (374 M in one of
ob007's), and DuckDB's GROUP BY ran out of memory on them, then failed in the
code that spills them to disk. Instead `sort` puts the windows' lines in order,
in the memory it is given and spilling to the job's $TMPDIR, and the lines of a
key, now next to each other, are added up one key at a time, with the windows'
total; `sort` must give back as many lines as it was given. DuckDB only writes
the summed lines to Parquet. It has misread tables at random (ticket 41): into
a row without a count or a key with a NUL byte, or by losing lines, which only
the total shows. So a contig's cube is kept only if it holds the windows'
total, with no such row.
"""
import argparse
import gzip
import os
import subprocess
import tempfile

from duckdb_tables import columns, connect, cores, half_memory_mb, identifier, literal, source, write

# Keys the scans write as whole numbers, kept as numbers in the cube.
INTEGER_KEYS = ("size", "k")


def opened(window):
    """A window table, gzipped or plain, opened for reading bytes."""
    return gzip.open(window, "rb") if str(window).endswith(".gz") else open(window, "rb")


def header(window):
    """A window table's first line, as written."""
    with opened(window) as lines:
        return lines.readline()


def feed(window, into):
    """
    Writes a window table's lines, without its header, into `into`, and
    returns how many. A NUL byte, or a last line without its newline, fails it.
    """
    count, last = 0, b"\n"
    with opened(window) as lines:
        lines.readline()
        while chunk := lines.read(1 << 24):
            if b"\0" in chunk:
                raise ValueError(f"{window}: a NUL byte")
            into.write(chunk)
            count += chunk.count(b"\n")
            last = chunk[-1:]
    if last != b"\n":
        raise ValueError(f"{window}: the last line has no newline")
    return count


def sort_options(folder):
    """sort's options: the job's cores, half its memory, spilling to `folder`."""
    options = ["-T", folder, f"--parallel={cores()}"]
    if half_memory_mb():
        options += ["-S", f"{half_memory_mb()}M"]
    return options


def write_sums(lines, out, fields):
    """
    Writes each key of `lines`, sorted by key, once to `out`, with the sum of
    its counts; returns how many lines there were and their total. A key is
    everything before the last tab. A key smaller than the one before means the
    lines weren't sorted as compared here, and lines of one key may not have met.
    """
    previous, read, total, count_of_key = None, 0, 0, 0
    for line in lines:
        read += 1
        key, _, text = line.rstrip(b"\n").rpartition(b"\t")
        tabs = line.count(b"\t")
        if tabs != fields - 1:
            raise ValueError(f"a line with {tabs + 1} fields, not {fields}: {line!r}")
        if not text.isdigit():
            raise ValueError(f"a count that isn't a whole number: {line!r}")
        count = int(text)
        total += count
        if key == previous:
            count_of_key += count
            continue
        if previous is not None:
            if key < previous:
                raise ValueError(f"the lines aren't sorted: {key!r} after {previous!r}")
            out.write(b"%s\t%d\n" % (previous, count_of_key))
        previous, count_of_key = key, count
    if previous is not None:
        out.write(b"%s\t%d\n" % (previous, count_of_key))
    return read, total


def added_up(windows, folder):
    """
    The windows' lines summed per key into a TSV in `folder`, sorted by key as
    text, and the windows' total. The windows must have the same columns, count
    the last. sort reads every line before it writes one, so the windows are
    fed in full before its output is read.
    """
    headers = {header(window) for window in windows}
    if len(headers) != 1:
        raise ValueError("the windows have different columns: " + " | ".join(sorted(map(repr, headers))))
    first = headers.pop()
    names = first.rstrip(b"\n").split(b"\t")   # a CRLF window fails here, at count\r
    if names[-1] != b"count":
        raise ValueError(f"the last column isn't count: {first!r}")
    summed = os.path.join(folder, "summed.tsv")
    # LC_ALL=C: byte order, as Python compares the keys.
    with subprocess.Popen(["sort", *sort_options(folder)], env=dict(os.environ, LC_ALL="C"),
                          stdin=subprocess.PIPE, stdout=subprocess.PIPE) as sort:
        try:
            fed = sum(feed(window, sort.stdin) for window in windows)
            sort.stdin.close()
            with open(summed, "wb") as out:
                out.write(first)
                read, total = write_sums(sort.stdout, out, len(names))
        except BrokenPipeError:
            raise RuntimeError(f"sort stopped reading the windows (exit code {sort.wait()})") from None
        except BaseException:
            sort.kill()
            raise
    if sort.returncode:
        raise RuntimeError(f"sort failed (exit code {sort.returncode})")
    if read != fed:
        raise RuntimeError(f"sort was given {fed:,} lines and gave back {read:,}")
    return summed, total


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


def typed(con, summed):
    """The query that reads the summed lines, in their order, `size`, `k` and the count as numbers."""
    table = source(summed)
    keys = [f"CAST({identifier(column)} AS INTEGER) AS {identifier(column)}" if column in INTEGER_KEYS
            else identifier(column) for column in columns(con, table) if column != "count"]
    count = identifier("count")
    # A count as BIGINT: Parquet would store a 128-bit integer as a double.
    return f"SELECT {', '.join(keys)}, CAST({count} AS BIGINT) AS {count} FROM {table}"


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

    if args.command == "sum":
        with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as folder:
            summed, total = added_up(args.windows, folder)
            con = connect()   # after sort is done: each is given half the job's memory
            write(con, typed(con, summed), args.output, check=holding(con, total))
    else:
        con = connect()
        write(con, joined(args.contigs, args.species, args.sample), args.output)


if __name__ == "__main__":
    main()
