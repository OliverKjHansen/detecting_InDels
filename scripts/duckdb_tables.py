"""
The big tables through DuckDB (ticket 40): read and written by the extension of
their names, .parquet, or TSV, plain or gzipped, and queried with the job's own
cores and memory, spilling to its temporary folder.
"""
import os

import duckdb


def connect():
    """
    A DuckDB connection with the job's cores and three quarters of its memory,
    spilling to its $TMPDIR. gwf asks Slurm for the cores with -c and the
    memory with --mem; outside a job, the machine's.
    """
    con = duckdb.connect()
    con.execute("SET enable_progress_bar = false")   # it would fill the job's log
    con.execute(f"SET threads = {os.environ.get('SLURM_CPUS_PER_TASK') or os.cpu_count()}")
    if os.environ.get("SLURM_MEM_PER_NODE"):   # MB
        con.execute(f"SET memory_limit = '{int(os.environ['SLURM_MEM_PER_NODE']) * 3 // 4}MB'")
    if os.environ.get("TMPDIR"):
        con.execute(f"SET temp_directory = {literal(os.environ['TMPDIR'])}")
    return con


def literal(text):
    """Text as an SQL string."""
    return "'" + str(text).replace("'", "''") + "'"


def identifier(column):
    """A column's name as an SQL identifier."""
    return '"' + column.replace('"', '""') + '"'


def source(*paths):
    """
    The tables at `paths`, one after another, as an SQL table: Parquet as
    stored, or TSV (plain or gzipped) as text exactly as written, an empty
    field as an empty string.
    """
    listed = "[" + ", ".join(literal(path) for path in paths) + "]"
    parquet = [str(path).endswith(".parquet") for path in paths]
    if all(parquet):
        return f"read_parquet({listed})"
    if any(parquet):
        raise ValueError(f"Parquet and TSV tables together: {', '.join(map(str, paths))}")
    # No field is \N, so none is read as NULL.
    return f"read_csv({listed}, delim='\t', header=true, all_varchar=true, quote='', escape='', nullstr='\\N')"


def columns(con, table):
    """The names of an SQL table's columns."""
    return [row[0] for row in con.execute(f"DESCRIBE SELECT * FROM {table}").fetchall()]


def write(con, query, path):
    """
    Writes the rows of `query` to `path`, Parquet or TSV by its name, to a
    temporary file first, renamed once complete.
    """
    path = str(path)
    if path.endswith(".parquet"):
        options = "FORMAT parquet, COMPRESSION zstd"
    else:
        options = "FORMAT csv, DELIMITER '\t', HEADER, QUOTE ''" + (", COMPRESSION gzip" if path.endswith(".gz") else "")
    con.execute(f"COPY ({query}) TO {literal(path + '.tmp')} ({options})")
    os.replace(path + ".tmp", path)
