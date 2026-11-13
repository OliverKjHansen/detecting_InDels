"""
The opportunity cube (ticket 40): scripts/opportunity_tables.py sums a contig's
window tables into one Parquet file, and joins a sample's contigs into one,
for DuckDB to filter and sum by any set of variables.
"""

import gzip
import os
import subprocess
import sys

import duckdb
import pytest

from duckdb_tables import write
from opportunity_tables import holding

SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts", "opportunity_tables.py")
HEADER = "contig\tindel_type\tsize\tcontext\ttissue\tcount\n"


def window(path, rows):
    """A scan's window table, gzipped, as the scans write it."""
    path.write_bytes(gzip.compress((HEADER + "".join("\t".join(map(str, row)) + "\n" for row in rows)).encode()))
    return path


def opportunity_tables(*arguments):
    return subprocess.run([sys.executable, SCRIPT, *map(str, arguments)], capture_output=True, text=True)


def rows(path):
    return duckdb.sql(f"SELECT * FROM read_parquet('{path}')").fetchall()


def test_a_contigs_windows_are_summed_per_key_into_parquet(tmp_path):
    windows = [window(tmp_path / "a.tsv.gz", [("c1", "insertion", 1, "A", "liver", 5), ("c1", "deletion", 2, "AC", "NA", 1)]),
               window(tmp_path / "b.tsv.gz", [("c1", "insertion", 1, "A", "liver", 2), ("c1", "deletion", 3, "ACG", "NA", 7)]),
               window(tmp_path / "c.tsv.gz", [])]   # a window without opportunities
    out = tmp_path / "c1.parquet"
    done = opportunity_tables("sum", out, *windows)
    assert done.returncode == 0, done.stderr
    assert sorted(rows(out)) == [("c1", "deletion", 2, "AC", "NA", 1), ("c1", "deletion", 3, "ACG", "NA", 7),
                                 ("c1", "insertion", 1, "A", "liver", 7)]


def test_counts_past_four_billion_add_up_exactly(tmp_path):
    # A contig's sums reach billions (ticket 41: 6,307,196,945 in one of carl's).
    windows = [window(tmp_path / "a.tsv.gz", [("c1", "insertion", 1, "A", "liver", 4_000_000_000)]),
               window(tmp_path / "b.tsv.gz", [("c1", "insertion", 1, "A", "liver", 4_000_000_001)])]
    out = tmp_path / "c1.parquet"
    done = opportunity_tables("sum", out, *windows)
    assert done.returncode == 0, done.stderr
    assert rows(out) == [("c1", "insertion", 1, "A", "liver", 8_000_000_001)]


def test_the_cubes_rows_are_sorted_by_key_as_text(tmp_path):
    # Methylation bins where one starts like another ("1" and "12", "-1" and "-10").
    header = b"contig\tmethylation_bin\tcount\n"
    a, b = tmp_path / "a.tsv.gz", tmp_path / "b.tsv.gz"
    a.write_bytes(gzip.compress(header + b"c1\t12\t1\nc1\t1\t20\n"))
    b.write_bytes(gzip.compress(header + b"c1\t1\t5\nc1\t-1\t3\nc1\t-10\t4\n"))
    out = tmp_path / "c1.parquet"
    done = opportunity_tables("sum", out, a, b)
    assert done.returncode == 0, done.stderr
    assert rows(out) == [("c1", "-1", 3), ("c1", "-10", 4), ("c1", "1", 25), ("c1", "12", 1)]


def test_the_cube_keeps_the_scans_columns_with_size_and_count_as_numbers(tmp_path):
    out = tmp_path / "c1.parquet"
    done = opportunity_tables("sum", out, window(tmp_path / "a.tsv.gz", [("c1", "insertion", 1, "A", "liver", 5)]))
    assert done.returncode == 0, done.stderr
    columns = duckdb.sql(f"DESCRIBE SELECT * FROM read_parquet('{out}')").fetchall()
    assert [(name, kind) for name, kind, *_ in columns] == [
        ("contig", "VARCHAR"), ("indel_type", "VARCHAR"), ("size", "INTEGER"), ("context", "VARCHAR"),
        ("tissue", "VARCHAR"), ("count", "BIGINT")]


def test_a_contig_without_opportunities_gets_an_empty_cube_with_the_columns(tmp_path):
    out = tmp_path / "c1.parquet"
    done = opportunity_tables("sum", out, window(tmp_path / "a.tsv.gz", []), window(tmp_path / "b.tsv.gz", []))
    assert done.returncode == 0, done.stderr
    assert rows(out) == []
    columns = duckdb.sql(f"DESCRIBE SELECT * FROM read_parquet('{out}')").fetchall()
    assert [(name, kind) for name, kind, *_ in columns] == [
        ("contig", "VARCHAR"), ("indel_type", "VARCHAR"), ("size", "INTEGER"), ("context", "VARCHAR"),
        ("tissue", "VARCHAR"), ("count", "BIGINT")]


def test_a_samples_contigs_are_joined_with_its_species_and_sample(tmp_path):
    contigs = []
    for contig, count in (("c1", 3), ("c2", 4)):
        contigs.append(tmp_path / f"{contig}.parquet")
        done = opportunity_tables("sum", contigs[-1],
                                  window(tmp_path / f"{contig}.tsv.gz", [(contig, "insertion", 1, "A", "liver", count)]))
        assert done.returncode == 0, done.stderr
    out = tmp_path / "merged.parquet"
    done = opportunity_tables("join", out, "--species", "ape", "--sample", "s1", *contigs)
    assert done.returncode == 0, done.stderr
    assert rows(out) == [("c1", "insertion", 1, "A", "liver", 3, "ape", "s1"),
                         ("c2", "insertion", 1, "A", "liver", 4, "ape", "s1")]


def test_a_merge_that_fails_leaves_no_file(tmp_path):
    # A window cut short, as a job killed while writing it leaves it.
    content = gzip.compress((HEADER + "c1\tinsertion\t1\tA\tliver\t5\n" * 20000).encode())
    cut = tmp_path / "cut.tsv.gz"
    cut.write_bytes(content[:len(content) // 2])
    out = tmp_path / "c1.parquet"
    done = opportunity_tables("sum", out, window(tmp_path / "a.tsv.gz", [("c1", "insertion", 1, "A", "liver", 5)]), cut)
    assert done.returncode != 0
    assert not out.exists()


def test_windows_with_different_columns_fail_the_sum_and_leave_no_file(tmp_path):
    # One window with a methylation bin where the other has a tissue.
    other = tmp_path / "b.tsv.gz"
    other.write_bytes(gzip.compress(b"contig\tindel_type\tsize\tcontext\tmethylation_bin\tcount\n"
                                    b"c1\tinsertion\t1\tA\t12\t5\n"))
    out = tmp_path / "c1.parquet"
    done = opportunity_tables("sum", out, window(tmp_path / "a.tsv.gz", [("c1", "insertion", 1, "A", "liver", 5)]), other)
    assert done.returncode != 0
    assert not list(tmp_path.glob("c1.parquet*"))


def test_a_cube_without_its_windows_total_is_not_written(tmp_path):
    # Only DuckDB misreading the windows makes such a cube (ticket 41), so the
    # check is given one directly.
    con = duckdb.connect()
    with pytest.raises(ValueError):
        write(con, "SELECT 'c1' AS contig, 4::BIGINT AS count", tmp_path / "c1.parquet", check=holding(con, 5))
    assert not list(tmp_path.glob("c1.parquet*"))


def test_a_cube_with_a_row_without_a_count_is_not_written(tmp_path):
    con = duckdb.connect()
    cube = "SELECT * FROM (VALUES ('c1', 5::BIGINT), ('c2', NULL)) AS t(contig, count)"
    with pytest.raises(ValueError):
        write(con, cube, tmp_path / "c1.parquet", check=holding(con, 5))
    assert not list(tmp_path.glob("c1.parquet*"))


def test_a_key_with_a_nul_byte_fails_the_sum_and_leaves_no_file(tmp_path):
    # In the pilot, DuckDB 2.0 made such a row from clean windows: a key cut
    # short by a NUL byte, its count lost (ticket 41).
    out = tmp_path / "c1.parquet"
    done = opportunity_tables("sum", out, window(tmp_path / "a.tsv.gz", [("c1", "insertion", 1, "A\0", "liver", 5)]))
    assert done.returncode != 0
    assert not list(tmp_path.glob("c1.parquet*"))


def test_a_count_that_isnt_a_whole_number_fails_the_sum_and_leaves_no_file(tmp_path):
    # DuckDB alone would read \N as no count, and 5.0 as 5.
    for count in ("\\N", "5.0"):
        out = tmp_path / "c1.parquet"
        done = opportunity_tables("sum", out, window(tmp_path / "a.tsv.gz", [("c1", "insertion", 1, "A", "liver", 5),
                                                                            ("c1", "insertion", 1, "A", "liver", count)]))
        assert done.returncode != 0, count
        assert not list(tmp_path.glob("c1.parquet*")), count


def test_workflow_4s_period_is_a_number_too(tmp_path):
    windows = tmp_path / "a.tsv.gz"
    windows.write_bytes(gzip.compress(b"contig\tcanonical_motif\tk\ttract_copies_bin\tcount\nc1\tAC\t2\t5-10\t4\n"))
    out = tmp_path / "c1.parquet"
    done = opportunity_tables("sum", out, windows)
    assert done.returncode == 0, done.stderr
    assert rows(out) == [("c1", "AC", 2, "5-10", 4)]


def test_a_join_that_fails_leaves_no_file(tmp_path):
    good = tmp_path / "c1.parquet"
    done = opportunity_tables("sum", good, window(tmp_path / "c1.tsv.gz", [("c1", "insertion", 1, "A", "liver", 3)]))
    assert done.returncode == 0, done.stderr
    cut = tmp_path / "c2.parquet"
    cut.write_bytes(good.read_bytes()[:len(good.read_bytes()) // 2])   # a contig's file cut short
    out = tmp_path / "merged.parquet"
    done = opportunity_tables("join", out, "--species", "ape", "--sample", "s1", good, cut)
    assert done.returncode != 0
    assert not out.exists()
