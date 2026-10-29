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


def test_the_cube_keeps_the_scans_columns_with_size_and_count_as_numbers(tmp_path):
    out = tmp_path / "c1.parquet"
    done = opportunity_tables("sum", out, window(tmp_path / "a.tsv.gz", [("c1", "insertion", 1, "A", "liver", 5)]))
    assert done.returncode == 0, done.stderr
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
