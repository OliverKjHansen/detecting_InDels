"""
How the mirror replay compares a local table with the cluster's: by content,
whatever the row order, the number format, or the columns only the cluster's
workflow appends.
"""

from compare_tables import compare


def write(path, text):
    path.write_text(text)
    return str(path)


def test_tables_in_a_different_row_order_are_the_same(tmp_path):
    local = write(tmp_path / "local.tsv", "#contig\tposition\tkmerQ\nc1\t10\t93\nc1\t20\t40\n")
    cluster = write(tmp_path / "cluster.tsv", "#contig\tposition\tkmerQ\nc1\t20\t40\nc1\t10\t93\n")
    assert compare(local, cluster).same


def test_numbers_written_in_a_different_format_are_the_same(tmp_path):
    # pandas writes an integer column holding NA as floats, and purities may
    # be printed to more or fewer decimals.
    local = write(tmp_path / "local.tsv", "np\tpurity\tsx\n20.0\t100.00\tNA\n7\t0.5\t1e-05\n")
    cluster = write(tmp_path / "cluster.tsv", "np\tpurity\tsx\n20\t100\tNA\n7.0\t0.50\t0.00001\n")
    assert compare(local, cluster).same


def test_columns_only_the_cluster_appends_are_ignored(tmp_path):
    # The workflows append species and sample after the scripts have run.
    local = write(tmp_path / "local.tsv", "contig\tposition\nc1\t10\n")
    cluster = write(tmp_path / "cluster.tsv", "contig\tposition\tspecies\tsample\nc1\t10\thuman\tchk\n")
    assert compare(local, cluster).same


def test_a_changed_value_is_a_difference_shown_on_both_sides(tmp_path):
    local = write(tmp_path / "local.tsv", "position\tsx\n10\t0\n20\t0\n")
    cluster = write(tmp_path / "cluster.tsv", "position\tsx\n10\t0\n20\t3\n")
    comparison = compare(local, cluster)
    assert not comparison.same
    assert (list(comparison.only_local), list(comparison.only_cluster)) == ([("20", "0")], [("20", "3")])


def test_a_row_written_twice_is_a_difference(tmp_path):
    local = write(tmp_path / "local.tsv", "position\n10\n10\n")
    cluster = write(tmp_path / "cluster.tsv", "position\n10\n")
    assert not compare(local, cluster).same


def test_a_column_only_the_local_table_has_is_a_difference(tmp_path):
    local = write(tmp_path / "local.tsv", "position\tstrand\n10\t+\n")
    cluster = write(tmp_path / "cluster.tsv", "position\n10\n")
    comparison = compare(local, cluster)
    assert not comparison.same
    assert comparison.local_columns == ["strand"]


def test_any_other_column_only_the_cluster_has_is_a_difference(tmp_path):
    # Local code that stops writing a column must not pass as "same".
    local = write(tmp_path / "local.tsv", "position\n10\n")
    cluster = write(tmp_path / "cluster.tsv", "position\tstrand\tspecies\n10\t+\thuman\n")
    comparison = compare(local, cluster)
    assert not comparison.same
    assert comparison.cluster_columns == ["strand"]


def test_large_integers_are_compared_exactly(tmp_path):
    local = write(tmp_path / "local.tsv", "count\n1000000001\n")
    cluster = write(tmp_path / "cluster.tsv", "count\n1000000002\n")
    assert not compare(local, cluster).same


def test_opportunity_rows_counting_nothing_are_ignored(tmp_path):
    # The scans once wrote a row for every kmerQ level, counted or not.
    local = write(tmp_path / "local.tsv", "context\tkmerQ\tcount\nAC\t93\t5\n")
    cluster = write(tmp_path / "cluster.tsv", "context\tkmerQ\tcount\nAC\t93\t5\nAC\t92\t0\nAC\t91\t0\n")
    assert compare(local, cluster).same
