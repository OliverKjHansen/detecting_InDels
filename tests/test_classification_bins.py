"""The germline/somatic bin edges live in two scripts; this checks they agree.

join_read_classification.sh puts a bin on every call, and
classification_read_counts.sh counts the sample's reads per bin so those calls
have a denominator. Each carries its own copy of the edges, so that either can
run without the other. A disagreement between them would silently misalign
numerator and denominator, which is what these tests are here to catch.
"""

import random
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JOIN = REPO / "scripts" / "join_read_classification.sh"
COUNTS = REPO / "scripts" / "classification_read_counts.sh"

GERMLINE = ["Sperm.Sperm", "Sperm.Spermatocytes"]
SOMATIC = ["Blood.Monocytes", "Bladder.SmoothMuscle", "Breast.BasalEpithelial",
           "Lung.AlveolarMacrophages", "Colon.RightEndocrine", "Thyroid.Epithelial"]

# The edges as documented in both scripts: bin(d) is the first value whose edge
# d falls below, and 12 above the last.
EDGES = [(-10, -12), (-8, -10), (-6, -8), (-4, -6), (-2, -4), (-1, -2), (0, -1),
         (1, 0), (2, 1), (4, 2), (6, 4), (8, 6), (10, 8)]


def bin_of(diff):
    """The bin a difference belongs in, independently of either script."""
    for edge, value in EDGES:
        if diff < edge:
            return value
    return 12


def write_likelihoods(path, diffs, cpgs):
    """A likelihood table whose reads have exactly the given differences.

    The somatic columns are fixed and one germline column carries the
    difference, so the score each script computes is known in advance.
    """
    columns = ["readID", "chr", "start"] + GERMLINE + SOMATIC + ["CpGs_onRead", "warning"]
    with open(path, "w") as out:
        out.write("\t".join(columns) + "\n")
        for i, (diff, cpg) in enumerate(zip(diffs, cpgs)):
            # max somatic is -20, so the germline maximum has to be -20 + diff
            germline = [f"{-20 + diff:.6f}", f"{-20 + diff - 5:.6f}"]
            somatic = ["-20.000000", "-25.000000", "-30.000000",
                       "-21.000000", "-40.000000", "-50.000000"]
            out.write("\t".join([f"read{i}", "c1", "1"] + germline + somatic
                                + [str(cpg), "NA"]) + "\n")


def write_calls(path, n_reads):
    with open(path, "w") as out:
        out.write("contig\tposition\tindel_string\treadID\n")
        for i in range(n_reads):
            out.write(f"c1\t1\tA\tread{i}\n")


def counts_from_script(tmp_path, likelihoods):
    out = tmp_path / "counts.tsv"
    subprocess.run(["bash", str(COUNTS), str(likelihoods), str(out)],
                   check=True, capture_output=True)
    rows = out.read_text().splitlines()[1:]
    return {int(f[0]): int(f[1]) for f in (row.split("\t") for row in rows)}


def counts_from_join(tmp_path, likelihoods, n_reads):
    calls = tmp_path / "calls.tsv"
    joined = tmp_path / "joined.tsv"
    write_calls(calls, n_reads)
    subprocess.run(["bash", str(JOIN), str(calls), str(likelihoods), str(joined)],
                   check=True, capture_output=True)
    lines = joined.read_text().splitlines()
    header = lines[0].split("\t")
    column = header.index("classification_bin")
    counts = {}
    for line in lines[1:]:
        value = line.split("\t")[column]
        counts[int(value)] = counts.get(int(value), 0) + 1
    return counts


# Differences that land on and just around every edge, where an off-by-one in
# either script would show.
BOUNDARY_DIFFS = [d + offset
                  for d, _ in EDGES
                  for offset in (-0.001, 0, 0.001)] + [-50, -11, 10, 10.001, 50]


def test_bins_agree_on_boundaries(tmp_path):
    likelihoods = tmp_path / "likelihoods.tsv"
    write_likelihoods(likelihoods, BOUNDARY_DIFFS, [10] * len(BOUNDARY_DIFFS))

    from_counts = counts_from_script(tmp_path, likelihoods)
    from_join = counts_from_join(tmp_path, likelihoods, len(BOUNDARY_DIFFS))

    assert from_counts == from_join


def test_bins_match_the_documented_edges(tmp_path):
    likelihoods = tmp_path / "likelihoods.tsv"
    write_likelihoods(likelihoods, BOUNDARY_DIFFS, [10] * len(BOUNDARY_DIFFS))

    expected = {}
    for diff in BOUNDARY_DIFFS:
        expected[bin_of(diff)] = expected.get(bin_of(diff), 0) + 1

    assert counts_from_script(tmp_path, likelihoods) == expected


def test_counts_and_cpgs_are_totals(tmp_path):
    random.seed(11)
    diffs = [random.uniform(-30, 30) for _ in range(500)]
    cpgs = [random.randint(0, 300) for _ in range(500)]
    likelihoods = tmp_path / "likelihoods.tsv"
    write_likelihoods(likelihoods, diffs, cpgs)

    out = tmp_path / "counts.tsv"
    subprocess.run(["bash", str(COUNTS), str(likelihoods), str(out)],
                   check=True, capture_output=True)
    rows = [row.split("\t") for row in out.read_text().splitlines()[1:]]

    assert sum(int(r[1]) for r in rows) == len(diffs)
    assert sum(int(r[2]) for r in rows) == sum(cpgs)
    # every read is counted under its own bin, with its own CpGs
    for row in rows:
        bin_value = int(row[0])
        mine = [c for d, c in zip(diffs, cpgs) if bin_of(d) == bin_value]
        assert int(row[1]) == len(mine)
        assert int(row[2]) == sum(mine)


def test_bins_are_in_order(tmp_path):
    random.seed(12)
    diffs = [random.uniform(-30, 30) for _ in range(200)]
    likelihoods = tmp_path / "likelihoods.tsv"
    write_likelihoods(likelihoods, diffs, [1] * len(diffs))

    out = tmp_path / "counts.tsv"
    subprocess.run(["bash", str(COUNTS), str(likelihoods), str(out)],
                   check=True, capture_output=True)
    bins = [int(row.split("\t")[0]) for row in out.read_text().splitlines()[1:]]

    assert bins == sorted(bins)


def test_missing_input_fails_loudly(tmp_path):
    done = subprocess.run(["bash", str(COUNTS), str(tmp_path / "nope.tsv"),
                           str(tmp_path / "out.tsv")], capture_output=True)
    assert done.returncode != 0
    assert not (tmp_path / "out.tsv").exists()


def test_unknown_column_names_fail(tmp_path):
    """A likelihood table without the germline columns must not yield empty counts."""
    path = tmp_path / "likelihoods.tsv"
    path.write_text("readID\tchr\tstart\tOther.Cells\tCpGs_onRead\n"
                    "read0\tc1\t1\t-20\t10\n")
    done = subprocess.run(["bash", str(COUNTS), str(path), str(tmp_path / "out.tsv")],
                          capture_output=True)
    assert done.returncode != 0
