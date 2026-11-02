"""
The compact read bins (ticket 44): scripts/read_bins.py writes a sample's read
classification table once more as a key per read name and a one-byte bin,
which every scan loads in place of the whole table.
"""

import os
import subprocess
import sys

import make_fixtures as mf
import targets
from read_bins import load

SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts", "read_bins.py")


def table(path, bins):
    mf.write_read_classification(str(path), bins)
    return path


def read_bins(*arguments):
    return subprocess.run([sys.executable, SCRIPT, *map(str, arguments)], capture_output=True, text=True)


def test_every_read_gets_its_bin_from_the_compact_table(tmp_path):
    names = [f"m84011_220902_175841_s1/{zmw}/ccs" for zmw in range(1000)]
    labels = ["-12", "-1", "0", "4", "12"]
    rows = [(name, labels[i % len(labels)]) for i, name in enumerate(names)]
    out = tmp_path / "read_bins.npz"
    done = read_bins(table(tmp_path / "read_classification.tsv", rows), out)
    assert done.returncode == 0, done.stderr
    bins = load(str(out))
    assert {name: bins.get(name, "NA") for name in names} == dict(rows)
    assert bins.get("m84011_220902_175841_s1/1000/ccs", "NA") == "NA"
    # An 8-byte key and a 1-byte bin per read, and a little for the file's own layout.
    assert out.stat().st_size < 9 * len(names) + 1000


def test_a_read_listed_twice_keeps_its_last_bin_as_the_text_table_did(tmp_path):
    out = tmp_path / "read_bins.npz"
    done = read_bins(table(tmp_path / "read_classification.tsv", [("r1", "4"), ("r2", "0"), ("r1", "-6")]), out)
    assert done.returncode == 0, done.stderr
    bins = load(str(out))
    assert (bins.get("r1", "NA"), bins.get("r2", "NA")) == ("-6", "0")


def test_more_bins_than_a_byte_holds_fail_and_leave_no_file(tmp_path):
    out = tmp_path / "read_bins.npz"
    done = read_bins(table(tmp_path / "read_classification.tsv", [(f"r{i}", str(i)) for i in range(257)]), out)
    assert done.returncode != 0
    assert not list(tmp_path.glob("read_bins.npz*"))


def test_a_sample_without_a_table_gives_every_read_na():
    assert load("NA").get("r1", "NA") == "NA"


def test_workflow_1_writes_a_methylation_samples_read_bins(tmp_path, monkeypatch):
    sample = targets.sample(tmp_path, monkeypatch, "ACGT" * 50)
    sample.configure("methylation_samples", **{targets.SPECIES: [targets.SAMPLE]})
    mf.write_read_classification(sample.path("read_classification"), [("r1", "4"), ("r2", "-6")])
    sample.run("precompute", f"read_bins_{targets.SPECIES}_{targets.SAMPLE}_diploid")
    bins = load(sample.path("read_bins"))
    assert [bins.get(name, "NA") for name in ("r1", "r2", "r3")] == ["4", "-6", "NA"]
