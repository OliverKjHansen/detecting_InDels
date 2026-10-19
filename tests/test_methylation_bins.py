"""
Methylation bins (ticket 15): for samples with per-read cell-type likelihoods,
calls and opportunities carry each read's germline/somatic bin, all 14 of
them, so per-bin denominators are counted rather than approximated. Reads
without a classification, and samples without likelihoods, get NA.
"""

import os
import subprocess

import make_fixtures as mf
from scans import nonrepetitive_scan, rows_of, tandem_repeat_scan

REFERENCE = mf.build_reference()
SCRIPTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")


def dataset(tmp_path, reads, bins):
    """A dataset whose read classification table gives these reads these bins."""
    clean = [mf.make_read(f"clean{i:02d}", REFERENCE) for i in range(4)]
    paths = mf.write_dataset(str(tmp_path), REFERENCE, reads + clean)
    paths["classification"] = str(tmp_path / "read_classification.tsv")
    with open(paths["classification"], "w") as out:
        out.write("readID\tmax_germline\tmax_somatic\tclassification_diff\tclassification_bin\tCpGs_onRead\n")
        for name, label in bins.items():
            out.write(f"{name}\t-20\t-21\t1\t{label}\t12\n")
    return paths


def bins_of(table):
    totals = {}
    for row in table:
        totals[row["methylation_bin"]] = totals.get(row["methylation_bin"], 0) + int(row.get("count", 1))
    return totals


def test_calls_carry_their_reads_bin_or_na(tmp_path):
    reads = [mf.make_read("germ", REFERENCE, "del", 500, 1), mf.make_read("soma", REFERENCE, "ins", 1003, 1),
             mf.make_read("none", REFERENCE, "del", 1500, 2)]
    calls, _, _ = nonrepetitive_scan(dataset(tmp_path, reads, {"germ": "4", "soma": "-6"}))
    by_read = {row["readID"]: row["methylation_bin"] for row in rows_of(calls)}
    assert by_read == {"germ": "4", "soma": "-6", "none": "NA"}


def test_per_bin_opportunities_add_up_to_the_samples_total(tmp_path):
    reads = [mf.make_read("germ", REFERENCE), mf.make_read("soma", REFERENCE)]
    paths = dataset(tmp_path, reads, {"germ": "4", "soma": "-6"})
    _, opportunities, _ = nonrepetitive_scan(paths)
    per_bin = bins_of(rows_of(opportunities))
    # Six reads, all clean copies of the contig: each counts the same positions.
    assert set(per_bin) == {"4", "-6", "NA"}
    assert per_bin["4"] == per_bin["-6"] == per_bin["NA"] / 4
    _, plain, _ = nonrepetitive_scan({**paths, "classification": "NA"})
    assert sum(per_bin.values()) == sum(int(row["count"]) for row in rows_of(plain))
    assert set(bins_of(rows_of(plain))) == {"NA"}


def test_the_repeat_scan_carries_the_bin_too(tmp_path):
    bed = str(tmp_path / "tracts.bed")
    mf.write_tr_bed(bed, [(2100, 2300, "AC", 2)])
    paths = dataset(tmp_path, [mf.make_read("germ", REFERENCE, "del", 2204, 2)], {"germ": "8"})
    calls, opportunities, _ = tandem_repeat_scan({"M1": bed}, paths)
    assert [row["methylation_bin"] for row in rows_of(calls["M1"])] == ["8"]
    assert set(bins_of(rows_of(opportunities["M1"]))) == {"8", "NA"}


def write_likelihoods(path, diffs):
    """A likelihood table in which read{i} scores diffs[i]: best germline minus best somatic."""
    with open(path, "w") as out:
        out.write("readID\tchr\tstart\tSperm.Sperm\tBlood.Monocytes\tCpGs_onRead\twarning\n")
        for i, diff in enumerate(diffs):
            out.write(f"read{i}\tc1\t1\t{-20 + diff}\t-20\t{10 + i}\tNA\n")


def test_the_per_read_table_bins_reads_as_the_classification_scripts_do(tmp_path):
    likelihoods = tmp_path / "likelihoods.tsv"
    write_likelihoods(likelihoods, [-11, -3, 0.5, 9.9, 15])
    arguments = ["Sperm.", "Blood.Monocytes", "-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10",
                 "-12,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,12"]
    table, counts = tmp_path / "reads.tsv", tmp_path / "counts.tsv"
    subprocess.run(["bash", os.path.join(SCRIPTS, "read_classification.sh"), str(likelihoods), str(table),
                    *arguments], check=True, capture_output=True)
    subprocess.run(["bash", os.path.join(SCRIPTS, "classification_read_counts.sh"), str(likelihoods), str(counts),
                    *arguments], check=True, capture_output=True)
    header, *rows = table.read_text().splitlines()
    assert header.split("\t") == ["readID", "max_germline", "max_somatic", "classification_diff",
                                  "classification_bin", "CpGs_onRead"]
    bins = {row.split("\t")[0]: row.split("\t")[4] for row in rows}
    assert bins == {"read0": "-12", "read1": "-4", "read2": "0", "read3": "8", "read4": "12"}
    per_bin = {line.split("\t")[0]: line.split("\t")[1] for line in counts.read_text().splitlines()[1:]}
    assert per_bin == {label: "1" for label in bins.values()}
