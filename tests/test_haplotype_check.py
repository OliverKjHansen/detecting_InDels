"""
The hap1/hap2 mismapping test (ticket 19). Each call-carrying read is aligned
to hap1 alone and to hap2 alone, and each call records whether its InDel is in
each alignment. A mutation is in the read, so it shows against both; a
germline difference between the haplotypes, on a read aligned to the wrong
haplotype in the diploid assembly, disappears against the read's own.

pbmm2 isn't run here: the fixture builder writes each read's alignment to
each haplotype, as the aligner would.
"""

import os
import subprocess
import sys

import pytest

import make_fixtures as mf
import targets
from scans import nonrepetitive_scan, rows_of

SCRIPTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")

HAP1 = mf.build_reference()
GERMLINE = 500                           # hap2 lacks hap1's base 500
HAP2 = HAP1[:GERMLINE] + HAP1[GERMLINE + 1:]
MUTATION = 1500                          # a 2 bp deletion of hap1's bases 1500-1501 in one read


def reads(reverse):
    """
    Per read, its alignment to the diploid assembly (where both reads sit on
    hap1's contig), to hap1 alone and to hap2 alone: the same read each time.
    """
    from_hap2 = [mf.make_read("from_hap2", HAP1, "del", GERMLINE, 1, reverse=reverse),
                 mf.make_read("from_hap2", HAP1, "del", GERMLINE, 1, reverse=reverse),
                 mf.make_read("from_hap2", HAP2, reverse=reverse)]
    # The mutated read starts after the germline difference, so against hap2
    # it differs only by the mutation.
    mutated = [mf.make_read("mutated", HAP1, "del", MUTATION, 2, start=1000, reverse=reverse),
               mf.make_read("mutated", HAP1, "del", MUTATION, 2, start=1000, reverse=reverse),
               mf.make_read("mutated", HAP2, "del", MUTATION - 1, 2, start=999, reverse=reverse)]
    for alignments in (from_hap2, mutated):
        assert len({alignment.query_sequence for alignment in alignments}) == 1
    return from_hap2, mutated


def write_records(path, alignments, contig):
    """The alignments as the workflows' locations.tsv: read, flag, contig, position, mapQ, CIGAR."""
    with open(path, "w") as out:
        for a in alignments:
            out.write(f"{a.query_name}\t{a.flag}\t{contig}\t{a.reference_start + 1}\t60\t{a.cigarstring}\n")
    return path


def checked_calls(tmp_path, reverse):
    from_hap2, mutated = reads(reverse)
    diploid = [from_hap2[0], mutated[0]] + [mf.make_read(f"clean{i:02d}", HAP1) for i in range(4)]
    calls, _, _ = nonrepetitive_scan(mf.write_dataset(str(tmp_path), HAP1, diploid))
    (tmp_path / "calls.tsv").write_text(calls)
    with open(tmp_path / "reads.fastq", "w") as out:
        for a in (from_hap2[0], mutated[0]):
            out.write(f"@{a.query_name}\n{a.get_forward_sequence()}\n+\n{'~' * a.query_length}\n")
    haploid = write_records(tmp_path / "haploid.tsv", [from_hap2[1], mutated[1]], "hap")
    hap1 = write_records(tmp_path / "hap1.tsv", [from_hap2[1], mutated[1]], "h1tg000001l")
    hap2 = write_records(tmp_path / "hap2.tsv", [from_hap2[2], mutated[2]], "h2tg000001l")
    done = subprocess.run([sys.executable, os.path.join(SCRIPTS, "filtering_haploid_mappping.py"),
                           str(tmp_path / "calls.tsv"), str(haploid), str(tmp_path / "reads.fastq"),
                           "--haplotype", f"hap1={hap1}", "--haplotype", f"hap2={hap2}"],
                          capture_output=True, text=True, check=True)
    return {row["readID"]: (row["hap1_check"], row["hap2_check"]) for row in rows_of(done.stdout)}


@pytest.mark.parametrize("reverse", [False, True])
def test_a_mutation_shows_against_both_haplotypes_and_a_germline_difference_against_one(tmp_path, reverse):
    assert checked_calls(tmp_path, reverse) == {"mutated": ("True", "True"), "from_hap2": ("True", "False")}


def test_a_read_missing_from_a_haplotypes_alignment_is_unchecked_there(tmp_path):
    from_hap2, mutated = reads(False)
    checks = checked_calls(tmp_path, False)
    # Rewrite hap2's records without the mutated read and check again.
    write_records(tmp_path / "hap2.tsv", [from_hap2[2]], "h2tg000001l")
    done = subprocess.run([sys.executable, os.path.join(SCRIPTS, "filtering_haploid_mappping.py"),
                           str(tmp_path / "calls.tsv"), str(tmp_path / "haploid.tsv"), str(tmp_path / "reads.fastq"),
                           "--haplotype", f"hap1={tmp_path / 'hap1.tsv'}", "--haplotype", f"hap2={tmp_path / 'hap2.tsv'}"],
                          capture_output=True, text=True, check=True)
    after = {row["readID"]: (row["hap1_check"], row["hap2_check"]) for row in rows_of(done.stdout)}
    assert checks["mutated"] == ("True", "True")
    assert after == {"mutated": ("True", "NA"), "from_hap2": ("True", "False")}


def test_precompute_splits_the_diploid_assembly_into_its_haplotypes(tmp_path, monkeypatch):
    sample = targets.sample(tmp_path, monkeypatch, "ACGT")
    sample.write("assembly", [">h1tg000001l", "AAAA", "CC", ">h2tg000001l rd:i:30", "GGGG", ">h1tg000002l", "TTTT"])
    sample.run("precompute", f"haplotypes_{targets.SPECIES}_{targets.SAMPLE}_diploid")
    with open(sample.path("haplotype_assembly", haplotype="hap1")) as hap1, \
         open(sample.path("haplotype_assembly", haplotype="hap2")) as hap2:
        assert (hap1.read(), hap2.read()) == (">h1tg000001l\nAAAA\nCC\n>h1tg000002l\nTTTT\n",
                                              ">h2tg000001l rd:i:30\nGGGG\n")
