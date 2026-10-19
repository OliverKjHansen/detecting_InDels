#!/usr/bin/env python3
"""
Generate the miniature dataset used by the local profile and the tests.

Everything is synthetic and deterministic, laid out to match the `layout`
templates in config.yaml so the `local` profile needs no special-casing:

    tests/fixtures/
      fasta/testspecies/testsample/diploid_assembly/testsample_diploid.fa[.fai]
      fasta/testspecies/testsample/haploid_assembly/testsample_haploid.fa[.fai]
      bam/testspecies/testsample/kinetics/testsample_kinetics_diploid.bam[.bai]
      coverage/testsample_avg_coverage.txt
      tissue/testspecies/testsample_sorted_tissue_classification.tsv
      microsat/testspecies/testsample/testcontig1_trs_tidy.tsv
      microsat/testspecies/testsample/testcontig1_unresolved_haplotypes.bed
      microsat/testspecies/testsample/testcontig1_filters.tsv   <- pre-merged mask
      tr_bed/testcontig1_M1.bed
      truth.tsv

The reference is built from four regions with known properties, and the reads
carry InDels planted at known positions, so the tests can assert an exact call
set rather than "something was produced".

    python tests/make_fixtures.py
"""

import array
import os
import random
import sys

import pysam

CONTIG = "testcontig1"
SPECIES = "testspecies"
SAMPLE = "testsample"

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")

# Reference layout ----------------------------------------------------------
#   [   0, 2000)  "ACGT" repeat  - complex, no homopolymers: InDels here are called
#   [2000, 2100)  poly-A         - InDels here are dropped by the homopolymer filter
#   [2100, 2300)  "CA" repeat    - microsatellite: masked out of workflow B,
#                                  and the target region for workflow C
#   [2300, 6000)  pseudo-random  - filler
COMPLEX_END = 2000
HOMOPOLYMER_START, HOMOPOLYMER_END = 2000, 2100
MICROSAT_START, MICROSAT_END = 2100, 2300
CONTIG_LENGTH = 6000

UNRESOLVED_START, UNRESOLVED_END = 5500, 5600

N_CLEAN_READS = 16
MAPPING_QUALITY = 60
BASE_QUALITY = 93  # PacBio HiFi maximum; keeps kmerQ above every threshold

# Planted InDels: (read_name, kind, position, length, expectation)
#
# `position` is the 0-based reference index of the first affected base. For a
# deletion that is the first deleted base; for an insertion the new bases go
# immediately after it.
PLANTED = [
    ("read_del1",     "del", 500,  1, "called by workflow B"),
    ("read_ins1",     "ins", 1003, 1, "called by workflow B"),
    ("read_del2",     "del", 1500, 2, "called by workflow B"),
    ("read_del_homo", "del", 2050, 1, "dropped: homopolymer context"),
    ("read_tr_del1",  "del", 2200, 1, "masked from B; not a multiple of k=2 in C"),
    ("read_tr_del2",  "del", 2204, 2, "masked from B; called by workflow C"),
]

TISSUES = ["liver", "brain"]


def build_reference():
    """Deterministic reference sequence with the four regions described above."""
    rng = random.Random(42)

    complex_region = ("ACGT" * ((COMPLEX_END // 4) + 1))[:COMPLEX_END]
    homopolymer = "A" * (HOMOPOLYMER_END - HOMOPOLYMER_START)
    microsat = ("CA" * (((MICROSAT_END - MICROSAT_START) // 2) + 1))[
        : MICROSAT_END - MICROSAT_START
    ]
    filler_length = CONTIG_LENGTH - MICROSAT_END
    filler = "".join(rng.choice("ACGT") for _ in range(filler_length))

    reference = complex_region + homopolymer + microsat + filler
    assert len(reference) == CONTIG_LENGTH, len(reference)
    return reference


def write_fasta(path, contig, sequence, width=60):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as handle:
        handle.write(f">{contig}\n")
        for i in range(0, len(sequence), width):
            handle.write(sequence[i : i + width] + "\n")
    pysam.faidx(path)


def make_read(name, reference, kind=None, position=None, length=0):
    """
    Build one AlignedSegment spanning the whole contig.

    Clean reads are a perfect match. Mutant reads carry a single InDel, encoded
    both in the sequence and in the CIGAR so pysam's pileup reports it.
    """
    read = pysam.AlignedSegment()

    if kind is None:
        sequence = reference
        cigar = [(0, len(reference))]
    elif kind == "del":
        sequence = reference[:position] + reference[position + length :]
        cigar = [
            (0, position),
            (2, length),
            (0, len(reference) - position - length),
        ]
    elif kind == "ins":
        # Choose an inserted base differing from both neighbours so the
        # insertion does not land in a homopolymer run and get filtered.
        left, right = reference[position], reference[position + 1]
        inserted = next(b for b in "ACGT" if b != left and b != right)
        sequence = (
            reference[: position + 1] + inserted * length + reference[position + 1 :]
        )
        cigar = [
            (0, position + 1),
            (1, length),
            (0, len(reference) - position - 1),
        ]
    else:
        raise ValueError(f"unknown kind: {kind}")

    read.query_name = name
    read.query_sequence = sequence
    read.flag = 0
    read.reference_id = 0
    read.reference_start = 0
    read.mapping_quality = MAPPING_QUALITY
    read.cigartuples = cigar
    # Must be assigned after query_sequence, which resets the qualities.
    read.query_qualities = array.array("B", [BASE_QUALITY] * len(sequence))

    # PacBio consensus tags the callers read.
    read.set_tag("np", 20, value_type="i")
    read.set_tag("ec", 20.0, value_type="f")
    read.set_tag("sm", array.array("B", [BASE_QUALITY] * len(sequence)))
    read.set_tag("sx", array.array("B", [BASE_QUALITY] * len(sequence)))
    # subread_array() reads this as (length, coverage) pairs.
    read.set_tag("sa", array.array("I", [len(sequence), 10]))

    return read


def write_bam(path, reference, reads):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": len(reference)}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for read in reads:
            out.write(read)
    pysam.index(path)


def write_lines(path, lines):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as handle:
        for line in lines:
            handle.write(line + "\n")


def main():
    reference = build_reference()

    # -- assemblies --------------------------------------------------------
    diploid = os.path.join(
        FIXTURES, "fasta", SPECIES, SAMPLE, "diploid_assembly", f"{SAMPLE}_diploid.fa"
    )
    haploid = os.path.join(
        FIXTURES, "fasta", SPECIES, SAMPLE, "haploid_assembly", f"{SAMPLE}_haploid.fa"
    )
    write_fasta(diploid, CONTIG, reference)
    write_fasta(haploid, CONTIG, reference)

    # -- reads -------------------------------------------------------------
    reads = [
        make_read(name, reference, kind, position, length)
        for name, kind, position, length, _ in PLANTED
    ]
    reads += [
        make_read(f"read_clean{i:02d}", reference) for i in range(N_CLEAN_READS)
    ]

    bam = os.path.join(
        FIXTURES, "bam", SPECIES, SAMPLE, "kinetics", f"{SAMPLE}_kinetics_diploid.bam"
    )
    write_bam(bam, reference, reads)

    # -- coverage ----------------------------------------------------------
    # Every read spans the whole contig, so depth equals the read count. The
    # callers keep positions within 0.5x-2x of this value.
    write_lines(
        os.path.join(FIXTURES, "coverage", f"{SAMPLE}_avg_coverage.txt"),
        [str(len(reads))],
    )

    # -- tissue classification --------------------------------------------
    write_lines(
        os.path.join(FIXTURES, "tissue", SPECIES, f"{SAMPLE}_sorted_tissue_classification.tsv"),
        [f"{read.query_name}\t{TISSUES[i % len(TISSUES)]}" for i, read in enumerate(reads)],
    )

    # -- masks -------------------------------------------------------------
    microsat_dir = os.path.join(FIXTURES, "microsat", SPECIES, SAMPLE)

    # Shaped like the real *_trs_tidy.tsv: the workflow only reads columns 1-3.
    write_lines(
        os.path.join(microsat_dir, f"{CONTIG}_trs_tidy.tsv"),
        [f"{CONTIG}\t{MICROSAT_START}\t{MICROSAT_END}\tCA\t2\t100"],
    )
    write_lines(
        os.path.join(microsat_dir, f"{CONTIG}_unresolved_haplotypes.bed"),
        [f"{CONTIG}\t{UNRESOLVED_START}\t{UNRESOLVED_END}"],
    )
    # Pre-merged union of the two, equivalent to what `bedtools merge` produces
    # in step 1. Committing it lets the tests exercise the callers without
    # bedtools installed locally.
    write_lines(
        os.path.join(microsat_dir, f"{CONTIG}_filters.tsv"),
        [
            f"{CONTIG}\t{MICROSAT_START}\t{MICROSAT_END}",
            f"{CONTIG}\t{UNRESOLVED_START}\t{UNRESOLVED_END}",
        ],
    )

    # -- tandem repeat BED for workflow C ---------------------------------
    # Columns per load_tr_bed_with_tree(): contig, start, end, canonical_motif,
    # k, copies, trf_purity, model_purity, model. First line is a header.
    write_lines(
        os.path.join(FIXTURES, "tr_bed", f"{CONTIG}_M1.bed"),
        [
            "contig\tstart\tend\tcanonical_motif\tk\tcopies\ttrf_purity\tmodel_purity\tmodel",
            f"{CONTIG}\t{MICROSAT_START}\t{MICROSAT_END}\tAC\t2\t100\t100.00\t100.00\tM1",
        ],
    )

    # -- truth manifest ----------------------------------------------------
    write_lines(
        os.path.join(FIXTURES, "truth.tsv"),
        ["read\tkind\tposition\tlength\texpectation"]
        + [f"{n}\t{k}\t{p}\t{l}\t{e}" for n, k, p, l, e in PLANTED],
    )

    print(f"Fixtures written to {FIXTURES}", file=sys.stderr)
    print(f"  contig {CONTIG}: {len(reference)} bp, {len(reads)} reads", file=sys.stderr)
    for name, kind, position, length, expectation in PLANTED:
        print(f"  {name:15s} {kind} {length}bp @ {position:5d}  {expectation}", file=sys.stderr)


if __name__ == "__main__":
    main()
