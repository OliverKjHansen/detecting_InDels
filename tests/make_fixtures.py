#!/usr/bin/env python3
"""
Generate the miniature dataset the tests run on.

conftest.py builds it into a temporary directory for every test session. Run
this file to write the same files to tests/fixtures/ to look at:

    python tests/make_fixtures.py

Everything is synthetic and deterministic, laid out like the cluster's inputs:

    fixtures/
      fasta/testspecies/testsample/diploid_assembly/testsample_diploid.fa[.fai]
      fasta/testspecies/testsample/haploid_assembly/testsample_haploid.fa[.fai]
      bam/testspecies/testsample/kinetics/testsample_kinetics_diploid.bam[.bai]
      coverage/testsample_avg_coverage.txt
      tissue/testspecies/testsample_sorted_tissue_classification.tsv
      microsat/testspecies/testsample/testcontig1_trs_tidy.tsv
      microsat/testspecies/testsample/testcontig1_unresolved_haplotypes.bed
      microsat/testspecies/testsample/testcontig1_filters.tsv   <- pre-merged mask
      tr_bed/testcontig1_M1.bed

The reference is built from four regions with known properties, and the reads
carry InDels planted at known positions (PLANTED), so the tests can assert an
exact call set rather than "something was produced".

Every fixture read is a forward-strand copy of the whole contig at mapQ 60 and
base quality 93. Tests that need the shapes real reads have build their own
dataset from the same pieces: make_read() for reverse-strand, partial,
soft-clipped or low-mapQ reads with chosen base qualities and subread tags,
write_tr_bed() for repeat tracts, which may overlap, and write_dataset() for
the BAM and the scan's other inputs.
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
SUBREADS = 10      # subreads covering each base, all of them matching

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


def make_read(name, reference, kind=None, position=None, length=0, *, start=0, end=None,
              reverse=False, mapq=MAPPING_QUALITY, soft_clip=(0, 0), substitutions=None,
              qualities=None, subreads=None):
    """
    One aligned read over reference[start:end], by default the whole contig,
    as pbmm2 writes it.

    Clean reads are a perfect match. A read can carry one InDel, encoded both
    in the sequence and in the CIGAR so pysam's pileup reports it: kind "del"
    or "ins", at `position`, the 0-based reference index of the first affected
    base. A deletion removes `length` bases from there; an insertion puts
    `length` new bases right after it. The anchor and the base after the gap
    must lie inside [start, end).

    reverse: aligned to the reverse strand (flag 16). As in a real BAM the
    sequence and base qualities are stored in reference orientation, but the
    per-base subread tags (sm, sx, sa) stay in the order the read was
    sequenced, so they run backwards along the alignment.

    mapq: the mapping quality.

    soft_clip: (left, right), the number of extra read bases soft-clipped off
    each end, in reference orientation. They match nothing in the contig.

    substitutions: the read's base at chosen reference positions, where it
    differs from the contig.

    qualities: base quality per reference position, with the key "inserted"
    for every inserted base. Other bases get BASE_QUALITY.

    subreads: (matching, mismatching, covering) subreads per reference
    position, with the key "inserted" for every inserted base. Other bases get
    SUBREADS covering subreads, all of them matching.
    """
    end = len(reference) if end is None else end
    if kind is not None:
        anchor = position if kind == "ins" else position - 1
        after = position + 1 if kind == "ins" else position + length
        if not start <= anchor < after < end:
            raise ValueError(f"{name}: the InDel at {position} is not inside the read [{start}, {end})")

    # The reference position of each read base, None for an inserted base.
    if kind is None:
        bases = list(range(start, end))
        cigar = [(0, end - start)]
    elif kind == "del":
        bases = list(range(start, position)) + list(range(position + length, end))
        cigar = [(0, position - start), (2, length), (0, end - position - length)]
    elif kind == "ins":
        bases = list(range(start, position + 1)) + [None] * length + list(range(position + 1, end))
        cigar = [(0, position + 1 - start), (1, length), (0, end - position - 1)]
    else:
        raise ValueError(f"unknown kind: {kind}")

    if kind == "ins":
        # Choose an inserted base differing from both neighbours so the
        # insertion does not land in a homopolymer run and get filtered.
        inserted = next(b for b in "ACGT" if b not in (reference[position], reference[position + 1]))
    substitutions = substitutions or {}
    sequence = "".join(inserted if p is None else substitutions.get(p, reference[p]) for p in bases)

    left, right = soft_clip
    if left:
        bases, sequence, cigar = ["clip"] * left + bases, clipped(left) + sequence, [(4, left)] + cigar
    if right:
        bases, sequence, cigar = bases + ["clip"] * right, sequence + clipped(right), cigar + [(4, right)]

    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = sequence
    read.flag = 16 if reverse else 0
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = mapq
    read.cigartuples = cigar
    # Must be assigned after query_sequence, which resets the qualities.
    qualities = qualities or {}
    read.query_qualities = array.array(
        "B", [qualities.get("inserted" if p is None else p, BASE_QUALITY) for p in bases])

    # PacBio consensus tags the callers read: per base, the subreads that
    # match the consensus (sm), that mismatch it (sx), and that cover it (sa,
    # run-length encoded as (length, coverage) pairs).
    subreads = subreads or {}
    per_base = [subreads.get("inserted" if p is None else p, (SUBREADS, 0, SUBREADS)) for p in bases]
    sm, sx, sa = (list(column) for column in zip(*per_base))
    if reverse:
        sm, sx, sa = sm[::-1], sx[::-1], sa[::-1]
    read.set_tag("np", 20, value_type="i")
    read.set_tag("ec", 20.0, value_type="f")
    read.set_tag("sm", array.array("B", sm))
    read.set_tag("sx", array.array("B", sx))
    read.set_tag("sa", array.array("I", run_length_encode(sa)))

    return read


def clipped(n):
    """n bases of sequence that belongs nowhere in the contig."""
    return ("TG" * n)[:n]


def run_length_encode(values):
    """[3, 3, 5] -> [2, 3, 1, 5]: (run length, value) pairs, as in the sa tag."""
    encoded = []
    for value in values:
        if encoded and encoded[-1] == value:
            encoded[-2] += 1
        else:
            encoded += [1, value]
    return encoded


def write_bam(path, reference, reads):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": len(reference)}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for read in sorted(reads, key=lambda read: read.reference_start):
            out.write(read)
    pysam.index(path)


def write_lines(path, lines):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as handle:
        for line in lines:
            handle.write(line + "\n")


def write_dataset(directory, reference, reads, mask=()):
    """
    Write a dataset for a test that needs its own reads: the assembly, the
    BAM, the mean coverage (the number of reads), a tissue label per read and
    a mask BED of the (start, end) intervals in `mask`. Returns the paths,
    keyed like conftest.fixture_paths.
    """
    paths = {
        "assembly": os.path.join(directory, "assembly.fa"),
        "bam": os.path.join(directory, "reads.bam"),
        "coverage": os.path.join(directory, "avg_coverage.txt"),
        "tissue": os.path.join(directory, "tissue.tsv"),
        "filters": os.path.join(directory, "mask.bed"),
    }
    write_fasta(paths["assembly"], CONTIG, reference)
    write_bam(paths["bam"], reference, reads)
    write_lines(paths["coverage"], [str(len(reads))])
    write_lines(paths["tissue"],
                [f"{read.query_name}\t{TISSUES[i % len(TISSUES)]}" for i, read in enumerate(reads)])
    write_lines(paths["filters"], [f"{CONTIG}\t{start}\t{end}" for start, end in mask])
    return paths


def write_tr_bed(path, tracts, model="M1"):
    """
    Write a repeat BED as workflow 2 does, header included. Each tract is
    (start, end, canonical_motif, k), optionally followed by its purity
    (default 100). Tracts may overlap.
    """
    lines = ["contig\tstart\tend\tcanonical_motif\tk\tcopies\ttrf_purity\tmodel_purity\tmodel"]
    for start, end, motif, k, *purity in tracts:
        purity = purity[0] if purity else 100.0
        lines.append(f"{CONTIG}\t{start}\t{end}\t{motif}\t{k}\t{(end - start) / k:g}"
                     f"\t{purity:.2f}\t{purity:.2f}\t{model}")
    write_lines(path, lines)


def build(fixtures):
    """Write the dataset under the directory `fixtures`."""
    reference = build_reference()

    # -- assemblies --------------------------------------------------------
    diploid = os.path.join(
        fixtures, "fasta", SPECIES, SAMPLE, "diploid_assembly", f"{SAMPLE}_diploid.fa"
    )
    haploid = os.path.join(
        fixtures, "fasta", SPECIES, SAMPLE, "haploid_assembly", f"{SAMPLE}_haploid.fa"
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
        fixtures, "bam", SPECIES, SAMPLE, "kinetics", f"{SAMPLE}_kinetics_diploid.bam"
    )
    write_bam(bam, reference, reads)

    # -- coverage ----------------------------------------------------------
    # Every read spans the whole contig, so depth equals the read count. The
    # callers keep positions within 0.5x-2x of this value.
    write_lines(
        os.path.join(fixtures, "coverage", f"{SAMPLE}_avg_coverage.txt"),
        [str(len(reads))],
    )

    # -- tissue classification --------------------------------------------
    write_lines(
        os.path.join(fixtures, "tissue", SPECIES, f"{SAMPLE}_sorted_tissue_classification.tsv"),
        [f"{read.query_name}\t{TISSUES[i % len(TISSUES)]}" for i, read in enumerate(reads)],
    )

    # -- masks -------------------------------------------------------------
    microsat_dir = os.path.join(fixtures, "microsat", SPECIES, SAMPLE)

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
    write_tr_bed(os.path.join(fixtures, "tr_bed", f"{CONTIG}_M1.bed"),
                 [(MICROSAT_START, MICROSAT_END, "AC", 2)])

    return reference, reads


def main():
    reference, reads = build(FIXTURES)
    print(f"Fixtures written to {FIXTURES}", file=sys.stderr)
    print(f"  contig {CONTIG}: {len(reference)} bp, {len(reads)} reads", file=sys.stderr)
    for name, kind, position, length, expectation in PLANTED:
        print(f"  {name:15s} {kind} {length}bp @ {position:5d}  {expectation}", file=sys.stderr)


if __name__ == "__main__":
    main()
