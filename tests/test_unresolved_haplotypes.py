"""
Unresolved-haplotype regions: a mismatch shared by enough reads at a high
enough base quality marks the span of those reads, which is then masked.
"""

import contextlib
import io

import make_fixtures as mf
import utilities

REFERENCE = mf.build_reference()


def regions(tmp_path, reads, min_base_quality=40, min_reads=2):
    """The regions resolve_haplotypes() reports for these reads plus ten clean ones."""
    clean = [mf.make_read(f"clean{i:02d}", REFERENCE) for i in range(10)]
    paths = mf.write_dataset(str(tmp_path), REFERENCE, reads + clean)
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        utilities.resolve_haplotypes(mf.CONTIG, 0, mf.CONTIG_LENGTH, paths["bam"], paths["assembly"],
                                     paths["coverage"], min_base_quality=min_base_quality, min_reads=min_reads)
    return [tuple(int(x) for x in line.split("\t")[1:]) for line in out.getvalue().splitlines()]


def other_base(position):
    return next(b for b in "ACGT" if b != REFERENCE[position])


def shared_mismatch(position, **qualities):
    """Two reads of different spans carrying the same mismatch."""
    change = {position: other_base(position)}
    return [mf.make_read("a", REFERENCE, start=1000, end=3000, substitutions=change, **qualities),
            mf.make_read("b", REFERENCE, start=1200, end=4000, substitutions=change, **qualities)]


def test_a_mismatch_shared_by_enough_reads_marks_their_span(tmp_path):
    assert regions(tmp_path, shared_mismatch(1500)) == [(1000, 4000)]


def test_fewer_reads_than_the_minimum_mark_nothing(tmp_path):
    assert regions(tmp_path, shared_mismatch(1500), min_reads=3) == []


def test_mismatches_below_the_base_quality_mark_nothing(tmp_path):
    assert regions(tmp_path, shared_mismatch(1500, qualities={1500: 30})) == []
    assert regions(tmp_path / "q30", shared_mismatch(1500, qualities={1500: 30}), min_base_quality=30) == \
        [(1000, 4000)]
