"""
Unit tests for the pure helper functions in the analysis code.

One of them behaves differently from what its name promises; that test is
marked BEHAVIOUR LOCK and asserts what the code does today rather than what it
arguably should do, so that changing it shows up as a deliberate test edit
instead of a silent shift in published mutation rates.
"""

import pytest

from parse_trf import canonicalise, hamming_pct_matches
from read_tags import decode_runs
from utilities import (
    homopolymer_check,
    reverse_complement,
)


# -- reverse_complement -----------------------------------------------------


@pytest.mark.parametrize(
    "sequence,expected",
    [("ATGC", "GCAT"), ("A", "T"), ("AAAA", "TTTT"), ("GT", "AC"), ("", "")],
)
def test_reverse_complement(sequence, expected):
    assert reverse_complement(sequence) == expected


# -- homopolymer_check ------------------------------------------------------


@pytest.mark.parametrize("sequence", ["AAA", "TTTT", "A", "GG"])
def test_homopolymer_check_true_for_actual_homopolymers(sequence):
    assert homopolymer_check(sequence) is True


@pytest.mark.parametrize("sequence", ["ACG", "ACGT", "AACGG"])
def test_homopolymer_check_false_when_two_or_more_breaks(sequence):
    assert homopolymer_check(sequence) is False


@pytest.mark.parametrize("sequence", ["AC", "AAC", "ACC", "AAACCC"])
def test_homopolymer_check_treats_single_break_sequences_as_homopolymers(sequence):
    """
    BEHAVIOUR LOCK.

    The function returns True (-> the call is discarded) for anything with at
    most one change in base identity, so a plain 2-mer like "AC" is filtered as
    if it were a homopolymer run. Far broader than the name suggests.
    """
    assert homopolymer_check(sequence) is True


# -- decode_runs: the sa tag ------------------------------------------------


def test_decode_runs_gives_each_base_the_coverage_of_its_run():
    # (length, coverage) pairs: the first run covers bases 0-9, the second 10-29.
    coverage = decode_runs([10, 3, 20, 7])
    assert (len(coverage), coverage[5], coverage[15]) == (30, 3, 7)


# -- parse_trf helpers ------------------------------------------------------


@pytest.mark.parametrize(
    "motif,expected",
    [("GT", "AC"), ("TG", "AC"), ("CA", "AC"), ("AC", "AC"), ("TTA", "AAT")],
)
def test_canonicalise_returns_smallest_rotation_of_either_strand(motif, expected):
    assert canonicalise(motif) == expected


@pytest.mark.parametrize(
    "motif,reverse_complement_motif",
    [("AG", "CT"), ("AC", "GT"), ("AAC", "GTT")],
)
def test_canonicalise_collapses_reverse_complements(motif, reverse_complement_motif):
    """The same repeat read from the other strand is one motif, not two."""
    assert canonicalise(motif) == canonicalise(reverse_complement_motif)


@pytest.mark.parametrize("motif", ["AT", "CG", "ATAT"])
def test_canonicalise_leaves_self_complementary_motifs_alone(motif):
    assert canonicalise(motif) == motif


def test_hamming_pct_matches_perfect_repeat():
    assert hamming_pct_matches("ATATAT", "AT") == 100.0


def test_hamming_pct_matches_with_a_mismatch():
    assert hamming_pct_matches("ATATAG", "AT") == pytest.approx(5 / 6 * 100)


def test_hamming_pct_matches_empty_sequence():
    assert hamming_pct_matches("", "AT") == 100.0
