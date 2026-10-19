"""
Unit tests for the pure helper functions in the analysis code.

One of them behaves differently from what its name promises; that test is
marked BEHAVIOUR LOCK and asserts what the code does today rather than what it
arguably should do, so that changing it shows up as a deliberate test edit
instead of a silent shift in published mutation rates.
"""

import pytest

from parse_trf import canonicalise, hamming_pct_matches
from filtering_haploid_mappping import (
    find_indel_position_in_read,
    find_next_op_with_distance,
    find_previous_op_with_distance,
)
from utilities import (
    homopolymer_check,
    reverse_complement,
    subread_array,
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


# -- subread_array ----------------------------------------------------------


def test_subread_array_returns_coverage_for_the_containing_block():
    # (length, coverage) pairs: block 0 covers [0,10), block 1 covers [10,30).
    assert subread_array([10, 3, 20, 7], 5) == 3
    assert subread_array([10, 3, 20, 7], 15) == 7


def test_subread_array_returns_none_past_the_end():
    assert subread_array([10, 3], 50) is None


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


# -- CIGAR walking ----------------------------------------------------------


def test_find_previous_op_excludes_deletions_from_the_distance():
    cigar = [("10", "M"), ("2", "D"), ("5", "M")]
    # Walking back from the last block: 5 (M) then the D contributes nothing.
    assert find_previous_op_with_distance(cigar, 2, "D") == 5


def test_find_next_op_returns_none_when_absent():
    assert find_next_op_with_distance([("10", "M")], 0, "I") is None


def test_find_indel_position_in_read_matching_deletion():
    found, distance, softclip, others = find_indel_position_in_read(
        "10M2D10M", 10, "deletion", "AB"
    )
    assert found is True
    assert distance == 0
    assert softclip is False
    assert others is None


def test_find_indel_position_in_read_matching_insertion():
    found, _, _, _ = find_indel_position_in_read("10M3I10M", 10, "insertion", "ABC")
    assert found is True


def test_find_indel_position_in_read_rejects_a_size_mismatch():
    """An InDel is only concordant if type AND size agree."""
    found, _, _, _ = find_indel_position_in_read("10M2D10M", 10, "deletion", "A")
    assert found is False


def test_find_indel_position_in_read_flags_softclipping():
    _, _, softclip, _ = find_indel_position_in_read(
        "5S10M2D10M", 15, "deletion", "AB"
    )
    assert softclip is True


def test_find_indel_position_in_read_reports_no_match_for_a_clean_cigar():
    found, distance, _, _ = find_indel_position_in_read(
        "100M", 50, "deletion", "A"
    )
    assert found is False
    assert distance is None
