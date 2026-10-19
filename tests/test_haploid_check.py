"""
The haploid check, scripts/filtering_haploid_mappping.py: is a call also in its
read's alignment to the haploid assembly?

Findings H1-H4 of documents/nonrepetitive_indel_review.md each have tests here,
built on the review's CIGAR table: a 40 bp read whose InDel follows read base
17, against haploid CIGARs that carry exactly that InDel, so the right answer
is "found, distance 0" unless the test says otherwise.
"""

import pytest

from filtering_haploid_mappping import check_call

# Nothing in it can slide at the calls below: the bases on either side of
# each InDel differ from the bases that would have to repeat.
READ = "GATCCAGTGCATGACTCAGTCGATCAAGCTGGTACACGTA"


def check(cigar, indel_type, indel_string, flag=0, diploid_direction="+", start=17, read=READ):
    """haploid_check and distance_to_indel_in_hap of one call."""
    return check_call(read, flag, cigar, start, indel_type, indel_string, diploid_direction)[:2]


# -- the columns --------------------------------------------------------------


def test_a_matching_deletion_passes_at_distance_0():
    assert check_call(READ, 0, "17M1D23M", 17, "deletion", "G", "+") == (True, 0, False, None)


def test_a_matching_insertion_passes_at_distance_0():
    assert check("17M2I21M", "insertion", "AG") == (True, 0)


@pytest.mark.parametrize("cigar,indel_type,indel_string", [
    ("17M2D23M", "deletion", "G"),     # 2 bp where the call has 1
    ("17M1I22M", "insertion", "AG"),   # 1 bp where the call has 2
    ("17M2I21M", "deletion", "TG"),    # an insertion where the call is a deletion
    ("40M", "deletion", "G"),          # nothing
])
def test_an_indel_of_another_size_or_type_or_none_fails(cigar, indel_type, indel_string):
    """A call is only concordant if type and size agree."""
    assert check(cigar, indel_type, indel_string) == (False, None)


def test_soft_clipping_and_other_operations_are_flagged():
    assert check_call(READ, 0, "5S12M1D23M", 17, "deletion", "G", "+") == (True, 0, True, None)
    assert check_call(READ, 0, "2H15M1D23M", 17, "deletion", "G", "+") == (True, 0, False, "HPN")


# -- H1: inserted bases count as read bases ----------------------------------


@pytest.mark.parametrize("cigar", ["5M2I10M1D23M", "5M2I10M1D10M1D13M"])
def test_an_insertion_before_the_call_does_not_put_the_walk_behind(cigar):
    """H1: before the fix the first failed and the second matched 8 bp away."""
    assert check(cigar, "deletion", "G") == (True, 0)


# -- H2: on a flipped read the event is mapped, not the anchor ---------------


@pytest.mark.parametrize("cigar,indel_type,indel_string", [
    ("23M1D17M", "deletion", "G"),
    ("10M1I12M1D17M", "deletion", "G"),
    ("10M2I11M1D17M", "deletion", "G"),
    ("22M1I17M", "insertion", "A"),
    ("21M2I17M", "insertion", "AG"),
])
def test_a_flipped_read_maps_the_event_to_the_other_strand(cigar, indel_type, indel_string):
    """
    H2: the read is forward on the diploid assembly and reverse on the haploid
    one. Mapping the anchor put deletions 1 bp off and lost insertions of 2 bp.
    """
    assert check(cigar, indel_type, indel_string, flag=16) == (True, 0)


def test_a_read_reverse_on_the_diploid_assembly_maps_the_same_way():
    """The deletion above, seen from a read aligned reverse to the diploid assembly instead."""
    assert check("17M1D23M", "deletion", "C", flag=0, diploid_direction="-", start=23) == (True, 0)


# -- H4: only the call's own event passes, and only it has a distance --------


def test_a_later_event_of_the_same_size_does_not_pass():
    """H4: the first deletion at or after the call was taken, 13 bp on."""
    assert check("30M1D10M", "deletion", "G") == (False, None)


def test_a_call_without_a_match_has_no_distance():
    """H4: a 2 bp deletion where the 1 bp call is fails, and was given distance 0."""
    assert check("17M2D23M", "deletion", "G") == (False, None)


# -- H3: an InDel matches anywhere in the stretch it can slide over ----------

# 40 bp reads with a CA repeat at read bases 12-21 (0-based), which lost or
# gained one copy. Both aligners put a gap at its leftmost place on their own
# assembly's strand, so on a flipped read they put it at opposite ends of the
# repeat: the diploid alignment after base 12, the haploid one at the far end.
FLANK = "GATCCAGTGCTG"
READ_DELETION = FLANK + "CA" * 5 + "TGGATCCGTAAGCTTGAC"    # the gap can be after base 12 to 22
READ_INSERTION = FLANK + "CA" * 6 + "TGGATCCGTAAGCTTG"     # the inserted CA can be anywhere in 12-23


@pytest.mark.parametrize("cigar,flag,expected", [
    ("12M2D28M", 0, (True, 0)),       # same strand, same place
    ("18M2D22M", 16, (True, 10)),     # flipped: the far end of the repeat, 10 bases on
    ("16M2D24M", 16, (False, None)),  # flipped, 2 bases past the repeat: another InDel
])
def test_a_deletion_in_a_repeat_matches_anywhere_it_can_slide(cigar, flag, expected):
    assert check(cigar, "deletion", "CA", flag=flag, start=12, read=READ_DELETION) == expected


def test_an_insertion_in_a_repeat_matches_anywhere_it_can_slide():
    assert check("16M2I22M", "insertion", "CA", flag=16, start=12, read=READ_INSERTION) == (True, 10)

