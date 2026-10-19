"""
What a read says about itself, read the same way for calls and opportunities:
PacBio's per-read and per-base tags, and the read's own InDels and mismatches
against the assembly.

The per-base subread tags are stored in the order the read was sequenced, and
pbmm2 leaves them so. The read's sequence and base qualities, though, are
stored along the alignment: reverse-complemented and reversed for a read
aligned to the reverse strand of the assembly. So on such a read the tags run
backwards along the alignment, and the tags of aligned base i sit at index
read_length - 1 - i (S1 in documents/nonrepetitive_indel_review.md). The
functions here return them in the alignment's order.
"""
from bisect import bisect_left, bisect_right

import numpy as np

from strata import Read


def decode_runs(runs):
    """sa's run-length encoding, (length, value) pairs, expanded to one value per base."""
    values = []
    for i in range(0, len(runs), 2):
        values.extend([runs[i + 1]] * runs[i])
    return values


def subread_tags(alignment):
    """
    (sm, sx, sa) per base, in the alignment's order: the subreads that match
    the consensus, that mismatch it, and that cover it. None if the read has
    no subread tags, as Sequel II reads don't.
    """
    if not alignment.has_tag("sa"):
        return None
    sm, sx = list(alignment.get_tag("sm")), list(alignment.get_tag("sx"))
    sa = decode_runs(alignment.get_tag("sa"))
    if alignment.is_reverse:
        sm, sx, sa = sm[::-1], sx[::-1], sa[::-1]
    return sm, sx, sa


def subread_support(tags):
    """Per base, the share of covering subreads that match (sm / sa); None without tags."""
    if tags is None:
        return None
    sm, _, sa = tags
    return [matching / covering if covering else 0.0 for matching, covering in zip(sm, sa)]


def tag(alignment, name):
    """A per-read tag's value, or None if the read doesn't have it."""
    return alignment.get_tag(name) if alignment.has_tag(name) else None


def describe(alignment, tissue):
    """The Read that read keys see."""
    return Read(alignment.query_name, tissue, alignment.query_length, alignment.is_reverse,
                alignment.mapping_quality, tag(alignment, "np"), tag(alignment, "rq"),
                tag(alignment, "fn"), tag(alignment, "rn"))


def events(alignment, reference):
    """
    The assembly positions of the read's own InDels and mismatches, sorted. An
    InDel is at its anchor, the base before the gap, as a call is. `reference`
    is the assembly from the read's first aligned base to its last.
    """
    found = []
    sequence, query, position = alignment.query_sequence, 0, alignment.reference_start
    for op, length in alignment.cigartuples:
        if op in (0, 7, 8):                  # aligned bases
            if op == 8:
                found.extend(range(position, position + length))
            elif op == 0:
                offset = position - alignment.reference_start
                read_bases = np.frombuffer(sequence[query:query + length].encode(), dtype=np.uint8)
                assembly_bases = np.frombuffer(reference[offset:offset + length].encode(), dtype=np.uint8)
                found.extend((np.nonzero(read_bases != assembly_bases)[0] + position).tolist())
            query, position = query + length, position + length
        elif op == 1:                        # insertion
            found.append(position - 1)
            query += length
        elif op == 2:                        # deletion
            found.append(position - 1)
            position += length
        elif op == 4:                        # soft clip
            query += length
        elif op == 3:                        # skipped assembly
            position += length
    return sorted(found)


def nearest_other_event(found, anchor):
    """Distance from the anchor to the nearest event not at the anchor itself, or None."""
    before, after = bisect_left(found, anchor), bisect_right(found, anchor)
    distances = []
    if before > 0:
        distances.append(anchor - found[before - 1])
    if after < len(found):
        distances.append(found[after] - anchor)
    return min(distances) if distances else None
