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
from collections import namedtuple

import numpy as np

from strata import Read

# One of a read's InDels, as pysam's pileup reports it: at its anchor, the
# aligned base before it (`anchor` on the assembly, `read_pos` in the read),
# with its length, a deletion or an insertion.
InDel = namedtuple("InDel", "deletion anchor read_pos length")


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


def describe(alignment, tissue, methylation_bin):
    """The Read that read keys see."""
    return Read(alignment.query_name, tissue, methylation_bin, alignment.query_length, alignment.is_reverse,
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


def indels(alignment):
    """
    The read's InDels, in the order of the alignment. Like pysam's pileup,
    it leaves out an InDel that doesn't follow an aligned base, such as an
    insertion right after a deletion.
    """
    found, query, position, after_aligned = [], 0, alignment.reference_start, False
    for op, length in alignment.cigartuples:
        if op in (1, 2) and after_aligned:
            found.append(InDel(op == 2, position - 1, query - 1, length))
        if op in (0, 7, 8, 1, 4):            # bases of the read
            query += length
        if op in (0, 7, 8, 2, 3):            # bases of the assembly
            position += length
        after_aligned = op in (0, 7, 8)
    return found


class Blocks:
    """Where a read's aligned bases are: its stretches of aligned bases, between InDels, clips and skips."""

    def __init__(self, alignment):
        self.positions, self.queries, self.lengths = [], [], []
        query, position = 0, alignment.reference_start
        for op, length in alignment.cigartuples:
            if op in (0, 7, 8):
                if self.positions and self.positions[-1] + self.lengths[-1] == position \
                        and self.queries[-1] + self.lengths[-1] == query:
                    self.lengths[-1] += length   # = and X runs make one stretch
                else:
                    self.positions.append(position)
                    self.queries.append(query)
                    self.lengths.append(length)
            if op in (0, 7, 8, 1, 4):
                query += length
            if op in (0, 7, 8, 2, 3):
                position += length

    def span(self, start, end):
        """
        The read bases aligned in [start, end) of the assembly, and those
        inserted between them, as a slice (first, last + 1) of the read. If a
        deletion takes all of [start, end), the two read bases on each side of
        it instead, as a call of that deletion is tested on. None if the
        stretch is before or after the read.
        """
        i = bisect_right(self.positions, start) - 1
        if i < 0 or self.positions[i] + self.lengths[i] <= start:
            i += 1   # start is before the read or in a gap: the next stretch's first base
            if i == len(self.positions):
                return None
            first = self.queries[i]
        else:
            first = self.queries[i] + start - self.positions[i]
        j = bisect_left(self.positions, end) - 1
        if j < 0:
            return None
        if j < i:   # all in the gap between stretches j and i
            return max(0, self.queries[j] + self.lengths[j] - 2), self.queries[i] + 2
        last = self.queries[j] + min(self.lengths[j], end - self.positions[j])
        return first, last


def nearest_outside(found, own, start, end):
    """
    Distance from [start, end) of the assembly to the nearest of the read's
    events, `found` (sorted), leaving out `own`, the anchors of the InDels the
    tract counts: 0 for an event inside, None if there is no other.
    """
    own = list(own)

    def other(position):
        """Whether the event at `position` is another's; else use up one of the tract's own there."""
        if position in own:
            own.remove(position)
            return False
        return True

    distances = []
    i = bisect_left(found, start)
    for position in found[i:]:
        if other(position):
            distances.append(max(0, position - (end - 1)))
            break
    for position in reversed(found[:i]):
        if other(position):
            distances.append(start - position)
            break
    return min(distances) if distances else None
