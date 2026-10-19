"""
PacBio's per-read and per-base tags, read the same way for calls and
opportunities.

The per-base subread tags are stored in the order the read was sequenced, and
pbmm2 leaves them so. The read's sequence and base qualities, though, are
stored along the alignment: reverse-complemented and reversed for a read
aligned to the reverse strand of the assembly. So on such a read the tags run
backwards along the alignment, and the tags of aligned base i sit at index
read_length - 1 - i (S1 in documents/nonrepetitive_indel_review.md). The
functions here return them in the alignment's order.
"""


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
