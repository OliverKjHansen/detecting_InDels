"""
The keys that split calls and opportunities alike, each defined once.

A rate is calls over opportunities, and a filter on the calls needs the same
filter on the opportunities, or it lowers the rate whatever it removes. So
every feature a call might be filtered on is a key here. Both scans write each
key as a column of every call and count the opportunities per value of it,
from the same definition, so adding a key adds it to both tables.

There are three kinds of key:

  read keys   a property of the read, computed once per read from a Read;
  site keys   a property of the bases a call or an opportunity is tested
              on, computed from a Site;
  tract keys  workflow 4 only: a property of a tract in one read that spans
              it, computed from a Tract.

A scan describes each read and each tested site with these records, the
same way for a call as for an opportunity. A key that needs more about a read
or a site adds a field here, and both scans have to fill it.
"""
import bisect
from collections import namedtuple

# What a read key sees: the read's name, its tissue label and its
# germline/somatic bin from the methylation model ("NA" without one), its
# length in bases, whether it is aligned to the reverse strand, its mapping
# quality to the diploid assembly, and its PacBio tags np (passes), rq
# (predicted accuracy), fn and rn (full passes on each strand of the
# molecule), None where the read has no such tag.
Read = namedtuple("Read", "name tissue methylation_bin length reverse mapq passes accuracy "
                          "forward_passes reverse_passes")

# What a site key sees, for the read bases tested: their base qualities, the
# share of covering subreads that match each of them (None if the read has no
# subread tags), the distance from the anchor to the nearer end of the read's
# alignment, and the distance from the anchor to the read's nearest other
# InDel or mismatch (None if it has none).
#
# Workflow 3 also fills (ticket 49), None elsewhere: the disagreeing subreads
# (sx) at the two read bases next to the event, L1 and R1 (None without subread
# tags); for a deletion, whether the deleted bases repeat right before or right
# after them (None for an insertion); and whether the base after the event is
# masked (only a deletion's can be: the mask test covers an insertion's).
Site = namedtuple("Site", "qualities support read_end_distance nearest_event "
                          "junction_sx repeats_next base_after_masked", defaults=(None, None, None))

# What a tract key sees, for one tract in one read that spans it: the tract's
# copies and purity (model_purity, %), and how many bases of the read's
# alignment lie before the tract and after it.
Tract = namedtuple("Tract", "copies purity left_flank right_flank")


def smrt_cell(read_name):
    """The SMRT Cell (movie) a read comes from: the first part of a PacBio read name."""
    return read_name.split("/", 1)[0]


def binned_or_missing(edges, missing="NA"):
    """As binned(), with `missing` for a value of None."""
    bin_of = binned(edges)
    return lambda value: missing if value is None else bin_of(value)


def binned(edges):
    """A function giving a value's bin among the edges: "<a", "a-b", ..., ">=z"."""
    labels = ([f"<{edges[0]:g}"] + [f"{low:g}-{high:g}" for low, high in zip(edges, edges[1:])]
              + [f">={edges[-1]:g}"])
    return lambda value: labels[bisect.bisect_right(edges, value)]


def keys(**parameters):
    """
    (read keys, site keys): lists of (column name, function). `parameters` are
    the `strata` block of config.yaml.
    """
    support_bin = binned(parameters["subread_support_bins"])
    length_bin = binned(parameters["read_length_bins"])
    distance_bin = binned(parameters["read_end_distance_bins"])
    mapq_bin = binned(parameters["mapq_bins"])
    passes_bin = binned_or_missing(parameters["passes_bins"])
    accuracy_bin = binned_or_missing(parameters["accuracy_bins"])
    strand_passes_bin = binned_or_missing(parameters["strand_passes_bins"])
    event_bin = binned_or_missing(parameters["nearby_event_bins"], missing="none")
    min_kmerq = parameters["min_kmerq"]
    read_keys = [
        ("tissue", lambda read: read.tissue),
        # All 14 bins of classify_reads.awk, so the germline threshold
        # can be chosen later; NA for reads or samples without a score.
        ("methylation_bin", lambda read: read.methylation_bin),
        ("read_length_bin", lambda read: length_bin(read.length)),
        # The strand of the assembly the read is aligned to.
        ("aligned_strand", lambda read: "-" if read.reverse else "+"),
        # The SMRT Cell: the movie, the first part of a PacBio read name.
        ("smrt_cell", lambda read: smrt_cell(read.name)),
        # The read's mapping quality to the diploid assembly (N1).
        ("mapQ_bin", lambda read: mapq_bin(read.mapq)),
        ("np_bin", lambda read: passes_bin(read.passes)),
        ("rq_bin", lambda read: accuracy_bin(read.accuracy)),
        # Full passes on the weaker strand of the molecule.
        ("strand_passes_bin", lambda read: strand_passes_bin(
            None if read.forward_passes is None or read.reverse_passes is None
            else min(read.forward_passes, read.reverse_passes))),
    ]
    site_keys = [
        # Whether every base tested reaches min_kmerq: pass/fail rather than
        # the 94 quality levels, to keep the opportunity table small.
        ("kmerQ_pass", lambda site: "pass" if min(site.qualities) >= min_kmerq else "fail"),
        # The lowest share of covering subreads that match, among the bases
        # tested; NA for reads without subread tags.
        ("subread_support", lambda site: "NA" if site.support is None else support_bin(min(site.support))),
        ("read_end_distance_bin", lambda site: distance_bin(site.read_end_distance)),
        # The read's nearest other InDel or mismatch; "none" if it has none.
        ("nearby_event_bin", lambda site: event_bin(site.nearest_event)),
    ]
    return read_keys, site_keys


def yes_no(value):
    return "NA" if value is None else ("yes" if value else "no")


def nonrepetitive_site_keys(**parameters):
    """
    Workflow 3's own site keys (ticket 49), on top of keys(): coarse, so the
    opportunity cube stays small; the calls carry the details as columns of
    their own (utilities.CALL_DETAIL_COLUMNS).
    """
    flank_bin = binned(parameters["flank_quality_bins"])
    return [
        # kmerQ by position (artefact_descriptions/kmerq_results.md): the
        # lowest base quality of the two bases next to the event, L1 and R1,
        # and of the two outside them, L2 and R2. site.qualities is L2 L1 R1 R2.
        ("inner_flank_q", lambda site: flank_bin(min(site.qualities[1], site.qualities[2]))),
        ("outer_flank_q", lambda site: flank_bin(min(site.qualities[0], site.qualities[3]))),
        # A subread that disagrees with the read at L1 or R1
        # (artefact_descriptions/false_positive_deletions.md); NA without tags.
        ("junction_disagreement", lambda site: yes_no(None if site.junction_sx is None else max(site.junction_sx) > 0)),
        ("repeats_next", lambda site: yes_no(site.repeats_next)),
        ("base_after_masked", lambda site: yes_no(site.base_after_masked)),
    ]


def tract_keys(**parameters):
    """Workflow 4's tract keys: a list of (column name, function). `parameters` as for keys()."""
    copies_bin = binned(parameters["tract_copies_bins"])
    purity_bin = binned(parameters["tract_purity_bins"])
    # A flank is the read's distance to its end, seen from the tract.
    flank_bin = binned(parameters["read_end_distance_bins"])
    return [
        ("tract_copies_bin", lambda tract: copies_bin(tract.copies)),
        ("tract_purity_bin", lambda tract: purity_bin(tract.purity)),
        ("left_flank_bin", lambda tract: flank_bin(tract.left_flank)),
        ("right_flank_bin", lambda tract: flank_bin(tract.right_flank)),
    ]


class Strata:
    """
    The keys of one scan: their column names, and their values for a read, a
    site or a tract. `nonrepetitive` gives workflow 3's (ticket 49): its own
    site keys added, the keys under `nonrepetitive_without` left out, and the
    bins under `nonrepetitive_bins` in place of the shared ones.
    """

    def __init__(self, nonrepetitive=False, **parameters):
        if nonrepetitive:
            parameters = {**parameters, **(parameters.get("nonrepetitive_bins") or {})}
        self.read_keys, self.site_keys = keys(**parameters)
        if nonrepetitive:
            without = set(parameters.get("nonrepetitive_without") or ())
            self.read_keys = [key for key in self.read_keys if key[0] not in without]
            self.site_keys = [key for key in self.site_keys if key[0] not in without] + nonrepetitive_site_keys(**parameters)
        self.columns = [name for name, _ in self.read_keys + self.site_keys]
        self.tract_keys = tract_keys(**parameters)
        self.tract_columns = [name for name, _ in self.tract_keys]

    def of_read(self, read):
        return tuple(key(read) for _, key in self.read_keys)

    def of_site(self, site):
        return tuple(key(site) for _, key in self.site_keys)

    def of_tract(self, tract):
        return tuple(key(tract) for _, key in self.tract_keys)
