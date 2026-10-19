"""
The keys that split calls and opportunities alike, each defined once.

A rate is calls over opportunities, and a filter on the calls needs the same
filter on the opportunities, or it lowers the rate whatever it removes. So
every feature a call might be filtered on is a key here. Both scans write each
key as a column of every call and count the opportunities per value of it,
from the same definition, so adding a key adds it to both tables.

There are two kinds of key:

  read keys   a property of the read, computed once per read from a Read;
  site keys   a property of the bases a call or an opportunity is tested
              on, computed from a Site.

A scan describes each read and each tested site with these two records, the
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
Site = namedtuple("Site", "qualities support read_end_distance nearest_event")


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
        ("smrt_cell", lambda read: read.name.split("/", 1)[0]),
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


class Strata:
    """The keys of one scan: their column names, and their values for a read or a site."""

    def __init__(self, **parameters):
        self.read_keys, self.site_keys = keys(**parameters)
        self.columns = [name for name, _ in self.read_keys + self.site_keys]

    def of_read(self, read):
        return tuple(key(read) for _, key in self.read_keys)

    def of_site(self, site):
        return tuple(key(site) for _, key in self.site_keys)
