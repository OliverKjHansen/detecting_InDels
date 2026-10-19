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

# What a read key sees: the read's name, its tissue label ("NA" without one),
# its length in bases, and whether it is aligned to the reverse strand.
Read = namedtuple("Read", "name tissue length reverse")

# What a site key sees, for the read bases tested: their base qualities, the
# share of covering subreads that match each of them (None if the read has no
# subread tags), and the distance from the anchor to the nearer end of the
# read's alignment.
Site = namedtuple("Site", "qualities support read_end_distance")


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
    min_kmerq = parameters["min_kmerq"]
    read_keys = [
        ("tissue", lambda read: read.tissue),
        ("read_length_bin", lambda read: length_bin(read.length)),
        # The strand of the assembly the read is aligned to.
        ("aligned_strand", lambda read: "-" if read.reverse else "+"),
        # The SMRT Cell: the movie, the first part of a PacBio read name.
        ("smrt_cell", lambda read: read.name.split("/", 1)[0]),
    ]
    site_keys = [
        # Whether every base tested reaches min_kmerq: pass/fail rather than
        # the 94 quality levels, to keep the opportunity table small.
        ("kmerQ_pass", lambda site: "pass" if min(site.qualities) >= min_kmerq else "fail"),
        # The lowest share of covering subreads that match, among the bases
        # tested; NA for reads without subread tags.
        ("subread_support", lambda site: "NA" if site.support is None else support_bin(min(site.support))),
        ("read_end_distance_bin", lambda site: distance_bin(site.read_end_distance)),
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
