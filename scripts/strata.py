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

# What a read key sees: the read's name and its tissue label ("NA" without one).
Read = namedtuple("Read", "name tissue")

# What a site key sees, for the read bases tested: their base qualities, and
# the share of covering subreads that match each of them (None if the read has
# no subread tags).
Site = namedtuple("Site", "qualities support")


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
    read_keys = [
        ("tissue", lambda read: read.tissue),
    ]
    site_keys = [
        # The lowest base quality among the bases tested.
        ("kmerQ", lambda site: min(site.qualities)),
        # The lowest share of covering subreads that match, among the bases
        # tested; NA for reads without subread tags.
        ("subread_support", lambda site: "NA" if site.support is None else support_bin(min(site.support))),
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
