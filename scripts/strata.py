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
from collections import namedtuple

# What a read key sees: the read's name and its tissue label ("NA" without one).
Read = namedtuple("Read", "name tissue")

# What a site key sees: the base qualities of the read bases tested.
Site = namedtuple("Site", "qualities")


def keys(**parameters):
    """
    (read keys, site keys): lists of (column name, function). `parameters` are
    the config values the keys depend on.
    """
    read_keys = [
        ("tissue", lambda read: read.tissue),
    ]
    site_keys = [
        # The lowest base quality among the bases tested.
        ("kmerQ", lambda site: min(site.qualities)),
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
