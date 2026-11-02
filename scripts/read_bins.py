"""
Every read's germline/somatic bin, compact (ticket 44). The scans of a
methylation sample look up the bins of their window's reads; reading the
sample's whole read classification table in every scan job made the jobs wait
on the storage (ticket 42). So workflow 1 writes the table once more, as a key
per read name, sorted, and the read's bin as a one-byte code: 9 bytes a read,
against the table's ~60.

    python read_bins.py <read_classification.tsv> <read_bins.npz>

A read's key is 64 bits of the BLAKE2b hash of its name, the same on every
machine and run, unlike Python's hash(). Two of a sample's ~10 M reads share a
key with a chance of about 3 in a million: the table is then refused, rather
than give one of them the other's bin. A read name listed twice keeps its last
row's bin, as the scans' reading of the text table did.
"""
import argparse
import contextlib
import hashlib
import os

import numpy


def key(name):
    """A read name's 64-bit key."""
    return int.from_bytes(hashlib.blake2b(name.encode(), digest_size=8).digest(), "little")


class ReadBins:
    """A sample's read bins, looked up by read name with get(), as in a dict."""

    def __init__(self, path):
        with numpy.load(path) as table:
            self.keys, self.codes = table["keys"], table["codes"]
            self.labels = [str(label) for label in table["labels"]]

    def get(self, name, default=None):
        wanted = numpy.uint64(key(name))
        i = numpy.searchsorted(self.keys, wanted)
        if i < len(self.keys) and self.keys[i] == wanted:
            return self.labels[self.codes[i]]
        return default


def load(path):
    """The bins of a sample's reads; for a sample without them ("NA"), none."""
    return {} if path == "NA" else ReadBins(path)


def write(table, out):
    """
    The compact form of read_classification.sh's table (the file `table`),
    written to `out` whole or not at all.
    """
    with open(table) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        read, label = header.index("readID"), header.index("classification_bin")
        label_of = {fields[read]: fields[label] for fields in (line.rstrip("\n").split("\t") for line in handle)}
    labels = sorted(set(label_of.values()))
    if len(labels) > 256:
        raise ValueError(f"{table}: {len(labels)} bins, more than a byte's codes")
    code_of = {label: code for code, label in enumerate(labels)}
    keys = numpy.fromiter((key(name) for name in label_of), dtype=numpy.uint64, count=len(label_of))
    codes = numpy.fromiter((code_of[label] for label in label_of.values()), dtype=numpy.uint8, count=len(label_of))
    order = numpy.argsort(keys)
    keys, codes = keys[order], codes[order]
    if (keys[1:] == keys[:-1]).any():
        raise ValueError(f"{table}: two read names share a key")
    try:
        with open(f"{out}.tmp", "wb") as handle:   # a file, so numpy doesn't add .npz to the name
            numpy.savez(handle, keys=keys, codes=codes, labels=numpy.array(labels, dtype=str))
    except BaseException:
        with contextlib.suppress(FileNotFoundError):
            os.remove(f"{out}.tmp")
        raise
    os.replace(f"{out}.tmp", out)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("table", help="read_classification.sh's table of the sample's reads")
    parser.add_argument("out", help="the compact table, .npz")
    args = parser.parse_args()
    write(args.table, args.out)


if __name__ == "__main__":
    main()
