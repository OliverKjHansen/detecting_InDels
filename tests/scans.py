"""
Run the scans on a dataset and parse what they write. Shared by the tests.

A dataset is a dict of paths with the keys conftest.fixture_paths uses:
assembly, bam, coverage, tissue and filters (the mask).
"""

import io

import utilities
from make_fixtures import CONTIG


def rows_of(text):
    """Parse a scan's TSV output into one dict per row."""
    lines = text.strip().split("\n")
    header = lines[0].lstrip("#").split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]


def nonrepetitive_scan(paths, window_start=None, window_end=None):
    """Run nonrepetitive_scan() into strings: (calls, opportunities, FASTQ)."""
    calls, opportunities, fastq = io.StringIO(), io.StringIO(), io.StringIO()
    utilities.nonrepetitive_scan(
        paths["filters"], CONTIG, paths["bam"], paths["assembly"],
        paths["coverage"], paths["tissue"], window_start, window_end,
        calls_out=calls, opportunities_out=opportunities, fastq_out=fastq,
    )
    return calls.getvalue(), opportunities.getvalue(), fastq.getvalue()


def tandem_repeat_scan(beds, paths, window_start=None, window_end=None):
    """Run tandem_repeat_scan() into strings: (calls per model, opportunities per model, FASTQ)."""
    calls = {model: io.StringIO() for model in beds}
    opportunities = {model: io.StringIO() for model in beds}
    fastq = io.StringIO()
    utilities.tandem_repeat_scan(
        beds, CONTIG, paths["bam"], paths["assembly"],
        paths["coverage"], paths["tissue"], window_start, window_end,
        calls_out=calls, opportunities_out=opportunities, fastq_out=fastq,
    )
    return ({model: out.getvalue() for model, out in calls.items()},
            {model: out.getvalue() for model, out in opportunities.items()},
            fastq.getvalue())
