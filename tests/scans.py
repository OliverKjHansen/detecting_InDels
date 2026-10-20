"""
Run the scans on a dataset and parse what they write. Shared by the tests.

A dataset is a dict of paths with the keys conftest.fixture_paths uses:
assembly, bam, coverage, tissue and filters (the mask).
"""

import io

import utilities
from make_fixtures import CONTIG

# The pipeline's parameters, as config.yaml sets them; a test can pass others.
COVERAGE_BAND = (0.5, 2.0)
INDEL_SIZES = (1, 3)
FLANK_WIDTH = 25
REACH = 6   # the longest period the repeat models annotate
MIN_TRACT_FLANK = 25
KEY_PARAMETERS = {"subread_support_bins": [0.8, 0.9, 1.0], "min_kmerq": 93,
                  "read_length_bins": [4000, 8000, 12000, 16000, 20000],
                  "read_end_distance_bins": [25, 100, 300, 1000],
                  "mapq_bins": [1, 20, 60], "passes_bins": [5, 10, 15, 20, 30],
                  "accuracy_bins": [0.99, 0.999, 0.9999], "strand_passes_bins": [2, 5, 10],
                  "nearby_event_bins": [5, 25, 100],
                  "tract_copies_bins": [5, 10, 20, 50], "tract_purity_bins": [90, 95, 99, 100]}


def rows_of(text):
    """Parse a scan's TSV output into one dict per row."""
    lines = text.strip().split("\n")
    header = lines[0].lstrip("#").split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]


def summed(tables):
    """
    Opportunity tables, as the scans write them, joined as the workflows join
    windows: {the values of every other column: summed count}, without zeros.
    """
    totals = {}
    for text in tables:
        for row in rows_of(text):
            key = tuple(value for column, value in row.items() if column != "count")
            totals[key] = totals.get(key, 0) + int(row["count"])
    return {key: count for key, count in totals.items() if count}


def nonrepetitive_scan(paths, window_start=None, window_end=None,
                       coverage_band=COVERAGE_BAND, indel_sizes=INDEL_SIZES):
    """Run nonrepetitive_scan() into strings: (calls, opportunities, FASTQ)."""
    calls, opportunities, fastq = io.StringIO(), io.StringIO(), io.StringIO()
    utilities.nonrepetitive_scan(
        paths["filters"], CONTIG, paths["bam"], paths["assembly"],
        paths["coverage"], paths["tissue"], window_start, window_end,
        calls_out=calls, opportunities_out=opportunities, fastq_out=fastq,
        coverage_band=coverage_band, indel_sizes=indel_sizes, key_parameters=KEY_PARAMETERS,
        flank_width=FLANK_WIDTH, classification_file=paths.get("classification", "NA"),
    )
    return calls.getvalue(), opportunities.getvalue(), fastq.getvalue()


def tandem_repeat_scan(beds, paths, window_start=None, window_end=None, coverage_band=COVERAGE_BAND):
    """Run tandem_repeat_scan() into strings: (calls per model, opportunities per model, FASTQ)."""
    calls = {model: io.StringIO() for model in beds}
    opportunities = {model: io.StringIO() for model in beds}
    fastq = io.StringIO()
    utilities.tandem_repeat_scan(
        beds, CONTIG, paths["bam"], paths["assembly"],
        paths["coverage"], paths["tissue"], window_start, window_end,
        calls_out=calls, opportunities_out=opportunities, fastq_out=fastq,
        coverage_band=coverage_band, key_parameters=KEY_PARAMETERS, flank_width=FLANK_WIDTH, reach=REACH,
        min_flank=MIN_TRACT_FLANK,
        classification_file=paths.get("classification", "NA"),
    )
    return ({model: out.getvalue() for model, out in calls.items()},
            {model: out.getvalue() for model, out in opportunities.items()},
            fastq.getvalue())
