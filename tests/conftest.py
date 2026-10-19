"""
Shared pytest setup.

Puts the repository root, scripts/ and workflows/ on sys.path so tests can
import the analysis modules and the workflows' helper, which use flat imports
(``from utilities import nonrepetitive_scan``) because they are normally run as
standalone scripts, and gwf puts a workflow's own folder on the path.
"""

import os
import sys

import pytest

import make_fixtures

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS_DIR = os.path.join(REPO_ROOT, "scripts")
WORKFLOW_DIR = os.path.join(REPO_ROOT, "workflows")

for path in (REPO_ROOT, SCRIPTS_DIR, WORKFLOW_DIR):
    if path not in sys.path:
        sys.path.insert(0, path)


@pytest.fixture(scope="session")
def fixtures(tmp_path_factory):
    """
    The miniature dataset of make_fixtures.py, built afresh for the session so
    a change to the generator can never leave stale fixtures behind.
    """
    directory = str(tmp_path_factory.mktemp("fixtures"))
    make_fixtures.build(directory)
    return directory


@pytest.fixture(scope="session")
def fixture_paths(fixtures):
    """The individual fixture files, keyed by role."""
    return {
        "assembly": os.path.join(
            fixtures, "fasta", "testspecies", "testsample",
            "diploid_assembly", "testsample_diploid.fa",
        ),
        "bam": os.path.join(
            fixtures, "bam", "testspecies", "testsample",
            "kinetics", "testsample_kinetics_diploid.bam",
        ),
        "coverage": os.path.join(fixtures, "coverage", "testsample_avg_coverage.txt"),
        "tissue": os.path.join(
            fixtures, "tissue", "testspecies",
            "testsample_sorted_tissue_classification.tsv",
        ),
        "filters": os.path.join(
            fixtures, "microsat", "testspecies", "testsample",
            "testcontig1_filters.tsv",
        ),
        "tr_bed": os.path.join(fixtures, "tr_bed", "testcontig1_M1.bed"),
    }
