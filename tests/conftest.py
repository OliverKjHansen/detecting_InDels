"""
Shared pytest setup.

Puts the repository root and scripts/ on sys.path so tests can import both
`config` and the analysis modules, which use flat imports
(``from utilities import small_indels``) because they are normally run as
standalone scripts.
"""

import os
import subprocess
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS_DIR = os.path.join(REPO_ROOT, "scripts")
FIXTURES = os.path.join(REPO_ROOT, "tests", "fixtures")

for path in (REPO_ROOT, SCRIPTS_DIR):
    if path not in sys.path:
        sys.path.insert(0, path)


@pytest.fixture(scope="session")
def fixtures():
    """Path to tests/fixtures/, generating it first if it is missing."""
    diploid = os.path.join(
        FIXTURES, "fasta", "testspecies", "testsample",
        "diploid_assembly", "testsample_diploid.fa",
    )
    if not os.path.exists(diploid):
        subprocess.run(
            [sys.executable, os.path.join(REPO_ROOT, "tests", "make_fixtures.py")],
            check=True,
        )
    return FIXTURES


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
