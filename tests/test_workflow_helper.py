"""
The helper every workflow loads config.yaml with: samples, paths built from
the `layout` templates, and the contigs long enough to analyse. Run against
the fixture config of workflow_config.py.
"""

import os

import pytest

from workflow_config import write_config


@pytest.fixture
def pipeline(tmp_path, monkeypatch):
    """The helper, loaded the way a workflow loads it: from the directory gwf runs in."""
    root = os.path.realpath(tmp_path)
    write_config(root)
    monkeypatch.chdir(root)
    from pipeline import Pipeline
    return Pipeline()


@pytest.fixture
def root(pipeline):
    return os.getcwd()


def test_paths_come_from_the_layout_templates(pipeline, root):
    assert pipeline.path("assembly", "ape", "s1") == f"{root}/fasta/ape/s1/diploid_assembly/s1_diploid.fa"
    # `output` is relative in the config, so it is made absolute from the
    # directory gwf runs in.
    assert pipeline.path("coverage", "ape", "s1") == f"{root}/output/precomute/ape/s1/avg_coverage.txt"
    assert pipeline.path("contig_fasta", "monkey", "s3", contig="c2") == \
        f"{root}/output/precomute/monkey/s3/contigs/c2.fa"


def test_an_override_replaces_the_standard_alignment(pipeline, root):
    assert pipeline.alignment("ape", "s2") == f"{root}/override/s2_diploid.bam"
    assert pipeline.alignment("ape", "s1") == f"{root}/bam/ape/s1/kinetics/s1_kinetics_diploid.bam"


def test_long_contigs_are_strictly_longer_than_the_minimum_in_index_order(pipeline):
    # c3 is exactly min_contig_length (1 Mb), so it is left out.
    assert list(pipeline.long_contigs("ape", "s1").items()) == [("c1", 12_000_000), ("c2", 1_500_000)]
    assert list(pipeline.long_contigs("monkey", "s3")) == ["c2"]


def test_a_sample_is_listed_in_a_table_only_under_its_own_species(pipeline):
    assert pipeline.is_listed("tissue_samples", "ape", "s1")
    assert not pipeline.is_listed("tissue_samples", "ape", "s2")
    assert not pipeline.is_listed("tissue_samples", "monkey", "s3")   # no monkey entry
    assert not pipeline.is_listed("no_such_table", "ape", "s1")


def test_a_templates_folder_needs_no_per_file_fields(pipeline, root):
    assert pipeline.folder("contig_fasta", "ape", "s1") == f"{root}/output/precomute/ape/s1/contigs"
    assert pipeline.folder("unresolved_window", "ape", "s1") == \
        f"{root}/output/precomute/ape/s1/unresolved_haplotypes/chunks"


def test_relative_paths_can_resolve_against_another_folder(tmp_path, monkeypatch):
    root = os.path.realpath(tmp_path / "elsewhere")
    os.makedirs(root)
    write_config(root)
    monkeypatch.chdir(tmp_path)
    from pipeline import Pipeline
    pipeline = Pipeline(f"{root}/config.yaml", base=root)
    assert pipeline.path("coverage", "ape", "s1") == f"{root}/output/precomute/ape/s1/avg_coverage.txt"


def test_a_folder_takes_the_fields_its_own_path_needs(pipeline, root):
    assert pipeline.folder("repetitive_calls", "ape", "s1", model="M1") == \
        f"{root}/output/repetitive_indels/ape/s1/M1/calls"
    with pytest.raises(ValueError):
        pipeline.folder("repetitive_calls", "ape", "s1")   # the model is part of the folder


def test_the_classification_scripts_get_their_cell_types_and_bins_from_the_config(pipeline):
    assert pipeline.classification_arguments().split(" ") == [
        "Sperm.",
        "Blood.Monocytes,Bladder.SmoothMuscle,Breast.BasalEpithelial,"
        "Lung.AlveolarMacrophages,Colon.RightEndocrine,Thyroid.Epithelial",
        "-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10",
        "-12,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,12",
    ]
