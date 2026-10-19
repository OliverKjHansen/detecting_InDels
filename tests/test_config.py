"""
Checks on config.yaml and the workflow sources.

These need no cluster storage and no gwf: they read the YAML and parse the
workflow files as text. The full dependency graph is validated on the cluster
with `gwf -f workflows/<name>.py status`, where the assemblies actually exist.
"""

import ast
import os

import pytest
import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKFLOW_DIR = os.path.join(REPO_ROOT, "workflows")
WORKFLOW_NAMES = ["precompute", "repetitive_dna", "indel_mutation_call", "repetitive_indel_call"]

REQUIRED_ROOTS = ["fasta", "bam", "tissue", "microsat", "output_root"]


@pytest.fixture(scope="module")
def cfg():
    with open(os.path.join(REPO_ROOT, "config.yaml")) as handle:
        return yaml.safe_load(handle)


def strip_comments(source):
    """Drop whole-line `#` comments, which legitimately discuss the old paths."""
    return "\n".join(
        line for line in source.splitlines() if not line.lstrip().startswith("#")
    )


@pytest.fixture(scope="module")
def workflow_sources():
    sources = {}
    for name in WORKFLOW_NAMES:
        with open(os.path.join(WORKFLOW_DIR, f"{name}.py")) as handle:
            sources[name] = handle.read()
    return sources


# -- config.yaml shape ------------------------------------------------------


def test_required_top_level_keys_are_present(cfg):
    for key in ("samples", "tandem_repeat_samples", "min_contig_length",
                "account", "paths", "tools", "layout", "remote"):
        assert key in cfg, f"config.yaml is missing '{key}'"


def test_every_path_root_the_workflows_need_is_defined(cfg):
    for root in REQUIRED_ROOTS:
        assert root in cfg["paths"], f"paths.{root} is missing"


def test_layout_templates_resolve_against_the_declared_roots(cfg):
    """
    Every placeholder in a layout template must be either a configured root or
    one of the per-target fields the workflows supply.
    """
    roots = {key: "ROOT" for key in cfg["paths"]}
    fields = {"species": "SPECIES", "sample": "SAMPLE", "contig": "CONTIG", "model": "MODEL",
              "assembly_type": cfg["assembly_type"]}

    for key, template in cfg["layout"].items():
        try:
            template.format(**roots, **fields)
        except KeyError as exc:
            pytest.fail(f"layout.{key} references undefined placeholder {exc}")


def test_output_templates_are_sample_keyed(cfg):
    """Everything the pipeline writes lives under {species}/{sample}."""
    for key, template in cfg["layout"].items():
        if template.startswith("{output_root}/"):
            assert "{species}/{sample}" in template, (
                f"layout.{key} is not keyed by species and sample"
            )


def test_tandem_repeat_samples_are_a_subset_of_samples(cfg):
    """
    A sample can be InDel-called without repeat annotation (workflow C skips
    it), but annotating a sample that is not in `samples` is a typo.
    """
    for species, samples in cfg["tandem_repeat_samples"].items():
        assert species in cfg["samples"], f"{species} is not in samples"
        assert set(samples) <= set(cfg["samples"][species]), (
            f"tandem_repeat_samples has entries for {species} not in samples"
        )


def test_sample_lists_are_lists(cfg):
    """A bare string would silently iterate character by character."""
    for table in ("samples", "tandem_repeat_samples"):
        for species, samples in cfg[table].items():
            assert isinstance(samples, list), f"{table}.{species} must be a list"


def test_min_contig_length_is_a_positive_int(cfg):
    assert isinstance(cfg["min_contig_length"], int)
    assert cfg["min_contig_length"] > 0


def test_remote_block_is_complete(cfg):
    for key in ("host", "workdir", "local_dir"):
        assert cfg["remote"].get(key), f"remote.{key} is not set"


# -- workflow sources -------------------------------------------------------


@pytest.mark.parametrize("name", WORKFLOW_NAMES)
def test_workflow_parses_as_python(name, workflow_sources):
    ast.parse(workflow_sources[name])


@pytest.mark.parametrize("name", WORKFLOW_NAMES)
def test_workflow_does_not_reach_into_another_checkout(name, workflow_sources):
    """
    The point of the migration: recipes must run the scripts in this repository.
    The originals pointed at ~/mutationalscanning/.../endjoining/scripts.
    """
    code = strip_comments(workflow_sources[name])
    assert "endjoining" not in code
    assert "microsatellites/scripts" not in code


@pytest.mark.parametrize("name", WORKFLOW_NAMES)
def test_workflow_reads_config_yaml_directly(name, workflow_sources):
    source = workflow_sources[name]
    assert "yaml.safe_load" in source
    assert "from config import" not in source
    assert "from common import" not in source


@pytest.mark.parametrize("name", ["repetitive_dna"])
def test_workflow_is_rooted_at_the_repository(name, workflow_sources):
    """gwf otherwise uses workflows/ as working_dir, putting output/ inside it."""
    assert "working_dir=" in strip_comments(workflow_sources[name])


@pytest.mark.parametrize("name", WORKFLOW_NAMES)
def test_workflow_hardcodes_no_samples(name, workflow_sources):
    """Samples come from config.yaml, not a sample_dict literal in the workflow."""
    assert "sample_dict = {" not in strip_comments(workflow_sources[name])


def test_repetitive_dna_reads_the_contigs_precompute_writes(workflow_sources):
    assert "{outdir}/contigs/{contig}.fa" in workflow_sources["precompute"]
    assert 'path("contig_fasta"' in workflow_sources["repetitive_dna"]


@pytest.mark.parametrize("name", WORKFLOW_NAMES)
def test_workflow_hardcodes_no_cluster_data_paths(name, workflow_sources):
    """Data locations belong in config.yaml, not inline in the workflow."""
    source = workflow_sources[name]
    assert "/faststorage/" not in source
    assert "DerivedData" not in source


def test_filter_bed_merges_both_mask_sources(workflow_sources):
    """
    Regression guard for incomplete.md 2.1: `bedtools merge` takes a single -i,
    so `merge -i - -i {unresolved_bed}` silently dropped the microsatellites.
    Both sources must reach one stream before the merge.
    """
    source = workflow_sources["indel_mutation_call"]
    assert "bedtools merge -i - >" in source
    assert "bedtools merge -i - -i" not in source
    assert "cut -f1-3 {unresolved_bed} >> " in source


def test_quality_plots_target_runs_the_r_script(workflow_sources):
    """
    Regression guard for incomplete.md 2.2: the target declared the PDF as an
    output while the Rscript call was commented out, so it could never complete.
    """
    source = workflow_sources["indel_mutation_call"]
    assert "Rscript {SCRIPTS}/plot_indels.R {merged_indels}" in source


# -- output layout conventions ----------------------------------------------
#
# The tree is sample-first and split in two: output/results/{species}/{sample}/
# holds the small tables that get fetched and plotted, output/work/{species}/
# {sample}/ holds everything regenerable. The directory carries species, sample,
# pipeline and model; filenames never repeat them.


@pytest.mark.parametrize("name", WORKFLOW_NAMES)
def test_output_roots_are_sample_keyed(name, workflow_sources):
    """
    Every output path descends from work/ or results/ keyed {species}/{sample}.
    The old layout keyed final/ and plots/ by {sample} alone, so two species
    sharing a sample name would have overwritten each other.
    """
    code = strip_comments(workflow_sources[name])
    for tree in ("work", "results"):
        for line in code.splitlines():
            if f'{{OUTPUT}}/{tree}/' in line:
                assert "{species}/{sample}" in line, (
                    f"{name}: {tree}/ path is not keyed by species and sample: {line.strip()}"
                )


@pytest.mark.parametrize("name", WORKFLOW_NAMES)
def test_no_legacy_output_directories(name, workflow_sources):
    """The stage-first tree and its tandem_repeat_ mirror are both gone."""
    code = strip_comments(workflow_sources[name])
    for legacy in ('{OUTPUT}/final', '{OUTPUT}/plots', '{OUTPUT}/small_indels',
                   '{OUTPUT}/opportunities', '{OUTPUT}/mutation_fastq',
                   '{OUTPUT}/haploid_mapping', '{OUTPUT}/tandem_repeat_'):
        assert legacy not in code, f"{name} still writes to {legacy}"


@pytest.mark.parametrize("name", WORKFLOW_NAMES)
def test_filenames_do_not_repeat_the_directory_keys(name, workflow_sources):
    """
    `tandem/M1/calls/{contig}.tsv`, not
    `tandem_repeat_indels/.../{contig}_M1_tr_indels.tsv`.
    """
    code = strip_comments(workflow_sources[name])
    assert "_tr_" not in code
    assert "_remapping.fastq" not in code
    assert "_small_indels" not in code
    assert "_haploid_combined" not in code


@pytest.mark.parametrize("name", WORKFLOW_NAMES)
def test_fasplit_root_is_gone(name, workflow_sources):
    """Per-contig FASTAs moved under work/{species}/{sample}/contigs/."""
    assert "fasplit/" not in workflow_sources[name].lower()


def test_fasplit_root_removed_from_config(cfg):
    assert "fasplit_root" not in cfg["paths"]


def test_results_leaf_names_match_across_pipelines(workflow_sources):
    """
    Both pipelines write the same four leaf names, which is what lets one
    recursive glob in R load every sample and both pipelines at once.
    """
    for name, prefix in [("indel_mutation_call", "results"),
                         ("repetitive_indel_call", "results")]:
        source = workflow_sources[name]
        for leaf in ("indels.tsv", "opportunities.tsv", "filter_report.tsv", "plots.pdf"):
            assert f"{{{prefix}}}/{leaf}" in source, f"{name} does not write {leaf}"


def test_large_intermediates_are_declared_as_outputs(workflow_sources):
    """
    The concatenated FASTQ and the remapped BAM used to be written but never
    declared, so GWF could neither track nor clean the two biggest artifacts.
    """
    for name in ("indel_mutation_call", "repetitive_indel_call"):
        source = workflow_sources[name]
        assert "outputs=[combined_fastq, remapped_bam, haploid_locations]" in source, (
            f"{name} does not declare its FASTQ and BAM"
        )


def test_filter_report_is_declared_as_an_output(workflow_sources):
    for name in ("indel_mutation_call", "repetitive_indel_call"):
        assert '{results}/filter_report.tsv"' in workflow_sources[name]


def test_opportunities_get_species_and_sample_columns(workflow_sources):
    """
    Both halves of a rate must be self-describing, or opportunities cannot be
    concatenated across samples the way the indel tables can.
    """
    for name in ("indel_mutation_call", "repetitive_indel_call"):
        source = workflow_sources[name]
        assert "{merged_opportunities}.tmp > {merged_opportunities}" in source, (
            f"{name} does not append species/sample to the opportunity table"
        )


@pytest.mark.parametrize("name", WORKFLOW_NAMES)
def test_no_misspelled_alignment_variable(name, workflow_sources):
    """The original carried `alignemt` in 8 places; resolving that was step one."""
    assert "alignemt" not in workflow_sources[name]
