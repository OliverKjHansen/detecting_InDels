"""
Checks on config.yaml that need no cluster.

A sample switched off in `samples` is still allowed in the other tables (its
age stays in `sample_age`, its tissue labels in `tissue_samples`), so these
check that every name is a known sample of a known species rather than that
every table is a subset of `samples`. A typo in a species or sample name is
what they catch: the workflows would otherwise skip the entry silently.
"""

import os

import pytest
import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SAMPLE_TABLES = ["samples", "tissue_samples", "methylation_samples", "sperm_samples"]


@pytest.fixture(scope="module")
def cfg():
    with open(os.path.join(REPO_ROOT, "config.yaml")) as handle:
        return yaml.safe_load(handle)


def known_samples(cfg, species):
    """Samples switched on in `samples`, plus those only on the sample sheet."""
    return set(cfg["samples"].get(species) or []) | set((cfg.get("sample_age") or {}).get(species) or {})


@pytest.mark.parametrize("table", SAMPLE_TABLES)
def test_every_sample_belongs_to_a_known_species(cfg, table):
    for species, samples in (cfg.get(table) or {}).items():
        assert species in cfg["samples"], f"{table}: unknown species {species}"
        # A bare string would be iterated character by character.
        assert isinstance(samples, list), f"{table}.{species} must be a list"
        assert len(samples) == len(set(samples)), f"{table}.{species} lists a sample twice"
        unknown = set(samples) - known_samples(cfg, species)
        assert not unknown, f"{table}.{species}: {sorted(unknown)} are not samples of {species}"


def test_every_age_belongs_to_a_known_species(cfg):
    for species, ages in cfg["sample_age"].items():
        assert species in cfg["samples"], f"sample_age: unknown species {species}"
        for sample, age in ages.items():
            assert isinstance(age, (int, float)) and age >= 0, f"sample_age.{species}.{sample}: {age!r}"


def test_override_paths_are_bams_of_switched_on_samples(cfg):
    for species, overrides in (cfg.get("alignment_overrides") or {}).items():
        for sample, path in overrides.items():
            assert sample in (cfg["samples"].get(species) or []), (
                f"alignment_overrides.{species}.{sample} is not in samples")
            assert isinstance(path, str) and path.endswith(".bam"), (
                f"alignment_overrides.{species}.{sample} is not a BAM: {path!r}")
            assert path.startswith(("/", "~/")), (
                f"alignment_overrides.{species}.{sample} is not an absolute path: {path!r}")
            assert not any(c.isspace() for c in path), (
                f"alignment_overrides.{species}.{sample} contains whitespace: {path!r}")


def test_layout_templates_resolve_against_the_declared_roots(cfg):
    """
    Every placeholder in a layout template is either a root under `paths` or
    one of the fields a workflow fills in.
    """
    roots = {key: "ROOT" for key in cfg["paths"]}
    fields = {"species": "SPECIES", "sample": "SAMPLE", "contig": "CONTIG", "model": "MODEL",
              "start": 0, "end": 1, "k": 2, "assembly_type": cfg["assembly_type"]}
    for key, template in cfg["layout"].items():
        try:
            template.format(**roots, **fields)
        except KeyError as exc:
            pytest.fail(f"layout.{key} references undefined placeholder {exc}")


def test_output_templates_are_sample_keyed(cfg):
    """Everything the pipeline writes lives under {species}/{sample}."""
    for key, template in cfg["layout"].items():
        if template.startswith("{output}/"):
            assert "{species}/{sample}" in template, f"layout.{key} is not keyed by species and sample"


def test_min_contig_length_is_a_positive_int(cfg):
    assert isinstance(cfg["min_contig_length"], int)
    assert cfg["min_contig_length"] > 0
