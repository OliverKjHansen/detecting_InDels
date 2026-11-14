"""
Checks on config.yaml that need no cluster.

A sample switched off in `samples` is still allowed in the other tables (its
age stays in `sample_age`, its tissue labels in `tissue_samples`), so these
check that every name is a known sample of a known species rather than that
every table is a subset of `samples`. A typo in a species or sample name is
what they catch: the workflows would otherwise skip the entry silently.
"""

import glob
import json
import os

import pytest
import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SAMPLE_TABLES = ["samples", "tissue_samples", "methylation_samples", "sperm_samples", "sequel_ii_samples"]


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


def test_new_data_share_the_assembly_of_a_switched_on_sample_of_their_species(cfg):
    # The new data use the sample's contigs and repeats, which only that
    # sample's own targets make (workflows 1 and 2).
    for species, shared in (cfg.get("assembly_of") or {}).items():
        assert species in cfg["samples"], f"assembly_of: unknown species {species}"
        assert isinstance(shared, dict), f"assembly_of.{species} must map new data to a sample"
        for new, sample in shared.items():
            assert new in known_samples(cfg, species), f"assembly_of.{species}: {new} is not a sample of {species}"
            assert sample != new, f"assembly_of.{species}.{new} names itself"
            assert sample not in shared, f"assembly_of.{species}.{new}: {sample} shares another's assembly itself"
            if new in (cfg["samples"].get(species) or []):
                assert sample in cfg["samples"][species], (
                    f"assembly_of.{species}.{new} is switched on, but {sample}, whose assembly it uses, is not")


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
              "haplotype": "HAPLOTYPE", "start": 0, "end": 1, "k": 2, "assembly_type": cfg["assembly_type"]}
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


# -- every key is read ---------------------------------------------------------


class Tracked(dict):
    """A config mapping that records the dotted names of the keys read from it."""

    def __init__(self, data, read, prefix=""):
        super().__init__(data)
        self.read, self.prefix = read, prefix

    def _wrap(self, key, value):
        name = f"{self.prefix}{key}"
        self.read.add(name)
        return Tracked(value, self.read, name + ".") if isinstance(value, dict) else value

    def __getitem__(self, key):
        return self._wrap(key, super().__getitem__(key))

    def get(self, key, default=None):
        return self._wrap(key, super().__getitem__(key)) if key in self else default

    def items(self):
        return [(key, self._wrap(key, value)) for key, value in super().items()]


def keys_to_check(cfg):
    """
    The top-level keys, and the keys inside settings blocks such as `layout`
    or `scan`. A table keyed by species (samples, sample_age, ...) counts as
    one key.
    """
    names = set()
    for key, value in cfg.items():
        names.add(key)
        if isinstance(value, dict) and not set(value) <= set(cfg["samples"]):
            names |= {f"{key}.{inner}" for inner in value}
    return names


def test_every_config_key_is_read_by_a_workflow_or_a_notebook(cfg, tmp_path, monkeypatch):
    import pipeline
    from test_workflow_snapshots import WORKFLOWS, load_targets
    from workflow_config import write_config

    root = os.path.realpath(tmp_path)
    write_config(root)
    read = set()
    load = pipeline.yaml.safe_load
    monkeypatch.setattr(pipeline.yaml, "safe_load", lambda handle: Tracked(load(handle), read))
    for name in WORKFLOWS:
        load_targets(name, root, monkeypatch)

    # The roots under `paths` are read by every Pipeline; a root is used if a
    # layout template that was read refers to it.
    read = {name for name in read if not name.startswith("paths.")}
    templates = [cfg["layout"][name.split(".", 1)[1]] for name in read if name.startswith("layout.")]
    read |= {f"paths.{root_name}" for root_name in cfg["paths"]
             if any(f"{{{root_name}}}" in template for template in templates)}

    notebooks = ""
    for notebook in glob.glob(os.path.join(REPO_ROOT, "plotting", "*.ipynb")):
        with open(notebook) as handle:
            notebooks += "".join("".join(cell["source"]) for cell in json.load(handle)["cells"]
                                 if cell["cell_type"] == "code")
    read |= {key for key in cfg if f'"{key}"' in notebooks}

    unread = sorted(keys_to_check(cfg) - read)
    assert not unread, f"config keys nothing reads: {unread}"


def test_the_final_filters_keep_values_their_keys_can_take(cfg):
    """
    A kept value no row can have would drop every row: each must be a label
    of the key's bins in `strata`, so changing the bins means changing the
    filter too.
    """
    import strata
    bins = {"read_end_distance_bin": "read_end_distance_bins", "read_length_bin": "read_length_bins",
            "mapQ_bin": "mapq_bins", "subread_support": "subread_support_bins",
            "tract_copies_bin": "tract_copies_bins", "tract_purity_bin": "tract_purity_bins"}
    bins |= {"inner_flank_q": "flank_quality_bins", "outer_flank_q": "flank_quality_bins"}
    yes_no = {"yes", "no", "NA"}
    labels = {"kmerQ_pass": {"pass", "fail"}, "junction_disagreement": yes_no, "repeats_next": yes_no,
              "base_after_masked": yes_no, "inserted_bases_pass": {"pass", "fail", "NA"},
              "haploid_check": {"True", "False", "NA"}, "hap1_check": {"True", "False", "NA"},
              "hap2_check": {"True", "False", "NA"}}
    for column, name in bins.items():
        edges = cfg["strata"][name]
        bin_of = strata.binned(edges)
        labels[column] = {bin_of(value) for value in [edges[0] - 1, *edges]}
    filters = cfg["final_filters"]
    for column, values in {**(filters["keep"] or {}), **filters["nonrepetitive_keep"],
                           **filters["nonrepetitive_calls_only"], **filters["repeat_keep"]}.items():
        assert column in labels, f"final_filters.{column}: add its labels to this test"
        assert set(map(str, values)) <= labels[column], f"final_filters.{column}: {values} are not all its labels"
