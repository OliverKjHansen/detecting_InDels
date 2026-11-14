"""
New data of a sample sequenced again: their own BAM, aligned to the sample's
assembly (config.yaml `assembly_of`). The new data are a sample of their own
from the reads on (coverage, SMRT Cells, unresolved haplotypes, the scans);
what follows from the assembly alone (its contigs, its haplotypes, workflow 2's
repeat annotation) is the original sample's, made once.

Loaded the way gwf loads the workflows, with the fixture config of
workflow_config.py plus ape/s1_new_data, sharing ape/s1's assembly; it has no
assembly index of its own.
"""

import os
from collections import Counter

import pytest
import yaml

from test_workflow_snapshots import WORKFLOWS, load_targets
from workflow_config import write_config

NEW = "s1_new_data"


@pytest.fixture
def root(tmp_path):
    root = os.path.realpath(tmp_path)
    config = write_config(root)
    config["samples"] = {**config["samples"], "ape": config["samples"]["ape"] + [NEW]}   # not workflow_config's own
    config["assembly_of"] = {"ape": {NEW: "s1"}}
    with open(os.path.join(root, "config.yaml"), "w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    return root


def of_new_data(targets):
    return {name: target for name, target in targets.items() if f"_ape_{NEW}" in name}


def test_new_data_are_scanned_with_their_own_reads_on_the_assembly_they_share(root, monkeypatch):
    scans = {name: target for name, target in of_new_data(load_targets("indel_mutation_call", root, monkeypatch)).items()
             if name.startswith("non_rep_Scan_")}
    assert scans
    for target in scans.values():
        assert "{ROOT}/fasta/ape/s1/diploid_assembly/s1_diploid.fa" in target["inputs"]
        assert f"{{ROOT}}/bam/ape/{NEW}/kinetics/{NEW}_kinetics_diploid.bam" in target["inputs"]
        assert f"{{ROOT}}/output/precompute/ape/{NEW}/avg_coverage.txt" in target["inputs"]


@pytest.mark.parametrize("workflow", WORKFLOWS)
def test_no_two_targets_make_the_same_file(root, monkeypatch, workflow):
    made = Counter(path for target in load_targets(workflow, root, monkeypatch).values() for path in target["outputs"])
    assert not sorted(path for path, times in made.items() if times > 1)


def test_the_shared_assembly_is_split_and_annotated_once_for_its_sample(root, monkeypatch):
    precompute = of_new_data(load_targets("precompute", root, monkeypatch))
    # From the reads: the new data's own.
    assert {"average_coverage", "smrt_cells", "unresolved_haplotype"} <= {name.split("_ape_")[0] for name in precompute}
    assert not [name for name in precompute if name.startswith(("faSplit_", "haplotypes_"))]
    assert not of_new_data(load_targets("repetitive_dna", root, monkeypatch))


def test_the_repeat_files_new_data_read_are_made_by_workflow_2(root, monkeypatch):
    made = {path for target in load_targets("repetitive_dna", root, monkeypatch).values() for path in target["outputs"]}
    read = {path for workflow in ("indel_mutation_call", "repetitive_indel_call")
            for target in of_new_data(load_targets(workflow, root, monkeypatch)).values()
            for path in target["inputs"] if "/repetitive_dna/" in path}
    assert read
    assert read <= made, sorted(read - made)[:5]
