"""
What each workflow would run, recorded so a refactor can show it changed nothing.

Each workflow is loaded the way gwf loads it, with the fixture config of
workflow_config.py, and every target's name, inputs, outputs, options and
command are compared with the snapshot in tests/snapshots/. Paths under the
fixture root are written as {ROOT}.

After a deliberate change, rewrite the snapshots and review their diff:

    UPDATE_SNAPSHOTS=1 .venv/bin/python -m pytest tests/test_workflow_snapshots.py
"""

import difflib
import json
import os
import runpy

import pytest

from workflow_config import REPO_ROOT, write_config

WORKFLOW_DIR = os.path.join(REPO_ROOT, "workflows")
SNAPSHOT_DIR = os.path.join(REPO_ROOT, "tests", "snapshots")
WORKFLOWS = ["precompute", "repetitive_dna", "indel_mutation_call",
             "repetitive_indel_call", "classification_denominator"]


def load_targets(name, root, monkeypatch):
    """Every target of workflows/{name}.py, loaded from `root` as gwf would run it."""
    monkeypatch.chdir(root)
    # gwf puts the workflow's folder on the path, so it can import its helpers.
    monkeypatch.syspath_prepend(WORKFLOW_DIR)
    gwf = runpy.run_path(os.path.join(WORKFLOW_DIR, f"{name}.py"))["gwf"]

    def relative(value):
        return value.replace(root, "{ROOT}").replace(REPO_ROOT, "{REPO}")

    return {
        target.name: {
            "inputs": [relative(path) for path in target.flattened_inputs()],
            "outputs": [relative(path) for path in target.flattened_outputs()],
            "options": target.options,
            "spec": relative(target.spec),
        }
        for target in gwf.targets.values()
    }


def describe_changes(recorded, current):
    """A readable account of how two snapshots differ."""
    lines = []
    for name in sorted(set(recorded) - set(current)):
        lines.append(f"removed: {name}")
    for name in sorted(set(current) - set(recorded)):
        lines.append(f"added:   {name}")
    for name in sorted(set(recorded) & set(current)):
        for field in ("inputs", "outputs", "options", "spec"):
            before, after = recorded[name][field], current[name][field]
            if before == after:
                continue
            lines.append(f"changed: {name} {field}")
            if field == "spec":
                lines += difflib.unified_diff(before.splitlines(), after.splitlines(), lineterm="", n=1)
            else:
                lines += [f"    was {before}", f"    now {after}"]
    return "\n".join(lines)


@pytest.mark.parametrize("name", WORKFLOWS)
def test_workflow_targets_match_the_snapshot(name, tmp_path, monkeypatch):
    root = os.path.realpath(tmp_path)
    write_config(root)
    current = load_targets(name, root, monkeypatch)
    assert current, f"{name} defines no targets for the fixture config"

    snapshot = os.path.join(SNAPSHOT_DIR, f"{name}.json")
    if os.environ.get("UPDATE_SNAPSHOTS"):
        os.makedirs(SNAPSHOT_DIR, exist_ok=True)
        with open(snapshot, "w") as handle:
            json.dump(current, handle, indent=1, sort_keys=True)
            handle.write("\n")
    if not os.path.exists(snapshot):
        pytest.fail(f"no snapshot for {name}; record one with UPDATE_SNAPSHOTS=1")

    with open(snapshot) as handle:
        recorded = json.load(handle)
    assert current == recorded, f"{name} differs from its snapshot:\n" + describe_changes(recorded, current)


def test_workflow_4_reads_the_long_period_window_beds_workflow_2_writes(tmp_path, monkeypatch):
    # Step 0 lists them itself, to read the sub-windows TRF couldn't finish.
    root = os.path.realpath(tmp_path)
    write_config(root)
    written = {path for target in load_targets("repetitive_dna", root, monkeypatch).values()
               for path in target["outputs"] if "/excluded/chunks/" in path}
    read = {path for name, target in load_targets("repetitive_indel_call", root, monkeypatch).items()
            if name.startswith("TRTracts_") for path in target["inputs"] if "/excluded/chunks/" in path}
    assert written and read == written


def test_workflow_3_merges_each_contig_in_its_own_job(tmp_path, monkeypatch):
    # So contigs merge side by side, and a failure costs one contig.
    root = os.path.realpath(tmp_path)
    write_config(root)
    targets = load_targets("indel_mutation_call", root, monkeypatch)
    for contig in ("c1", "c2"):
        merge = targets[f"non_rep_MergeWindows_{contig}_ape_s1_diploid"]
        assert merge["inputs"] and all(os.path.basename(path).startswith(f"{contig}_") for path in merge["inputs"])
        assert sorted(os.path.basename(path) for path in merge["outputs"]) == \
            [f"{contig}.fastq.gz", f"{contig}.tsv", f"{contig}.tsv.gz"]
