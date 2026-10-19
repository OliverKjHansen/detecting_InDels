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
