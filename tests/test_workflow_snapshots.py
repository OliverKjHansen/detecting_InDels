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


@pytest.mark.parametrize("budget, groups", [(500_000_000, [["c1", "c2"]]), (10_000_000, [["c1"], ["c2"]])])
def test_workflow_3_merges_each_contig_once_in_groups_up_to_the_budget(tmp_path, monkeypatch, budget, groups):
    # Contigs adding up to scan.merge_bp_per_job are merged in one job (ticket
    # 49). The opportunities are summed in jobs of their own (ticket 40), so
    # changing how leaves the calls and reads, and the remapping after them,
    # as they are. ape/s1's contigs: c1 (12 Mb) and c2 (1.5 Mb).
    import yaml
    root = os.path.realpath(tmp_path)
    write_config(root)
    with open(os.path.join(root, "config.yaml")) as handle:
        config = yaml.safe_load(handle)
    config["scan"]["merge_bp_per_job"] = budget
    with open(os.path.join(root, "config.yaml"), "w") as handle:
        yaml.safe_dump(config, handle)
    targets = load_targets("indel_mutation_call", root, monkeypatch)
    for kind, suffixes in (("non_rep_MergeWindows_", [".fastq.gz", ".tsv"]), ("non_rep_MergeOpportunities_", [".parquet"])):
        jobs = {name: target for name, target in targets.items() if name.startswith(kind) and name.endswith("_ape_s1_diploid")}
        made = [sorted(os.path.basename(path) for path in target["outputs"]) for target in jobs.values()]
        assert sorted(made) == sorted(sorted(f"{contig}{suffix}" for contig in group for suffix in suffixes)
                                      for group in groups)
        for target in jobs.values():
            contigs = {os.path.basename(path).split(".")[0] for path in target["outputs"]}
            assert target["inputs"] and all(os.path.basename(path).split("_")[0] in contigs for path in target["inputs"])
            if kind == "non_rep_MergeOpportunities_":
                assert all("/opportunities/chunks/" in path for path in target["inputs"])


def test_workflow_3_sums_a_tissue_samples_opportunities_in_groups_up_to_its_own_budget(tmp_path, monkeypatch):
    # pink_eric's sums of 500 Mb of contigs timed out at 4 h (ticket 55); its
    # calls and reads merged in time. ape/s1 is a tissue sample: c1 (12 Mb)
    # and c2 (1.5 Mb).
    import yaml
    root = os.path.realpath(tmp_path)
    write_config(root)
    with open(os.path.join(root, "config.yaml")) as handle:
        config = yaml.safe_load(handle)
    config["scan"].update(merge_bp_per_job=500_000_000, tissue_sum_bp_per_job=12_000_000)
    with open(os.path.join(root, "config.yaml"), "w") as handle:
        yaml.safe_dump(config, handle)
    targets = load_targets("indel_mutation_call", root, monkeypatch)

    def groups(kind):
        return sorted(sorted(os.path.basename(path).split(".")[0] for path in target["outputs"])
                      for name, target in targets.items() if name.startswith(kind) and name.endswith("_ape_s1_diploid"))

    assert groups("non_rep_MergeWindows_") == [["c1", "c1", "c2", "c2"]]   # calls and reads of each
    assert groups("non_rep_MergeOpportunities_") == [["c1"], ["c2"]]


def test_workflow_4_sums_the_opportunities_in_their_own_job(tmp_path, monkeypatch):
    root = os.path.realpath(tmp_path)
    write_config(root)
    targets = load_targets("repetitive_indel_call", root, monkeypatch)
    merge = targets["TRMergeWindows_c1_ape_s1"]
    assert not any(path.endswith(".parquet") for path in merge["outputs"])
    opportunities = targets["TRMergeOpportunities_c1_ape_s1"]
    assert sorted(path.split("/")[-3:] for path in opportunities["outputs"]) == [
        [model, "opportunities", "c1.parquet"] for model in ("M0", "M1", "M2")]
    assert opportunities["inputs"] and all("/opportunities/chunks/c1_" in path for path in opportunities["inputs"])


def test_workflow_3_rates_each_sample_by_filter_from_its_candidate_tables_with_final_filterings_steps(tmp_path, monkeypatch):
    # Ticket 56: the calls and opportunities before final filtering; final
    # filtering's own steps, then the funnel's (filter_rates in config.yaml).
    import shlex
    root = os.path.realpath(tmp_path)
    config = write_config(root)
    targets = load_targets("indel_mutation_call", root, monkeypatch)
    rates = targets["non_rep_FilterRates_ape_s1_diploid"]
    final = targets["non_rep_FinalFiltering_ape_s1_diploid"]
    assert sorted(path.split("/", 1)[1] for path in rates["outputs"]) == [
        "results/nonrepetitive_indels/ape/s1/rates_by_filter.tsv", "results/nonrepetitive_indels/ape/s1/rates_funnel.tsv"]
    candidates = [path for path in final["inputs"] if "/merged/" in path]
    assert len(candidates) == 2 and set(candidates) <= set(rates["inputs"])

    def options(spec, name):
        words = shlex.split(spec)
        return [words[i + 1] for i, word in enumerate(words) if word == name]

    for name in ("--keep", "--calls-only", "--platform", "--capped-quality", "--smrt-cells"):
        assert options(rates["spec"], name) == options(final["spec"], name), name
    assert options(rates["spec"], "--step") == [f"{column}={','.join(map(str, values))}" for column, values
                                                in config["filter_rates"]["nonrepetitive_steps"].items()]
    assert options(rates["spec"], "--step")
