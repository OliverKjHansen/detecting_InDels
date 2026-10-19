"""
The mirror exporter end to end, on a fake cluster checkout: the fixture data put
where workflows 1 and 2 leave their outputs, and workflows 3 and 4 run on it
with the pipeline's own code. Windows exported from it in one run must replay
locally to the same tables, and a replay must notice when they differ.
"""

import contextlib
import io
import os
import runpy
import shutil
import subprocess
import sys
from functools import partial

import pysam
import pytest
import yaml

import make_fixtures as mf
import utilities
from pipeline import Pipeline
from workflow_config import REPO_ROOT

SCRIPTS = os.path.join(REPO_ROOT, "scripts")
TRACTS = {"M0": [(2150, 2300, "AC", 2)],
          "M1": [(2100, 2300, "AC", 2)],
          "M2": [(2100, 2300, "AC", 2), (2190, 2250, "ACAG", 4, 90.0)]}
MODELS = list(TRACTS)   # the models config.yaml lists
SAMPLES = [("testspecies", "testsample"), ("otherspecies", "othersample")]
CHUNK = f"{mf.CONTIG}_0_{mf.CONTIG_LENGTH}"   # the contig is shorter than one scan window


def run_script(script, *args, stdout=None):
    """
    Run one of the pipeline's filtering scripts in this process, as from the
    command line: much faster than a new Python importing pandas each time.
    """
    argv, sys.argv = sys.argv, [script, *args]
    try:
        with contextlib.redirect_stdout(stdout or io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            runpy.run_path(os.path.join(SCRIPTS, script), run_name="__main__")
    finally:
        sys.argv = argv


def write_haploid_records(bam, calls, path):
    """
    locations.tsv for the reads carrying a call: the fixture's haploid
    assembly is the diploid one, so each read aligns to it as it does there.
    """
    with open(calls) as handle:
        column = handle.readline().rstrip("\n").split("\t").index("readID")
        names = {line.split("\t")[column] for line in handle}
    with pysam.AlignmentFile(bam) as reads, open(path, "w") as out:
        for read in sorted(reads.fetch(mf.CONTIG), key=lambda read: read.query_name):
            if read.query_name in names:
                out.write(f"{read.query_name}\t{read.flag}\t{read.reference_name}\t"
                          f"{read.reference_start + 1}\t{read.mapping_quality}\t{read.cigarstring}\n")


def haploid_check(calls, locations, reads, merged, species, sample):
    """The haploid-check step, with species and sample appended as the workflow's awk does."""
    checked = io.StringIO()
    run_script("filtering_haploid_mappping.py", calls, locations, reads, stdout=checked)
    header, *rows = checked.getvalue().splitlines()
    with open(merged, "w") as out:
        out.write(header + "\tspecies\tsample\n")
        out.writelines(f"{row}\t{species}\t{sample}\n" for row in rows)


def final_filtering(merged, final_calls, report, *extra):
    run_script("filtering_small_indels.py", merged, "-o", final_calls,
               "--keep", "kmerQ_pass=pass", *map(str, extra), "--report", report)


def made(path):
    """The path, with its folder created."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


def place(source, target):
    os.makedirs(os.path.dirname(target), exist_ok=True)
    shutil.copy(source, target)


@pytest.fixture(scope="module")
def cluster(tmp_path_factory, fixture_paths):
    """A checkout on the cluster that has run workflows 1-4 for both samples."""
    root = os.path.realpath(tmp_path_factory.mktemp("cluster"))
    # A checkout: the replays run its scripts unless told otherwise.
    os.symlink(SCRIPTS, f"{root}/scripts")
    os.symlink(os.path.join(REPO_ROOT, "workflows"), f"{root}/workflows")
    with open(os.path.join(REPO_ROOT, "config.yaml")) as handle:
        config = yaml.safe_load(handle)
    config.update({
        "samples": {species: [sample] for species, sample in SAMPLES},
        "tissue_samples": {"testspecies": ["testsample"]},
        "methylation_samples": {}, "sperm_samples": {}, "sample_age": {}, "alignment_overrides": {},
        "paths": {"fasta": f"{root}/fasta", "bam": f"{root}/bam", "read_labels": f"{root}/labels",
                  "classification": f"{root}/classification",
                  "output": "output", "results": "results", "plots": "plots"},
        "min_contig_length": 1000,
    })
    with open(f"{root}/config.yaml", "w") as handle:
        yaml.safe_dump(config, handle)
    pipeline = Pipeline(f"{root}/config.yaml", base=root)
    filters = pipeline.config["final_filters"]
    coverage_band = pipeline.config["scan"]["coverage_band"]

    for species, sample in SAMPLES:
        path = partial(pipeline.path, species=species, sample=sample)
        assembly, bam, coverage = path("assembly"), pipeline.alignment(species, sample), path("coverage")
        tissue = path("tissue") if pipeline.is_listed("tissue_samples", species, sample) else "NA"

        # Workflows 1 and 2.
        for suffix in ("", ".fai"):
            place(fixture_paths["assembly"] + suffix, assembly + suffix)
        for suffix in ("", ".bai"):
            place(fixture_paths["bam"] + suffix, bam + suffix)
        place(fixture_paths["coverage"], coverage)
        place(fixture_paths["filters"].replace("_filters.tsv", "_unresolved_haplotypes.bed"),
              path("unresolved", contig=mf.CONTIG))
        if tissue != "NA":
            place(fixture_paths["tissue"], tissue)
        for model in MODELS:
            mf.write_tr_bed(made(path("tandem_repeats", contig=mf.CONTIG, model=model)), TRACTS[model], model)

        # Workflow 3.
        mask = made(path("nonrepetitive_mask", contig=mf.CONTIG))
        place(fixture_paths["filters"], mask)
        window = dict(contig=mf.CONTIG, start=0, end=mf.CONTIG_LENGTH)
        with open(made(path("nonrepetitive_window_calls", **window)), "w") as calls, \
                open(made(path("nonrepetitive_window_opportunities", **window)), "w") as opportunities, \
                open(made(path("nonrepetitive_window_reads", **window)), "w") as fastq:
            utilities.nonrepetitive_scan(mask, mf.CONTIG, bam, assembly, coverage, tissue, 0, mf.CONTIG_LENGTH,
                                         calls_out=calls, opportunities_out=opportunities, fastq_out=fastq,
                                         coverage_band=coverage_band,
                                         indel_sizes=pipeline.config["scan"]["indel_sizes"],
                                         key_parameters=pipeline.config["strata"],
                                         flank_width=pipeline.config["scan"]["flank_width"])
        calls = made(path("nonrepetitive_calls", contig=mf.CONTIG))
        place(path("nonrepetitive_window_calls", **window), calls)
        locations = made(path("nonrepetitive_haploid_locations"))
        write_haploid_records(bam, calls, locations)
        merged = made(path("nonrepetitive_merged_calls"))
        haploid_check(calls, locations, path("nonrepetitive_window_reads", **window), merged, species, sample)
        final_filtering(merged, made(path("nonrepetitive_results_calls")), path("nonrepetitive_results_report"))

        # Workflow 4.
        with contextlib.ExitStack() as files:
            def opened(name, **fields):
                return files.enter_context(open(made(path(name, **window, **fields)), "w"))
            utilities.tandem_repeat_scan(
                {model: path("tandem_repeats", contig=mf.CONTIG, model=model) for model in MODELS},
                mf.CONTIG, bam, assembly, coverage, tissue, 0, mf.CONTIG_LENGTH,
                calls_out={model: opened("repetitive_window_calls", model=model) for model in MODELS},
                opportunities_out={model: opened("repetitive_window_opportunities", model=model) for model in MODELS},
                fastq_out=opened("repetitive_window_reads"),
                coverage_band=coverage_band, key_parameters=pipeline.config["strata"],
                flank_width=pipeline.config["scan"]["flank_width"])
        locations = made(path("repetitive_haploid_locations"))
        all_calls = f"{os.path.dirname(locations)}/all_calls.tsv"
        with open(all_calls, "w") as out:
            for i, model in enumerate(MODELS):
                calls = made(path("repetitive_calls", contig=mf.CONTIG, model=model))
                place(path("repetitive_window_calls", **window, model=model), calls)
                with open(calls) as handle:
                    out.writelines(line for j, line in enumerate(handle) if i == 0 or j > 0)
        write_haploid_records(bam, all_calls, locations)
        for model in MODELS:
            merged = made(path("repetitive_merged_calls", model=model))
            haploid_check(path("repetitive_calls", contig=mf.CONTIG, model=model), locations,
                          path("repetitive_window_reads", **window), merged, species, sample)
            final_filtering(merged, made(path("repetitive_results_calls", model=model)),
                            path("repetitive_results_report", model=model),
                            "--min-read-end-distance", filters["min_read_end_distance"])
    return root


@pytest.fixture(scope="module")
def windows(cluster):
    """Both samples exported in one run: one with a scan-job window, one with --window-size."""
    subprocess.run([sys.executable, os.path.join(SCRIPTS, "export_local_subset.py"),
                    f"testspecies/testsample:{mf.CONTIG}:0-{mf.CONTIG_LENGTH}", "otherspecies/othersample",
                    "--window-size", "3000", "--out", "mirror"], cwd=cluster, check=True, capture_output=True)
    return {
        "scan job": f"{cluster}/mirror/testspecies/testsample/{CHUNK}",
        # The fixture's unique calls are at 500-1500, so the busiest 3 kb window is the first.
        "smaller": f"{cluster}/mirror/otherspecies/othersample/{mf.CONTIG}_0_3000",
    }


def replay(window, repo=None):
    """
    Run a window's run_local.sh: (exit code, output). Without `repo` it runs
    the scripts of the checkout the window sits in, as it does by default.
    """
    env = {key: value for key, value in os.environ.items() if key != "REPO"}
    env["PYTHON"] = sys.executable
    if repo:
        env["REPO"] = repo
    done = subprocess.run(["bash", f"{window}/run_local.sh"], capture_output=True, text=True, env=env)
    return done.returncode, done.stdout + done.stderr


def test_every_exported_window_replays_both_workflows_to_the_clusters_tables(windows):
    for kind, window in windows.items():
        code, output = replay(window)
        assert code == 0, output
        assert "DIFFERENT" not in output
        same = [line for line in output.splitlines() if line.startswith("same:")]
        # Calls, calls after the haploid check and final calls, for workflow 3
        # and each model of workflow 4; opportunities too for a whole scan job.
        assert len(same) == (4 * 4 if kind == "scan job" else 3 * 4), output


def test_a_replay_that_differs_from_the_cluster_fails(windows, tmp_path):
    # The folder's name says which window it holds.
    window = shutil.copytree(windows["scan job"], tmp_path / CHUNK)
    final = window / "cluster" / "repetitive_indels" / "M1" / "indels.tsv"
    header, *rows = final.read_text().splitlines()
    assert rows, "the fixture should leave final calls inside the repeat tract"
    final.write_text("\n".join([header] + rows[1:]) + "\n")
    code, output = replay(window, repo=REPO_ROOT)
    assert code == 1
    assert [line for line in output.splitlines() if line.startswith("DIFFERENT")] == \
        ["DIFFERENT: workflow 4 M1 final calls"]


def test_the_help_says_the_mirror_is_for_software_testing_only():
    done = subprocess.run([sys.executable, os.path.join(SCRIPTS, "export_local_subset.py"), "--help"],
                          capture_output=True, text=True, check=True)
    assert "software testing only" in " ".join(done.stdout.split())
