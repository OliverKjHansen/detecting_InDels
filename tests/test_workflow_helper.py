"""
The helper every workflow loads config.yaml with: samples, paths built from
the `layout` templates, and the contigs long enough to analyse. Run against
the fixture config of workflow_config.py.
"""

import gzip
import os
import subprocess

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
    assert pipeline.path("coverage", "ape", "s1") == f"{root}/output/precompute/ape/s1/avg_coverage.txt"
    assert pipeline.path("contig_fasta", "monkey", "s3", contig="c2") == \
        f"{root}/output/precompute/monkey/s3/contigs/c2.fa"


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
    assert pipeline.folder("contig_fasta", "ape", "s1") == f"{root}/output/precompute/ape/s1/contigs"
    assert pipeline.folder("unresolved_window", "ape", "s1") == \
        f"{root}/output/precompute/ape/s1/unresolved_haplotypes/chunks"


def test_relative_paths_can_resolve_against_another_folder(tmp_path, monkeypatch):
    root = os.path.realpath(tmp_path / "elsewhere")
    os.makedirs(root)
    write_config(root)
    monkeypatch.chdir(tmp_path)
    from pipeline import Pipeline
    pipeline = Pipeline(f"{root}/config.yaml", base=root)
    assert pipeline.path("coverage", "ape", "s1") == f"{root}/output/precompute/ape/s1/avg_coverage.txt"


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


def test_window_counts_are_summed_on_every_other_column(tmp_path):
    import subprocess
    from pipeline import summed_counts
    first, second, out = tmp_path / "a.tsv", tmp_path / "b.tsv", tmp_path / "sum.tsv"
    first.write_text("contig\tcontext\tmapQ\ttissue\tkmerQ\tcount\nc1\tA\t60\tliver\t93\t5\nc1\tA\t60\tbrain\t93\t1\n")
    second.write_text("contig\tcontext\tmapQ\ttissue\tkmerQ\tcount\nc1\tA\t60\tliver\t93\t2\nc1\tC\t5\tliver\t40\t7\n")
    subprocess.run(summed_counts([first, second], out), shell=True, check=True)
    header, *rows = out.read_text().splitlines()
    assert header == "contig\tcontext\tmapQ\ttissue\tkmerQ\tcount"
    assert sorted(rows) == ["c1\tA\t60\tbrain\t93\t1", "c1\tA\t60\tliver\t93\t7", "c1\tC\t5\tliver\t40\t7"]


def test_the_key_parameters_are_the_strata_block_as_one_shell_word(pipeline):
    import json
    import shlex
    words = shlex.split(f"--key-parameters {pipeline.key_parameters()}")
    assert len(words) == 2
    assert json.loads(words[1]) == pipeline.config["strata"]


def test_final_filtering_gets_the_platform_and_the_filters_of_its_workflow():
    import shlex
    from pipeline import final_filtering_options, shell
    config = {"sequel_ii_samples": {"ape": ["s2"]},
              "final_filters": {"keep": {"kmerQ_pass": ["pass"]}, "repeat_keep": {"read_end_distance_bin": [">=1000"]},
                                "capped_quality": 40}}
    revio = final_filtering_options(config, "ape", "s1")
    assert revio == ["--platform", "Revio", "--capped-quality", "40", "--keep", "kmerQ_pass=pass"]
    sequel = final_filtering_options(config, "ape", "s2", repeats=True)
    assert sequel[:2] == ["--platform", "Sequel II"]
    assert sequel[-2:] == ["--keep", "read_end_distance_bin=>=1000"]
    # The shell sees each argument as one word, so ">=" is no redirection.
    assert shlex.split(shell(sequel)) == sequel


def test_gzipped_window_counts_are_summed_into_a_gzipped_table(tmp_path):
    from pipeline import summed_counts
    header = "contig\tcontext\ttissue\tcount\n"
    windows = [tmp_path / "a.tsv.gz", tmp_path / "b.tsv.gz", tmp_path / "c.tsv.gz"]
    windows[0].write_bytes(gzip.compress((header + "c1\tA\tliver\t5\nc1\tC\tliver\t1\nc1\tA\tliver\t2\n").encode()))
    windows[1].write_bytes(gzip.compress((header + "c1\tC\tbrain\t7\nc1\tA\tliver\t3\n").encode()))
    windows[2].write_bytes(gzip.compress(header.encode()))   # a window without opportunities
    out = tmp_path / "sum.tsv.gz"
    subprocess.run(summed_counts(windows, out), shell=True, check=True)
    header_line, *rows = gzip.decompress(out.read_bytes()).decode().splitlines()
    assert header_line == header.strip()
    assert sorted(rows) == ["c1\tA\tliver\t10", "c1\tC\tbrain\t7", "c1\tC\tliver\t1"]


def as_gwf_runs_it(command, tmp_path):
    """Run a target's command as gwf's Slurm script does: bash, stopping at the first failing command."""
    marker = tmp_path / "went_on"
    return subprocess.run(["bash", "-c", f"set -e\n{command}\ntouch {marker}"], capture_output=True), marker.exists()


def joins(tmp_path, windows):
    """Every command that joins window or contig files, writing to tmp_path, by name."""
    from pipeline import joined_reads, summed_counts, tables_with_sample
    return {"summed_counts": summed_counts(windows, tmp_path / "sum.tsv.gz"),
            "joined_reads": joined_reads(windows, tmp_path / "reads.fastq.gz"),
            "tables_with_sample": tables_with_sample(windows, tmp_path / "merged.tsv.gz", "ape", "s1")}


@pytest.mark.parametrize("join", ["summed_counts", "joined_reads", "tables_with_sample"])
def test_a_join_stops_the_job_when_a_window_cant_be_read(tmp_path, join):
    # A window cut short, as a job killed while writing it leaves it.
    content = gzip.compress(("contig\tcount\n" + "".join(f"c1\t{i}\n" for i in range(20000))).encode())
    good, cut = tmp_path / "good.tsv.gz", tmp_path / "cut.tsv.gz"
    good.write_bytes(gzip.compress(b"contig\tcount\nc1\t1\n"))
    cut.write_bytes(content[:len(content) // 2])
    done, went_on = as_gwf_runs_it(joins(tmp_path, [good, cut])[join], tmp_path)
    assert done.returncode != 0 and not went_on


@pytest.mark.parametrize("join", ["summed_counts", "tables_with_sample"])
def test_a_join_of_large_windows_runs_to_the_end(tmp_path, join):
    # Bigger than a pipe holds, so a reader that stops after the header
    # would kill the decompression of the rest.
    rows = "".join(f"c1\tk{i % 1000}\t1\n" for i in range(300000))
    windows = [tmp_path / "a.tsv.gz", tmp_path / "b.tsv.gz"]
    for window in windows:
        window.write_bytes(gzip.compress(("contig\tkey\tcount\n" + rows).encode()))
    done, went_on = as_gwf_runs_it(joins(tmp_path, windows)[join], tmp_path)
    assert done.returncode == 0 and went_on, done.stderr
