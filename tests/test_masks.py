"""
What the InDel workflows leave out, built by the workflows' own targets on a
sample with one short contig: the excluded group, workflow 3's mask, the
tracts workflow 4 measures, and the summaries of what is left out.

Each test writes the files a target reads where the workflow before it leaves
them, runs the target's command with bash, as gwf would, and reads what it
wrote. The tests need no cluster and no external binaries: TRF and samtools
are stand-ins that print canned output.
"""

import pytest
import yaml

import targets
from targets import CONTIG, SAMPLE, SPECIES

MODELS = ["M0", "M1", "M2"]   # the models config.yaml lists
# Every base differs from the next, so the contig has no homopolymer run but
# those a test puts in.
BACKGROUND = "ACGT" * 250


def with_runs(*runs):
    """The background contig with homopolymer runs put in: (start, base, length) each."""
    sequence = list(BACKGROUND)
    for start, base, length in runs:
        sequence[start:start + length] = base * length
    return "".join(sequence)


def sample(tmp_path, monkeypatch, sequence=BACKGROUND):
    return targets.sample(tmp_path, monkeypatch, sequence)


def trf_row(start, end, motif):
    """A line of TRF's -ngs output for a pure tract at 0-based [start, end) of the window."""
    sequence = (motif * (end - start))[:end - start]
    return (f"{start + 1} {end} {len(motif)} {(end - start) / len(motif):.1f} {len(motif)} 100 0 {2 * (end - start)}"
            f" 25 25 25 25 2.00 {motif} {sequence} AAAAA CCCCC")


def annotate_excluded(sample, trf_rows=(), timed_out=False):
    """
    Workflow 2's excluded group for the contig, from its homopolymer runs and
    what TRF reports for the window holding the contig, the `trf_rows`. The
    stand-in TRF prints them and records its arguments in trf_arguments.txt.
    The stand-in timeout runs TRF, or, with timed_out, stops it as time runs
    out (status 124).
    """
    with open(f"{sample.root}/trf.txt", "w") as handle:
        handle.writelines(f"{line}\n" for line in ["@window", *trf_rows])
    sample.stand_in("trf", f'#!/bin/sh\necho "$@" > {sample.root}/trf_arguments.txt\ncat {sample.root}/trf.txt\n')
    sample.stand_in("timeout", '#!/bin/sh\nexit 124\n' if timed_out else '#!/bin/sh\nshift\nexec "$@"\n')
    sample.run("repetitive_dna", f"long_period_trf_{CONTIG}_0_{SPECIES}_{SAMPLE}_diploid")
    sample.run("repetitive_dna", f"excluded_{CONTIG}_{SPECIES}_{SAMPLE}_diploid")


def annotate_repeats(sample, tracts=(), models=MODELS):
    """
    Workflow 2's tracts of the contig, resolved from per-k BEDs that hold
    these tracts, (start, end, motif, k) each, for the models listed, and none
    for the others.
    """
    for model in MODELS:
        for k in sample.pipeline.config["repeats"]["periods"]:
            rows = [tract for tract in tracts if tract[3] == k] if model in models else []
            if model == "M1":   # TRF writes one per window
                sample.write_tracts("repeats_per_k_window", rows, model=model, contig=CONTIG, k=k, start=0)
            else:
                sample.write_tracts("repeats_per_k", rows, model=model, contig=CONTIG, k=k)
        sample.run("repetitive_dna", f"{model}_merge_k_{CONTIG}_{SPECIES}_{SAMPLE}_diploid")


def write_unresolved(sample, regions):
    """Workflow 1's unresolved-haplotype regions of the contig, (start, end) each."""
    sample.write("unresolved", [f"{CONTIG}\t{start}\t{end}" for start, end in regions], contig=CONTIG)


def measure_tracts(sample, unresolved=()):
    """The tracts workflow 4 measures, per model, with these unresolved-haplotype regions."""
    write_unresolved(sample, unresolved)
    sample.run("repetitive_indel_call", f"TRTracts_{SPECIES}_{SAMPLE}")
    return {model: sample.regions("repetitive_tracts", contig=CONTIG, model=model) for model in MODELS}


def build_mask(sample, unresolved=()):
    """Workflow 3's mask, with these unresolved-haplotype regions."""
    write_unresolved(sample, unresolved)
    sample.run("indel_mutation_call", f"non_rep_Mask_{SPECIES}_{SAMPLE}_diploid")
    return sample.regions("nonrepetitive_mask", contig=CONTIG)


# -- the excluded group (decisions Q8, Q16) ------------------------------------


def test_workflow_3_masks_a_6_bp_homopolymer_and_not_a_5_bp_one(tmp_path, monkeypatch):
    s = sample(tmp_path, monkeypatch, with_runs((200, "A", 6), (401, "C", 5)))
    annotate_excluded(s)
    annotate_repeats(s)
    assert build_mask(s) == [(200, 206)]


def test_workflow_3_masks_tracts_with_periods_above_the_models(tmp_path, monkeypatch):
    s = sample(tmp_path, monkeypatch)
    annotate_excluded(s, [trf_row(100, 145, "ACGTTGCAT"),   # period 9
                          trf_row(500, 530, "AACCGT")])     # period 6: the models' (M1's)
    annotate_repeats(s)
    assert build_mask(s) == [(100, 145)]


def test_long_periods_are_found_with_m1s_scoring_up_to_the_configured_period(tmp_path, monkeypatch):
    s = sample(tmp_path, monkeypatch)
    annotate_excluded(s)
    with open(f"{s.root}/trf_arguments.txt") as handle:
        arguments = handle.read().split()
    repeats = s.pipeline.config["repeats"]
    assert arguments[1:8] == repeats["trf_params"].split() + [str(repeats["long_period_max"])]


def test_a_short_period_reported_at_a_multiple_of_its_period_is_not_a_long_period(tmp_path, monkeypatch):
    # TRF can report an (AC)n tract again at period 8, as ACACACAC. It is M1's
    # tract, so it must not exclude that tract from workflow 4.
    s = sample(tmp_path, monkeypatch)
    annotate_excluded(s, [trf_row(300, 340, "ACACACAC"), trf_row(600, 640, "AAAAAAAA")])
    assert s.regions("long_periods", contig=CONTIG) == []


def test_workflow_4_drops_every_tract_that_overlaps_an_excluded_region_whole(tmp_path, monkeypatch):
    s = sample(tmp_path, monkeypatch, with_runs((200, "A", 6), (401, "C", 5)))
    annotate_excluded(s, [trf_row(600, 650, "ACGTTGCAT")])
    tracts = [(150, 201, "AC", 2),    # its last base is the 6 bp run's first
              (206, 240, "AG", 2),    # starts where the run ends
              (380, 420, "AC", 2),    # holds the 5 bp run
              (640, 700, "AAG", 3)]   # overlaps the period-9 tract
    annotate_repeats(s, tracts)
    assert measure_tracts(s) == {model: [(206, 240), (380, 420)] for model in MODELS}


def test_the_excluded_bases_are_summarised_per_sample(tmp_path, monkeypatch):
    s = sample(tmp_path, monkeypatch, with_runs((200, "A", 6), (401, "C", 5), (520, "G", 8)))
    annotate_excluded(s, [trf_row(100, 145, "ACGTTGCAT"), trf_row(130, 160, "ACGTTGCATA"),
                          trf_row(515, 540, "ACGTTGCAT")])
    s.run("repetitive_dna", f"excluded_summary_{SPECIES}_{SAMPLE}_diploid")
    # Regions are merged intervals, and the share is of the contigs' 1000 bp.
    assert [(row["class"], row["regions"], row["bp"], row["share"]) for row in s.table("excluded_summary")] == [
        ("homopolymer", "2", "14", "0.014"),   # [200, 206), [520, 528)
        ("long_period", "2", "85", "0.085"),   # [100, 160), [515, 540)
        ("long_period_timed_out", "0", "0", "0"),
        ("excluded", "3", "91", "0.091"),      # their union
    ]


def test_a_window_trf_cant_finish_in_time_is_excluded_whole(tmp_path, monkeypatch):
    # In a satellite array TRF can run for days: the window's core, here the
    # whole contig, is excluded whatever TRF had found, and recorded.
    s = sample(tmp_path, monkeypatch, with_runs((200, "A", 6)))
    annotate_excluded(s, [trf_row(100, 145, "ACGTTGCAT")], timed_out=True)
    assert s.regions("long_periods", contig=CONTIG) == [(0, 1000)]
    assert s.regions("long_periods_timed_out", contig=CONTIG) == [(0, 1000)]
    assert s.regions("excluded", contig=CONTIG) == [(0, 1000)]
    s.run("repetitive_dna", f"excluded_summary_{SPECIES}_{SAMPLE}_diploid")
    rows = {row["class"]: (row["regions"], row["bp"], row["share"]) for row in s.table("excluded_summary")}
    assert rows["long_period_timed_out"] == ("1", "1000", "1")


def test_the_tracts_workflow_4_drops_are_summarised_per_sample_and_model(tmp_path, monkeypatch):
    s = sample(tmp_path, monkeypatch, with_runs((200, "A", 6)))
    annotate_excluded(s)
    annotate_repeats(s, [(150, 201, "AC", 2), (206, 240, "AG", 2), (380, 420, "AC", 2)])
    measure_tracts(s)
    assert [(row["model"], row["tracts"], row["bp"], row["tracts_measured"], row["bp_measured"])
            for row in s.table("repetitive_tracts_summary")] == [(model, "3", "125", "2", "74") for model in MODELS]


# -- tracts dropped by overlap resolution (R9) ---------------------------------
#
# Resolution keeps one tract per base per model: the AC tract, of the smaller
# period, and drops the AAG tract that overlaps it whole, bases 140-190 too.
OVERLAPPING = [(100, 140, "AC", 2), (130, 190, "AAG", 3)]


@pytest.mark.parametrize("model", MODELS)
def test_workflow_3_masks_every_base_of_a_tract_that_overlap_resolution_drops(tmp_path, monkeypatch, model):
    s = sample(tmp_path, monkeypatch)
    annotate_excluded(s)
    annotate_repeats(s, OVERLAPPING, models=[model])
    assert build_mask(s) == [(100, 190)]


def test_workflow_4_still_measures_one_tract_per_base_per_model(tmp_path, monkeypatch):
    s = sample(tmp_path, monkeypatch)
    annotate_excluded(s)
    annotate_repeats(s, OVERLAPPING)
    assert measure_tracts(s) == {model: [(100, 140)] for model in MODELS}


def test_the_bases_of_dropped_tracts_are_summarised_per_sample_and_model(tmp_path, monkeypatch):
    s = sample(tmp_path, monkeypatch)
    annotate_repeats(s, OVERLAPPING + [(300, 330, "AC", 2)], models=["M1"])
    for model in MODELS:
        s.run("repetitive_dna", f"cat_{model}_{SPECIES}_{SAMPLE}_diploid")
    # Before resolution 100-190 and 300-330 are in a tract; 140-190 is in no
    # tract after.
    assert [(row["model"], row["bp_before"], row["bp_kept"], row["bp_dropped"], row["share_dropped"])
            for model in MODELS for row in s.table("repeats_dropped_summary", model=model)] == [
        ("M0", "0", "0", "0", "0"),
        ("M1", "120", "70", "50", "0.416667"),
        ("M2", "0", "0", "0", "0"),
    ]


# -- unresolved haplotype regions (decisions Q19, Q32) -------------------------


def test_the_share_of_the_genome_the_unresolved_regions_mask_is_written_per_sample(tmp_path, monkeypatch):
    s = sample(tmp_path, monkeypatch)
    write_unresolved(s, [(100, 300), (500, 550)])
    s.run("precompute", f"unresolved_summary_{SPECIES}_{SAMPLE}_diploid")
    assert [(row["regions"], row["bp"], row["contigs_bp"], row["share"])
            for row in s.table("unresolved_summary")] == [("2", "250", "1000", "0.25")]


def test_workflow_4_drops_every_tract_that_overlaps_an_unresolved_region_whole(tmp_path, monkeypatch):
    s = sample(tmp_path, monkeypatch)
    annotate_excluded(s)
    annotate_repeats(s, [(100, 140, "AC", 2), (290, 320, "AG", 2), (400, 440, "AC", 2)])
    # The second region only touches the third tract.
    assert measure_tracts(s, unresolved=[(200, 300), (440, 460)]) == \
        {model: [(100, 140), (400, 440)] for model in MODELS}


def test_workflow_3_still_masks_the_unresolved_regions(tmp_path, monkeypatch):
    s = sample(tmp_path, monkeypatch)
    annotate_excluded(s)
    annotate_repeats(s)
    assert build_mask(s, unresolved=[(600, 700)]) == [(600, 700)]


def test_trf_runs_on_each_sub_window_in_turn_and_only_a_slow_one_is_excluded(tmp_path, monkeypatch):
    s = sample(tmp_path, monkeypatch)
    with open(f"{s.root}/config.yaml") as handle:
        config = yaml.safe_load(handle)
    config["repeats"].update({"long_period_window_size": 400, "long_period_overlap": 10})
    with open(f"{s.root}/config.yaml", "w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    # samtools records the stretch it fetches; the second TRF run times out.
    s.stand_in("samtools", f'#!/bin/sh\necho "$3" >> {s.root}/fetched.txt\necho ">window"\n')
    second_times_out = f'#!/bin/sh\nn=$(cat {s.root}/runs 2>/dev/null || echo 0)\nn=$((n + 1))\n' \
                       f'echo $n > {s.root}/runs\nif [ $n -eq 2 ]; then exit 124; fi\nshift\nexec "$@"\n'
    with open(f"{s.root}/trf.txt", "w") as handle:
        handle.write("@window\n")
    s.stand_in("trf", f'#!/bin/sh\ncat {s.root}/trf.txt\n')
    s.stand_in("timeout", second_times_out)
    s.run("repetitive_dna", f"long_period_trf_{CONTIG}_0_{SPECIES}_{SAMPLE}_diploid")
    s.run("repetitive_dna", f"excluded_{CONTIG}_{SPECIES}_{SAMPLE}_diploid")
    with open(f"{s.root}/fetched.txt") as handle:
        assert handle.read().split() == [f"{CONTIG}:1-410", f"{CONTIG}:391-810", f"{CONTIG}:791-1000"]
    assert s.regions("long_periods_timed_out", contig=CONTIG) == [(400, 800)]
    assert s.regions("excluded", contig=CONTIG) == [(400, 800)]
