"""
Workflow 3's own keys and call columns (ticket 49): kmerQ by position, the
subreads disagreeing at the junction, whether a deletion's bases repeat next
to them, whether the base after it is masked, and the calls' details.
"""
import itertools

import pytest

import strata
import utilities
from scans import KEY_PARAMETERS, nonrepetitive_scan, rows_of

KEYS = ["inner_flank_q", "outer_flank_q", "junction_disagreement", "repeats_next", "base_after_masked"]


def values(qualities, junction=None, repeats=None, after_masked=None):
    site = strata.Site(qualities, None, 1000, None, junction, repeats, after_masked)
    return {name: key(site) for name, key in strata.nonrepetitive_site_keys(**KEY_PARAMETERS)}


def test_kmerq_by_position_bins_the_inner_and_the_outer_flanks():
    # site.qualities is L2 L1 R1 R2
    assert values([60, 93, 93, 59])["inner_flank_q"] == ">=93"
    assert values([60, 93, 93, 59])["outer_flank_q"] == "<60"
    assert values([93, 92, 93, 70])["inner_flank_q"] == "60-93"
    assert values([93, 92, 93, 70])["outer_flank_q"] == "60-93"


def test_the_proposed_kmerq_rule_is_a_choice_of_bins():
    """L1 and R1 at Q93, L2 and R2 at Q60 or more: inner >=93, outer 60-93 or >=93."""
    for q in itertools.product([0, 59, 60, 92, 93], repeat=4):
        v = values(list(q))
        proposed = v["inner_flank_q"] == ">=93" and v["outer_flank_q"] in ("60-93", ">=93")
        assert proposed == (min(q[1], q[2]) >= 93 and min(q[0], q[3]) >= 60), q


def test_junction_disagreement_is_any_disagreeing_subread_at_l1_or_r1():
    assert values([93] * 4, junction=(0, 0))["junction_disagreement"] == "no"
    assert values([93] * 4, junction=(0, 2))["junction_disagreement"] == "yes"
    assert values([93] * 4)["junction_disagreement"] == "NA"


def test_the_yes_no_keys_are_na_where_they_dont_apply():
    assert values([93] * 4)["repeats_next"] == "NA"
    assert values([93] * 4, repeats=True)["repeats_next"] == "yes"
    assert values([93] * 4, after_masked=False)["base_after_masked"] == "no"
    assert values([93] * 4, after_masked=True)["base_after_masked"] == "yes"


def test_only_workflow_3_has_its_keys():
    assert not set(KEYS) & set(strata.Strata(**KEY_PARAMETERS).columns)
    assert set(KEYS) <= set(strata.Strata(nonrepetitive=True, **KEY_PARAMETERS).columns)


def test_workflow_3_leaves_out_keys_and_bins_its_own_way_and_workflow_4_doesnt():
    import yaml
    from conftest import REPO_ROOT
    with open(f"{REPO_ROOT}/config.yaml") as handle:
        parameters = yaml.safe_load(handle)["strata"]
    read = strata.Read("r", "NA", "NA", 15000, False, 25, 15, 0.9995, 6, 3)
    site = strata.Site([93] * 4, None, 2000, 30, (0, 0), None, False)
    for nonrepetitive, expected in ((True, {"np_bin": "12-17", "mapQ_bin": "20-60", "nearby_event_bin": "20-50"}),
                                    (False, {"np_bin": "15-18", "mapQ_bin": "20-30", "nearby_event_bin": "25-35"})):
        keys = strata.Strata(nonrepetitive=nonrepetitive, **parameters)
        values = dict(zip(keys.columns, keys.of_read(read) + keys.of_site(site)))
        assert {name: values[name] for name in expected} == expected
        assert ("strand_passes_bin" in keys.columns) != nonrepetitive
    np_bin = strata.Strata(nonrepetitive=True, **parameters)
    labels = {dict(zip(np_bin.columns, np_bin.of_read(read._replace(passes=p))))["np_bin"] for p in (3, 7, 12, 17, 40)}
    assert labels == {"<7", "7-12", "12-17", ">=17"}


@pytest.mark.parametrize("sequence, start, k, expected", [
    ("GGACACTT", 2, 2, True),    # AC, then AC
    ("ACAC", 2, 2, True),        # AC before
    ("GTACGT", 2, 2, False),
    ("ACGTT", 0, 2, False),      # nothing before
    ("CACAG", 1, 3, False),      # ACA: too few bases on either side
    ("TTAC", 2, 2, False),       # nothing after
    ("ATGATGC", 0, 3, True),
])
def test_repeats_next(sequence, start, k, expected):
    assert utilities.repeats_next(sequence, start, k) == expected


@pytest.fixture
def mask(tmp_path):
    bed = tmp_path / "mask.bed"
    bed.write_text("c\t10\t20\nc\t50\t60\n")   # bases 10..19 and 50..59
    return utilities.Mask(str(bed))


@pytest.mark.parametrize("first, last, limit, expected", [
    (18, 22, 100, 0),     # a masked base inside
    (20, 23, 100, 1),     # 19, right before
    (46, 49, 100, 1),     # 50, right after
    (25, 28, 100, 6),     # 19 is 6 before 25
    (30, 32, 5, None),    # 11 and 18 away: beyond the limit
])
def test_the_distance_to_the_mask(mask, first, last, limit, expected):
    assert mask.distance(first, last, limit) == expected


@pytest.fixture(scope="module")
def scanned(fixture_paths):
    calls, opportunities, _ = nonrepetitive_scan(fixture_paths)
    return rows_of(calls), rows_of(opportunities)


def test_every_call_and_opportunity_carries_the_keys_and_every_call_the_details(scanned):
    calls, opportunities = scanned
    assert calls and opportunities
    assert all(key in row for row in calls + opportunities for key in KEYS)
    assert all(column in row for row in calls for column in utilities.CALL_DETAIL_COLUMNS)


def test_kmerq_pass_is_both_flank_bins_at_q93(scanned):
    calls, opportunities = scanned
    for row in calls + opportunities:
        assert (row["kmerQ_pass"] == "pass") == (row["inner_flank_q"] == row["outer_flank_q"] == ">=93")


def test_a_calls_flank_bins_and_junction_come_from_its_detail_columns(scanned):
    flank_bin = strata.binned(KEY_PARAMETERS["flank_quality_bins"])
    for row in scanned[0]:
        assert row["inner_flank_q"] == flank_bin(min(int(row["qual_L1"]), int(row["qual_R1"])))
        assert row["outer_flank_q"] == flank_bin(min(int(row["qual_L2"]), int(row["qual_R2"])))
        assert int(row["kmerQ"]) == min(int(row[f"qual_{b}"]) for b in ("L2", "L1", "R1", "R2"))
        if row["sx_L1"] == "NA":
            assert row["junction_disagreement"] == "NA"
        else:
            assert (row["junction_disagreement"] == "yes") == (max(int(row["sx_L1"]), int(row["sx_R1"])) > 0)


def test_insertions_have_neither_a_masked_base_after_nor_a_repeat_key(scanned):
    calls, opportunities = scanned
    insertions = [row for row in calls + opportunities if row["indel_type"] == "insertion"]
    assert insertions
    assert all(row["base_after_masked"] == "no" and row["repeats_next"] == "NA" for row in insertions)


def test_an_insertion_calls_inserted_bases_pass_at_the_minimum_quality(scanned):
    for row in scanned[0]:
        if row["indel_type"] == "insertion":
            expected = "pass" if int(row["inserted_min_quality"]) >= KEY_PARAMETERS["min_inserted_quality"] else "fail"
            assert row["inserted_bases_pass"] == expected
        else:
            assert row["inserted_min_quality"] == row["inserted_bases_pass"] == "NA"


def test_a_deletion_call_repeats_next_by_the_assembly(scanned):
    for row in scanned[0]:
        if row["indel_type"] == "deletion":
            # flank_sequence: 25 bases, the anchor, then 25; the deleted bases start at 26
            flank, k = row["flank_sequence"].upper(), len(row["indel_string"])
            assert (row["repeats_next"] == "yes") == utilities.repeats_next(flank, 26, k)


def test_the_base_after_a_deletion_masked_is_recorded_not_filtered(fixture_paths, tmp_path):
    """Masking the base right after a deletion call keeps the call, and records it, on both tables."""
    calls, _, _ = nonrepetitive_scan(fixture_paths)
    deletion = next(row for row in rows_of(calls) if row["indel_type"] == "deletion")
    anchor, k = int(deletion["position"]) - 1, len(deletion["indel_string"])   # position is 1-based
    with open(fixture_paths["filters"]) as handle:
        bed = handle.read()
    masked = tmp_path / "mask.bed"
    contig = bed.split("\t", 1)[0] if bed else deletion["contig"]
    masked.write_text(bed + f"{contig}\t{anchor + k + 1}\t{anchor + k + 2}\n")
    calls, opportunities, _ = nonrepetitive_scan({**fixture_paths, "filters": str(masked)})
    same = [row for row in rows_of(calls)
            if (row["readID"], row["position"], row["indel_string"]) ==
               (deletion["readID"], deletion["position"], deletion["indel_string"])]
    assert len(same) == 1
    assert same[0]["base_after_masked"] == "yes" and same[0]["mask_distance"] == "0"
    deletions = [row for row in rows_of(opportunities) if row["indel_type"] == "deletion"]
    assert any(row["base_after_masked"] == "yes" for row in deletions)
    for row in rows_of(calls):
        if row["indel_type"] == "deletion":
            assert (row["base_after_masked"] == "yes") == (row["mask_distance"] == "0")
