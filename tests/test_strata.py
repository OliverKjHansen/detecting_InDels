"""
Calls and opportunities are split by one set of keys, scripts/strata.py, so a
key added there is a column of both tables, in both scans.
"""

import pytest

import strata
from scans import nonrepetitive_scan, rows_of, tandem_repeat_scan


@pytest.fixture
def dummy_keys(monkeypatch):
    """
    One more read key and one more site key: the read's name, and the lowest
    quality tested with the exact distance to the read end, which splits the
    opportunities finer than any key already there.
    """
    keys = strata.keys

    def more_keys(**parameters):
        read_keys, site_keys = keys(**parameters)
        return (read_keys + [("dummy_read", lambda read: read.name)],
                site_keys + [("dummy_site", lambda site: f"{min(site.qualities)}/{site.read_end_distance}")])

    monkeypatch.setattr(strata, "keys", more_keys)


def tables(fixture_paths):
    calls, opportunities, _ = nonrepetitive_scan(fixture_paths)
    tr_calls, tr_opportunities, _ = tandem_repeat_scan({"M1": fixture_paths["tr_bed"]}, fixture_paths)
    return {"workflow 3 calls": rows_of(calls), "workflow 3 opportunities": rows_of(opportunities),
            "workflow 4 calls": rows_of(tr_calls["M1"]), "workflow 4 opportunities": rows_of(tr_opportunities["M1"])}


def test_a_key_added_once_is_a_column_of_both_tables_in_both_scans(dummy_keys, fixture_paths):
    for name, rows in tables(fixture_paths).items():
        assert rows, name
        assert all(row["dummy_read"].startswith("read_") for row in rows), name
        # The site key saw the bases kmerQ_pass tests.
        assert all((int(row["dummy_site"].split("/")[0]) >= 93) == (row["kmerQ_pass"] == "pass")
                   for row in rows), name
        if "calls" in name:   # a call carries its own read's value
            assert all(row["dummy_read"] == row["readID"] for row in rows), name


def test_a_new_key_splits_the_opportunities_without_changing_their_total(dummy_keys, fixture_paths, monkeypatch):
    split = tables(fixture_paths)
    monkeypatch.undo()
    plain = tables(fixture_paths)
    for name in ("workflow 3 opportunities", "workflow 4 opportunities"):
        assert sum(int(row["count"]) for row in split[name]) == sum(int(row["count"]) for row in plain[name])
        assert len(split[name]) > len(plain[name])
