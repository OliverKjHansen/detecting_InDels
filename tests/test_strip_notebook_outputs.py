"""Committed notebooks carry no outputs; local copies keep theirs.

strip_notebook_outputs.py is a git clean filter: git pipes a notebook through it
on its way into the index, so what gets committed has no outputs or execution
counts, while the file on disk keeps both. Saved plots made the three analysis
notebooks 3-8 MB each, nearly all of it output.
"""

import copy
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
STRIP = REPO / "scripts" / "strip_notebook_outputs.py"

NOTEBOOK = {
    "cells": [
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": ["# Rates outside repeats\n"],
        },
        {
            "cell_type": "code",
            "execution_count": 7,
            "metadata": {"vscode": {"languageId": "r"}},
            "outputs": [
                {"name": "stdout", "output_type": "stream", "text": ["60486\n"]},
                {
                    "data": {"image/png": "iVBORw0KGgo="},
                    "metadata": {},
                    "output_type": "display_data",
                },
            ],
            "source": ["nrow(calls)\n"],
        },
    ],
    "metadata": {"kernelspec": {"display_name": "R", "language": "R", "name": "ir"}},
    "nbformat": 4,
    "nbformat_minor": 5,
}


def run_filter(text):
    """The script as git runs it: notebook text in, notebook text out."""
    result = subprocess.run([sys.executable, str(STRIP)], input=text,
                            capture_output=True, text=True, check=True)
    return result.stdout


def stripped_notebook():
    return json.loads(run_filter(json.dumps(NOTEBOOK)))


def test_outputs_and_execution_counts_are_removed():
    code = stripped_notebook()["cells"][1]

    assert code["outputs"] == []
    assert code["execution_count"] is None
    assert code["source"] == ["nrow(calls)\n"]
    assert code["metadata"] == {"vscode": {"languageId": "r"}}


def test_stripping_twice_changes_nothing():
    """git re-runs the filter on already clean files; a second pass must be a no-op."""
    notebook = copy.deepcopy(NOTEBOOK)
    notebook["cells"][0]["source"] = ["Target: about 1×10⁻⁹ per bp\n"]

    once = run_filter(json.dumps(notebook))

    assert run_filter(once) == once
    assert "1×10⁻⁹" in once


def test_markdown_and_notebook_metadata_are_untouched():
    stripped = stripped_notebook()

    assert stripped["cells"][0] == NOTEBOOK["cells"][0]
    assert stripped["metadata"] == NOTEBOOK["metadata"]
    assert stripped["nbformat"] == 4
    assert stripped["nbformat_minor"] == 5


def git(repo, *args):
    return subprocess.run(["git", "-C", str(repo), *args], capture_output=True,
                          text=True, check=True).stdout


def test_git_stages_notebooks_without_outputs_and_leaves_disk_alone(tmp_path):
    """The repo's own .gitattributes, with the filter set up as the README says."""
    git(tmp_path, "init", "-q")
    (tmp_path / ".gitattributes").write_bytes((REPO / ".gitattributes").read_bytes())
    git(tmp_path, "config", "filter.strip-notebook-outputs.clean",
        f'"{sys.executable}" "{STRIP}"')
    on_disk = json.dumps(NOTEBOOK, indent=1)
    (tmp_path / "plotting").mkdir()
    (tmp_path / "plotting" / "rates.ipynb").write_text(on_disk)

    git(tmp_path, "add", "plotting/rates.ipynb")

    staged = json.loads(git(tmp_path, "show", ":plotting/rates.ipynb"))
    assert staged["cells"][1]["outputs"] == []
    assert staged["cells"][1]["execution_count"] is None
    assert (tmp_path / "plotting" / "rates.ipynb").read_text() == on_disk
    # The outputs on disk don't count as unstaged changes.
    assert git(tmp_path, "diff", "--name-only") == ""
