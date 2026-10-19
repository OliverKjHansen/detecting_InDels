"""Strip outputs and execution counts from a Jupyter notebook.

Reads a notebook on stdin and writes it to stdout without outputs or execution
counts. Used as a git clean filter (see .gitattributes and the README), so
committed notebooks carry no outputs while the copies on disk keep them.

Standard library only, so the filter runs under any python3 without the
project environment.
"""

import json
import sys


def strip(notebook):
    """Empty every code cell's outputs and execution count, in place."""
    for cell in notebook.get("cells", []):
        if cell.get("cell_type") == "code":
            cell["outputs"] = []
            cell["execution_count"] = None
    return notebook


def main():
    """Notebook in on stdin, stripped notebook out on stdout."""
    notebook = json.load(sys.stdin)
    # The layout nbformat itself writes, so a re-saved notebook diffs cleanly.
    sys.stdout.write(json.dumps(strip(notebook), indent=1, sort_keys=True,
                                ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
