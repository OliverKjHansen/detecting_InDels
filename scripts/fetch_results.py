#!/usr/bin/env python3
"""
Pull finished results down from the cluster for local plotting in R.

Only `output/results/` comes down — the filtered InDel tables, the opportunity
tables, the filter reports and the diagnostic PDFs. `output/work/` holds the
per-contig intermediates, the FASTQs and the BAMs; those stay on the cluster.

    python scripts/fetch_results.py                     # everything
    python scripts/fetch_results.py --dry-run           # show what would transfer
    python scripts/fetch_results.py marmoset            # one species
    python scripts/fetch_results.py marmoset/lambiek    # one sample
    python scripts/fetch_results.py --delete            # mirror exactly

The local tree mirrors the remote one, so `results/{species}/{sample}/...` here
is the same path as `output/results/{species}/{sample}/...` there.

Host, remote working directory and local destination come from the `remote`
block in config.yaml.
"""

import argparse
import os
import subprocess
import sys

import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# The subtree of output_root that holds results. Everything else is work/.
RESULTS_SUBDIR = "results"


def load_config():
    with open(os.path.join(REPO_ROOT, "config.yaml")) as handle:
        return yaml.safe_load(handle)


def build_command(host, workdir, output_root, local_dir, subset="",
                  dry_run=False, delete=False):
    """Assemble the rsync invocation."""
    base = f"{workdir.rstrip('/')}/{output_root}/{RESULTS_SUBDIR}"
    if subset:
        remote = f"{host}:{base}/{subset.strip('/')}/"
        local = os.path.join(local_dir, subset.strip("/")) + os.sep
    else:
        remote = f"{host}:{base}/"
        local = local_dir.rstrip("/") + os.sep

    command = ["rsync", "-avz", "--partial", "--human-readable"]
    if dry_run:
        command.append("--dry-run")
    if delete:
        command.append("--delete")
    command += [remote, local]
    return command


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("subset", nargs="?", default="",
                        help="limit to a species ('marmoset') or a sample "
                             "('marmoset/lambiek'); default is everything")
    parser.add_argument("--dry-run", action="store_true",
                        help="show what would be transferred, copy nothing")
    parser.add_argument("--delete", action="store_true",
                        help="remove local files that no longer exist remotely")
    args = parser.parse_args()

    cfg = load_config()
    remote = cfg.get("remote", {})
    for key in ("host", "workdir", "local_dir"):
        if not remote.get(key):
            sys.exit(f"remote.{key} is not set in config.yaml; cannot fetch results")

    output_root = cfg["paths"]["output_root"]
    local_dir = remote["local_dir"]

    command = build_command(
        remote["host"], remote["workdir"], output_root, local_dir, args.subset,
        dry_run=args.dry_run, delete=args.delete,
    )

    destination = os.path.join(local_dir, args.subset) if args.subset else local_dir
    if not args.dry_run:
        os.makedirs(destination, exist_ok=True)

    print(f"  $ {' '.join(command)}", file=sys.stderr)
    result = subprocess.run(command)
    if result.returncode != 0:
        sys.exit(
            f"rsync failed (exit {result.returncode}). Check that you can reach "
            f"{remote['host']} and that the pipeline has produced results there."
        )

    if not args.dry_run:
        print(f"\nDone. Results are in {destination}/ — ready to load in R:\n"
              f'  files <- list.files("{local_dir}", "^indels\\\\.tsv$", '
              f"recursive = TRUE, full.names = TRUE)",
              file=sys.stderr)


if __name__ == "__main__":
    main()
