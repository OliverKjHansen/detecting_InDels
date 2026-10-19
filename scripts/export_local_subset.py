"""
Copy windows of samples' aligned reads from the cluster, with everything
workflows 3 and 4 need to run on them, so the pipeline can be replayed locally.
The mirror is for software testing only: it shows whether the code does on a
laptop what it does on the cluster. No data conclusions are drawn from it.

Run it on the cluster, from the repository root, once workflows 1-4 have run
for the samples. One run exports several samples, for example a Revio sperm
sample, a tissue-labelled one, a methylation one and a Sequel II one:

    python scripts/export_local_subset.py human/chk:h1tg000001l:0-10000000 \\
        chimp/yoran human/ob006 chimp/carl --window-size 1000000

Each sample is species/sample, optionally followed by :contig:start-end, a
window that is 0-based and half-open like the scan's windows. A sample given
without one gets the window of --window-size bp holding the most unique calls
in its results/. A window that is exactly one scan job (it starts at a
multiple of 10 Mb and is 10 Mb long, or runs to the contig's end) is compared
with the cluster call by call and opportunity by opportunity. A smaller one
exports and replays faster but is compared call by call only, because the
cluster counts opportunities per scan job.

Every alignment that overlaps a window is kept, whatever its flags, because
that is what the scans see: the coverage test compares each position's depth
with the sample's mean, so the reads cannot be thinned.

A table missing on the cluster stops the export at that sample; the samples
before it are exported. The final opportunity tables and filter reports are
not compared, because the cluster writes them per sample, not per window.

It needs only pysam and pyyaml. A 10 Mb window of a 30x sample is roughly
0.5 GB without the per-base kinetics tags, which are dropped unless
--keep-kinetics. Each window gets its own folder,
mirror/{species}/{sample}/{contig}_{start}_{end}/. Copy it to the same place
in a local checkout:

    rsync -av CLUSTER:REPO/mirror/ mirror/

and run its run_local.sh there. The folder has the cluster's layout, and its
config.yaml points into it with relative paths:

  config.yaml        the cluster's config for this one sample
  fasta/             the window's contig; the sequence outside the window and a
                     100 kb margin is N, so coordinates match the cluster
  bam/               every alignment overlapping the window, indexed
  output/precomute/  mean coverage, unresolved regions, tissue labels
  output/repetitive_dna/
                     the M0, M1 and M2 tracts near the window
  output/nonrepetitive_indels/, output/repetitive_indels/
                     workflow 3's mask near the window, and both workflows'
                     haploid alignment records (locations.tsv) of the window's
                     reads, so the haploid check runs without remapping
  cluster/           the cluster's tables for the window, per workflow (and per
                     model for workflow 4): the scan's calls and
                     opportunities, the calls after the haploid check, and the
                     final calls
  run_local.sh       runs scripts/replay_window.py on the window: reruns
                     workflows 3 and 4 and compares each table with the
                     cluster's (scripts/compare_tables.py)
"""
import argparse
import os
import shutil
import sys
import time
from collections import Counter
from functools import partial

import pysam
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "workflows"))
from pipeline import Pipeline  # noqa: E402

KINETICS_TAGS = ("fi", "ri", "fp", "rp", "ip", "pw")
MARGIN = 100_000              # real sequence kept on each side of the window

STARTED = time.time()

# run_local.sh: the replay itself is built from the current code by
# replay_window.py, so a window outlives changes to the scripts' arguments.
REPLAY_SCRIPT = """#!/usr/bin/env bash
# Replays workflows 3 and 4 on this window with the scripts of the checkout
# REPO, by default the one this folder sits in (mirror/{species}/{sample}/
# under the repository root), and compares each table with the cluster's.
# Exits 1 if any table differs. For software testing only. PYTHON is the
# Python to run them with, one with pysam and pandas.
set -euo pipefail
cd "$(dirname "$0")"
REPO="${REPO:-$(cd ../../../.. && pwd)}"
exec "${PYTHON:-python}" "$REPO/scripts/replay_window.py" .
"""


def log(message):
    print(f"[{time.time() - STARTED:5.0f} s] {message}", file=sys.stderr, flush=True)


def parse_sample(spec):
    """'species/sample[:contig:start-end]' -> (species, sample, window or None)."""
    name, _, window = spec.partition(":")
    species, sample = name.split("/")
    if not window:
        return species, sample, None
    contig, span = window.rsplit(":", 1)
    start, end = (int(x.replace(",", "")) for x in span.split("-"))
    return species, sample, (contig, start, end)


def busiest_window(final_calls, contigs, size):
    """The window of `size` bp holding the most unique calls in the final results."""
    counts = Counter()
    with open(final_calls) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        contig, position, unique = (header.index(c) for c in ("contig", "position", "unique"))
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[unique] == "unique" and fields[contig] in contigs:
                counts[(fields[contig], (int(fields[position]) - 1) // size)] += 1
    if not counts:
        sys.exit(f"No unique calls in {final_calls}; give a window.")
    (contig, index), n = counts.most_common(1)[0]
    log(f"window with the most unique calls: {n} calls")
    return contig, index * size, min((index + 1) * size, contigs[contig])


def copy_overlapping(source, target, contig, start, end, header):
    """Copy the lines of a BED-like file on `contig` that overlap [start, end)."""
    if not os.path.exists(source):
        return 0
    kept = 0
    os.makedirs(os.path.dirname(target), exist_ok=True)
    with open(source) as handle, open(target, "w") as out:
        if header:
            out.write(handle.readline())
        for line in handle:
            fields = line.split("\t", 3)
            if fields[0] == contig and int(fields[1]) < end and int(fields[2]) > start:
                out.write(line)
                kept += 1
    return kept


def copy_calls_in_window(source, target, contig, start, end):
    """Copy the rows of a call table whose anchor lies in the window."""
    if not os.path.exists(source):
        return None
    kept = 0
    os.makedirs(os.path.dirname(target), exist_ok=True)
    with open(source) as handle, open(target, "w") as out:
        first = handle.readline()
        out.write(first)
        header = first.rstrip("\n").lstrip("#").split("\t")
        contig_column, position = header.index("contig"), header.index("position")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[contig_column] == contig and start <= int(fields[position]) - 1 < end:
                out.write(line)
                kept += 1
    return kept


def copy_rows_of_reads(source, target, names, header=False):
    """Copy the rows of a table keyed by read name in its first column, for the window's reads."""
    os.makedirs(os.path.dirname(target), exist_ok=True)
    with open(source) as handle, open(target, "w") as out:
        if header:
            out.write(handle.readline())
        out.writelines(line for line in handle if line.split("\t", 1)[0] in names)


def copy_cluster_tables(source, species, sample, workflow, compare, contig, start, end, is_scan_job, **fields):
    """
    The cluster's tables for the window, from one workflow ("nonrepetitive" or
    "repetitive", and for workflow 4 one model): the scan's calls and, for a
    whole scan job, opportunities; the calls after the haploid check; the
    final calls.
    """
    path = partial(source.path, species=species, sample=sample, **fields)
    os.makedirs(compare, exist_ok=True)
    if is_scan_job:
        for table in ("calls", "opportunities"):
            shutil.copy(path(f"{workflow}_window_{table}", contig=contig, start=start, end=end), f"{compare}/{table}.tsv")
    else:
        copy_calls_in_window(path(f"{workflow}_calls", contig=contig), f"{compare}/calls.tsv", contig, start, end)
    copy_calls_in_window(path(f"{workflow}_merged_calls"), f"{compare}/merged_indels.tsv", contig, start, end)
    copy_calls_in_window(path(f"{workflow}_results_calls"), f"{compare}/indels.tsv", contig, start, end)


def mirror_config(config, species, sample, root):
    """The cluster's config, cut down to the one sample, with paths into the mirror folder."""
    def only_this(table):
        listed = (config.get(table) or {}).get(species) or []
        return {species: [sample]} if sample in listed else {}

    mirrored = dict(config)
    mirrored.update({
        "samples": {species: [sample]},
        "tissue_samples": only_this("tissue_samples"),
        "methylation_samples": only_this("methylation_samples"),
        "sperm_samples": only_this("sperm_samples"),
        "sample_age": {species: {sample: config["sample_age"][species][sample]}}
                      if sample in ((config.get("sample_age") or {}).get(species) or {}) else {},
        "alignment_overrides": {},
        # Relative, so the folder works wherever it is copied to.
        "paths": {name: name for name in config["paths"]},
    })
    with open(f"{root}/config.yaml", "w") as out:
        yaml.safe_dump(mirrored, out, sort_keys=False)
    return Pipeline(f"{root}/config.yaml", base=root)


def export_window(source, species, sample, contig, start, end, out, keep_kinetics):
    """Export one window of one sample to {out}/{species}/{sample}/{contig}_{start}_{end}/."""
    contigs = source.long_contigs(species, sample)
    if contig not in contigs or not 0 <= start < end <= contigs[contig]:
        sys.exit(f"{species}/{sample}: {contig}:{start}-{end} is not inside a contig longer than min_contig_length")
    scan_window = source.config["scan"]["window_size"]
    is_scan_job = start % scan_window == 0 and end == min(start + scan_window, contigs[contig])
    root = os.path.abspath(os.path.join(out, species, sample, f"{contig}_{start}_{end}"))
    os.makedirs(root, exist_ok=True)
    target = mirror_config(source.config, species, sample, root)
    has_workflow4 = os.path.exists(source.path("repetitive_haploid_locations", species, sample))
    log(f"{species}/{sample} {contig}:{start}-{end}"
        f" ({'one scan job' if is_scan_job else 'not a scan job: calls compared, opportunities not'})")
    if not has_workflow4:
        log("  workflow 4 has not run for this sample: exporting workflow 3 only")

    def both(name, **fields):
        return source.path(name, species, sample, **fields), target.path(name, species, sample, **fields)

    # The reads: every alignment overlapping the window, sorted as fetched.
    log("  alignments")
    alignment = target.alignment(species, sample)
    os.makedirs(os.path.dirname(alignment), exist_ok=True)
    names, records = set(), 0
    with pysam.AlignmentFile(source.alignment(species, sample), "rb") as bam, \
            pysam.AlignmentFile(alignment, "wb", template=bam) as copy:
        for read in bam.fetch(contig, start, end):
            if not keep_kinetics:
                for tag in KINETICS_TAGS:
                    if read.has_tag(tag):
                        read.set_tag(tag, None)
            copy.write(read)
            names.add(read.query_name)
            records += 1
    pysam.index(alignment)
    log(f"  {records:,} alignment records of {len(names):,} reads")

    # The contig at full length, so coordinates match; N away from the window.
    log("  assembly")
    source_assembly, assembly = both("assembly")
    os.makedirs(os.path.dirname(assembly), exist_ok=True)
    length = contigs[contig]
    left, right = max(0, start - MARGIN), min(length, end + MARGIN)
    with pysam.FastaFile(source_assembly) as fasta:
        sequence = fasta.fetch(contig, left, right)
    with open(assembly, "w") as handle:
        handle.write(f">{contig}\n")
        full_lines, rest = divmod(left, 60)
        handle.writelines("N" * 60 + "\n" for _ in range(full_lines))
        body = "N" * rest + sequence + "N" * (length - right)
        handle.writelines(body[i:i + 60] + "\n" for i in range(0, len(body), 60))
    pysam.faidx(assembly)

    # Workflows 1 and 2: coverage, unresolved regions, tissue labels, repeat tracts.
    log("  precomputed inputs and repeat tracts")
    os.makedirs(target.folder("coverage", species, sample), exist_ok=True)
    shutil.copy(*both("coverage"))
    with open(target.path("long_contigs", species, sample), "w") as handle:
        handle.write(f"{contig}\t0\t{length}\n")
    copy_overlapping(*both("unresolved", contig=contig), contig, left, right, header=False)
    tissue = "NA"
    if source.is_listed("tissue_samples", species, sample):
        source_tissue, tissue = both("tissue")
        copy_rows_of_reads(source_tissue, tissue, names, header=True)
    source_classification, classification = both("read_classification")
    if source.is_listed("methylation_samples", species, sample) and os.path.exists(source_classification):
        copy_rows_of_reads(source_classification, classification, names, header=True)
    models = source.config["repeats"]["models"]
    for model in models:
        copy_overlapping(*both("tandem_repeats", contig=contig, model=model), contig, left, right, header=True)

    # Workflows 3 and 4: the mask, and the haploid alignments of the window's reads.
    log("  mask and haploid alignment records")
    copy_overlapping(*both("nonrepetitive_mask", contig=contig), contig, left, right, header=False)
    for workflow in ["nonrepetitive"] + (["repetitive"] if has_workflow4 else []):
        copy_rows_of_reads(*both(f"{workflow}_haploid_locations"), names)

    # The cluster's tables for the window, to compare a local run with.
    log("  the cluster's tables for the window")
    copy_cluster_tables(source, species, sample, "nonrepetitive", f"{root}/cluster/nonrepetitive_indels",
                        contig, start, end, is_scan_job)
    if has_workflow4:
        for model in models:
            copy_cluster_tables(source, species, sample, "repetitive", f"{root}/cluster/repetitive_indels/{model}",
                                contig, start, end, is_scan_job, model=model)

    with open(f"{root}/run_local.sh", "w") as handle:
        handle.write(REPLAY_SCRIPT)
    os.chmod(f"{root}/run_local.sh", 0o755)

    sizes = {os.path.relpath(os.path.join(d, f), root): os.path.getsize(os.path.join(d, f))
             for d, _, files in os.walk(root) for f in files}
    with open(f"{root}/README.txt", "w") as handle:
        handle.write(f"Made by scripts/export_local_subset.py on {time.strftime('%Y-%m-%d %H:%M')}\n"
                     f"command: {' '.join(sys.argv)}\n"
                     f"window: {contig}:{start}-{end} ({'one scan job' if is_scan_job else 'not a scan job'})\n"
                     f"BAM: {source.alignment(species, sample)}\n"
                     f"alignment records: {records:,} of {len(names):,} reads\n"
                     f"kinetics kept: {keep_kinetics}\n"
                     f"workflow 4 exported: {has_workflow4}\n"
                     "For software testing only: no data conclusions are drawn from the mirror.\n\n"
                     "files (MB):\n")
        handle.writelines(f"  {name}: {size / 1e6:.1f}\n" for name, size in sorted(sizes.items()))
    log(f"  done: {sum(sizes.values()) / 1e6:.0f} MB in {root}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("samples", nargs="+", metavar="species/sample[:contig:start-end]",
                        help="the samples to export, each with its window or none")
    parser.add_argument("--window-size", type=int,
                        help="size of the window chosen for a sample given without one (default: one scan job)")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--out", default="mirror",
                        help="each window goes to {out}/{species}/{sample}/{contig}_{start}_{end}/")
    parser.add_argument("--keep-kinetics", action="store_true", help="keep the per-base kinetics tags")
    args = parser.parse_args()

    source = Pipeline(args.config)
    window_size = args.window_size or source.config["scan"]["window_size"]
    for spec in args.samples:
        species, sample, window = parse_sample(spec)
        if window is None:
            window = busiest_window(source.path("nonrepetitive_results_calls", species, sample),
                                    source.long_contigs(species, sample), window_size)
        export_window(source, species, sample, *window, args.out, args.keep_kinetics)


if __name__ == "__main__":
    main()
