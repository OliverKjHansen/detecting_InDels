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
  run_local.sh       reruns workflows 3 and 4 on the window and compares each
                     table with the cluster's (scripts/compare_tables.py)
"""
import argparse
import os
import shutil
import sys
import time
from collections import Counter

import pysam
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "workflows"))
from pipeline import Pipeline  # noqa: E402

KINETICS_TAGS = ("fi", "ri", "fp", "rp", "ip", "pw")
SCAN_WINDOW = 10_000_000      # workflows 3 and 4 scan each contig in windows of this size
MARGIN = 100_000              # real sequence kept on each side of the window
MODELS = ["M0", "M1", "M2"]   # the repeat models of workflow 2

# The final filters, as workflows 3 and 4 apply them.
MIN_KMERQ = 93
MIN_READ_END_DISTANCE = 100   # workflow 4 only

STARTED = time.time()


def log(message):
    print(f"[{time.time() - STARTED:5.0f} s] {message}", file=sys.stderr, flush=True)


def workflow_dirs(pipeline, species, sample):
    """Workflows 3 and 4's working and results folders, built as the workflows build them."""
    output, results = pipeline.roots["output"], pipeline.roots["results"]
    return {
        "nonrepetitive": f"{output}/nonrepetitive_indels/{species}/{sample}",
        "nonrepetitive_results": f"{results}/nonrepetitive_indels/{species}/{sample}",
        "repetitive": f"{output}/repetitive_indels/{species}/{sample}",
        "repetitive_results": f"{results}/repetitive_indels/{species}/{sample}",
    }


def parse_sample(spec):
    """'species/sample[:contig:start-end]' -> (species, sample, window or None)."""
    name, _, window = spec.partition(":")
    species, sample = name.split("/")
    if not window:
        return species, sample, None
    contig, span = window.rsplit(":", 1)
    start, end = (int(x.replace(",", "")) for x in span.split("-"))
    return species, sample, (contig, start, end)


def busiest_window(results, contigs, size):
    """The window of `size` bp holding the most unique calls in the final results."""
    counts = Counter()
    with open(f"{results}/indels.tsv") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        contig, position, unique = (header.index(c) for c in ("contig", "position", "unique"))
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[unique] == "unique" and fields[contig] in contigs:
                counts[(fields[contig], (int(fields[position]) - 1) // size)] += 1
    if not counts:
        sys.exit(f"No unique calls in {results}/indels.tsv; give a window.")
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


def copy_cluster_tables(work, results, compare, contig, start, end, is_scan_job):
    """
    The cluster's tables for the window, from one workflow (one model of
    workflow 4): the scan's calls and, for a whole scan job, opportunities;
    the calls after the haploid check; the final calls.
    """
    os.makedirs(compare, exist_ok=True)
    chunk = f"{contig}_{start}_{end}"
    if is_scan_job:
        for table in ("calls", "opportunities"):
            shutil.copy(f"{work}/{table}/chunks/{chunk}.tsv", f"{compare}/{table}.tsv")
    else:
        copy_calls_in_window(f"{work}/calls/{contig}.tsv", f"{compare}/calls.tsv", contig, start, end)
    copy_calls_in_window(f"{work}/merged/indels.tsv", f"{compare}/merged_indels.tsv", contig, start, end)
    copy_calls_in_window(f"{results}/indels.tsv", f"{compare}/indels.tsv", contig, start, end)


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
    is_scan_job = start % SCAN_WINDOW == 0 and end == min(start + SCAN_WINDOW, contigs[contig])
    root = os.path.abspath(os.path.join(out, species, sample, f"{contig}_{start}_{end}"))
    os.makedirs(root, exist_ok=True)
    target = mirror_config(source.config, species, sample, root)
    source_dirs, target_dirs = workflow_dirs(source, species, sample), workflow_dirs(target, species, sample)
    has_workflow4 = os.path.exists(f"{source_dirs['repetitive']}/haploid/locations.tsv")
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
    for model in MODELS:
        copy_overlapping(*both("tandem_repeats", contig=contig, model=model), contig, left, right, header=True)

    # Workflows 3 and 4: the mask, and the haploid alignments of the window's reads.
    log("  mask and haploid alignment records")
    copy_overlapping(f"{source_dirs['nonrepetitive']}/mask/{contig}.bed",
                     f"{target_dirs['nonrepetitive']}/mask/{contig}.bed", contig, left, right, header=False)
    workflows = ["nonrepetitive"] + (["repetitive"] if has_workflow4 else [])
    for workflow in workflows:
        copy_rows_of_reads(f"{source_dirs[workflow]}/haploid/locations.tsv",
                           f"{target_dirs[workflow]}/haploid/locations.tsv", names)

    # The cluster's tables for the window, to compare a local run with.
    log("  the cluster's tables for the window")
    copy_cluster_tables(source_dirs["nonrepetitive"], source_dirs["nonrepetitive_results"],
                        f"{root}/cluster/nonrepetitive_indels", contig, start, end, is_scan_job)
    if has_workflow4:
        for model in MODELS:
            copy_cluster_tables(f"{source_dirs['repetitive']}/{model}", f"{source_dirs['repetitive_results']}/{model}",
                                f"{root}/cluster/repetitive_indels/{model}", contig, start, end, is_scan_job)

    relative = lambda path: os.path.relpath(path, root)
    with open(f"{root}/run_local.sh", "w") as handle:
        handle.write(replay_script(
            contig, start, end, out, species, sample, is_scan_job, has_workflow4,
            bam=relative(alignment), assembly=relative(assembly),
            coverage=relative(target.path("coverage", species, sample)),
            tissue=relative(tissue) if tissue != "NA" else "NA",
            mask=relative(f"{target_dirs['nonrepetitive']}/mask/{contig}.bed"),
            beds={model: relative(target.path("tandem_repeats", species, sample, contig=contig, model=model))
                  for model in MODELS},
            locations={workflow: relative(f"{target_dirs[workflow]}/haploid/locations.tsv") for workflow in workflows}))
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


def replay_script(contig, start, end, out, species, sample, is_scan_job, has_workflow4,
                  bam, assembly, coverage, tissue, mask, beds, locations):
    """run_local.sh: replay workflows 3 and 4 on the window and compare with the cluster."""
    script = f"""#!/usr/bin/env bash
# Replays workflows 3 and 4 on {contig}:{start}-{end} and compares each table
# with the cluster's, by content (scripts/compare_tables.py). Exits 1 if any
# table differs. For software testing only.
#
# REPO is the checkout whose scripts/ to run; by default the one this folder
# sits in ({out}/{species}/{sample}/ under the repository root). PYTHON is the
# Python to run them with, one with pysam and pandas. The scripts' own messages
# go to local/replay.log.
set -euo pipefail
cd "$(dirname "$0")"
REPO="${{REPO:-$(cd ../../../.. && pwd)}}"
PYTHON="${{PYTHON:-python}}"
mkdir -p local
: > local/replay.log
trap 'echo "replay failed; see $(pwd)/local/replay.log" >&2' ERR
status=0
compare() {{ "$PYTHON" "$REPO/scripts/compare_tables.py" "$1" "$2" --label "$3" || status=1; }}

# Workflow 3, step 2: the scan.
W3=local/nonrepetitive_indels
mkdir -p $W3
"$PYTHON" "$REPO/scripts/nonrepetitive_scan.py" {mask} {contig} {bam} {assembly} {coverage} {tissue} {start} {end} \\
    $W3/calls.tsv $W3/opportunities.tsv $W3/reads.fastq 2>> local/replay.log
# Step 5: the haploid check, against the cluster's haploid alignments.
"$PYTHON" "$REPO/scripts/filtering_haploid_mappping.py" $W3/calls.tsv {locations['nonrepetitive']} \\
    > $W3/merged_indels.tsv 2>> local/replay.log
# Step 7: final filtering.
"$PYTHON" "$REPO/scripts/filtering_small_indels.py" $W3/merged_indels.tsv -o $W3/indels.tsv \\
    --haploid-check True --min-quality {MIN_KMERQ} --report $W3/filter_report.tsv 2>> local/replay.log
"$PYTHON" "$REPO/scripts/filtering_opportunities.py" $W3/opportunities.tsv -o $W3/opportunities_q{MIN_KMERQ}.tsv \\
    --min-kmerq {MIN_KMERQ} 2>> local/replay.log

C3=cluster/nonrepetitive_indels
compare $W3/calls.tsv $C3/calls.tsv "workflow 3 calls"
"""
    if is_scan_job:
        script += 'compare $W3/opportunities.tsv $C3/opportunities.tsv "workflow 3 opportunities"\n'
    script += """compare $W3/merged_indels.tsv $C3/merged_indels.tsv "workflow 3 calls after the haploid check"
compare $W3/indels.tsv $C3/indels.tsv "workflow 3 final calls"
"""
    if has_workflow4:
        model_args = " \\\n    ".join(f"--model {model} {beds[model]} $W4/{model}/calls.tsv $W4/{model}/opportunities.tsv"
                                      for model in MODELS)
        script += f"""
# Workflow 4, step 1: the scan, every model in one pass.
W4=local/repetitive_indels
mkdir -p {" ".join(f"$W4/{model}" for model in MODELS)}
"$PYTHON" "$REPO/scripts/tandem_repeat_scan.py" {contig} {bam} {assembly} {coverage} {tissue} {start} {end} $W4/reads.fastq \\
    {model_args} 2>> local/replay.log
for model in {" ".join(MODELS)}; do
    # Step 4: the haploid check, against the cluster's haploid alignments.
    "$PYTHON" "$REPO/scripts/filtering_haploid_mappping.py" $W4/$model/calls.tsv {locations['repetitive']} \\
        > $W4/$model/merged_indels.tsv 2>> local/replay.log
    # Step 6: final filtering.
    "$PYTHON" "$REPO/scripts/filtering_small_indels.py" $W4/$model/merged_indels.tsv -o $W4/$model/indels.tsv \\
        --haploid-check True --min-quality {MIN_KMERQ} --min-read-end-distance {MIN_READ_END_DISTANCE} \\
        --report $W4/$model/filter_report.tsv 2>> local/replay.log
    "$PYTHON" "$REPO/scripts/filtering_opportunities.py" $W4/$model/opportunities.tsv \\
        -o $W4/$model/opportunities_q{MIN_KMERQ}.tsv --min-kmerq {MIN_KMERQ} 2>> local/replay.log

    C4=cluster/repetitive_indels/$model
    compare $W4/$model/calls.tsv $C4/calls.tsv "workflow 4 $model calls"
"""
        if is_scan_job:
            script += '    compare $W4/$model/opportunities.tsv $C4/opportunities.tsv "workflow 4 $model opportunities"\n'
        script += """    compare $W4/$model/merged_indels.tsv $C4/merged_indels.tsv "workflow 4 $model calls after the haploid check"
    compare $W4/$model/indels.tsv $C4/indels.tsv "workflow 4 $model final calls"
done
"""
    return script + "\nexit $status\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("samples", nargs="+", metavar="species/sample[:contig:start-end]",
                        help="the samples to export, each with its window or none")
    parser.add_argument("--window-size", type=int, default=SCAN_WINDOW,
                        help="size of the window chosen for a sample given without one (default: one scan job)")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--out", default="mirror",
                        help="each window goes to {out}/{species}/{sample}/{contig}_{start}_{end}/")
    parser.add_argument("--keep-kinetics", action="store_true", help="keep the per-base kinetics tags")
    args = parser.parse_args()

    source = Pipeline(args.config)
    for spec in args.samples:
        species, sample, window = parse_sample(spec)
        if window is None:
            window = busiest_window(workflow_dirs(source, species, sample)["nonrepetitive_results"],
                                    source.long_contigs(species, sample), args.window_size)
        export_window(source, species, sample, *window, args.out, args.keep_kinetics)


if __name__ == "__main__":
    main()
