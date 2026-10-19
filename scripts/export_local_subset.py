"""
Copy one window of a sample's aligned reads from the cluster, with everything
workflow 3 needs to run on it, so the pipeline can be mirrored locally.

Run it on the cluster, from the repository root, once workflows 1-3 have run
for the sample:

    python scripts/export_local_subset.py human chk h1tg000001l:0-10000000

The window is 0-based and half-open, like the scan's windows. A window that is
exactly one scan job (it starts at a multiple of 10 Mb and is 10 Mb long, or
runs to the contig's end) can be compared with the cluster's output call by
call and opportunity by opportunity. Leave the window out to take the scan
window with the most unique calls in results/.

Every alignment that overlaps the window is kept, whatever its flags, because
that is what the scan sees: the coverage test compares each position's depth
with the sample's mean, so the reads cannot be thinned.

It needs only pysam and pyyaml. A 10 Mb window of a 30x sample is roughly
0.5 GB without the per-base kinetics tags, which are dropped unless
--keep-kinetics. Each window gets its own folder,
mirror/{species}/{sample}/{contig}_{start}_{end}/. Copy it to the same place
in a local checkout:

    rsync -av CLUSTER:REPO/mirror/human/chk/ mirror/human/chk/

and run its run_local.sh there. The folder has the cluster's layout, so the
pipeline's scripts run on it unchanged:

  config.yaml        paths pointing into this folder
  fasta/             the window's contig; the sequence outside the window and a
                     100 kb margin is N, so coordinates match the cluster
  bam/               every alignment overlapping the window, indexed
  output/precomute/  mean coverage, unresolved regions, tissue labels
  output/repetitive_dna/
                     the M0, M1 and M2 tracts near the window
  output/nonrepetitive_indels/
                     the mask near the window, and the haploid alignment records
                     (locations.tsv) of the window's reads, so the haploid check
                     runs without remapping
  cluster/           the cluster's results for the window, to compare with: the
                     scan job's calls and opportunities, the calls after the
                     haploid check, and the final calls
  run_local.sh       reruns workflow 3 on the window and compares
"""
import argparse
import os
import shutil
import sys
import time
from collections import Counter
from os.path import expanduser

import pysam
import yaml

KINETICS_TAGS = ("fi", "ri", "fp", "rp", "ip", "pw")
SCAN_WINDOW = 10_000_000      # workflow 3 scans each contig in windows of this size
MARGIN = 100_000              # real sequence kept on each side of the window
STARTED = time.time()


def log(message):
    print(f"[{time.time() - STARTED:5.0f} s] {message}", file=sys.stderr, flush=True)


def cluster_paths(config, species, sample):
    """The cluster's files for the sample, built as in workflows/indel_mutation_call.py."""
    fasta = expanduser(config["paths"]["fasta"])
    bam = expanduser(config["paths"]["bam"])
    output = os.path.abspath(config["paths"]["output"])
    results = os.path.abspath(config["paths"]["results"])
    kind = config["assembly_type"]
    overrides = (config.get("alignment_overrides") or {}).get(species) or {}
    tissue_samples = (config.get("tissue_samples") or {}).get(species) or []
    precompute = f"{output}/precomute/{species}/{sample}"
    return {
        "assembly": f"{fasta}/{species}/{sample}/{kind}_assembly/{sample}_{kind}.fa",
        "alignment": expanduser(overrides.get(
            sample, f"{bam}/{species}/{sample}/kinetics/{sample}_kinetics_{kind}.bam")),
        "precompute": precompute,
        "tissue": f"{precompute}/sorted_tissue_classification.tsv" if sample in tissue_samples else None,
        "repeats": f"{output}/repetitive_dna/{species}/{sample}",
        "workdir": f"{output}/nonrepetitive_indels/{species}/{sample}",
        "results": f"{results}/nonrepetitive_indels/{species}/{sample}",
    }


def mirror_paths(root, species, sample, kind):
    """The same files under the mirror folder."""
    precompute = f"{root}/output/precomute/{species}/{sample}"
    return {
        "assembly": f"{root}/fasta/{species}/{sample}/{kind}_assembly/{sample}_{kind}.fa",
        "alignment": f"{root}/bam/{species}/{sample}/kinetics/{sample}_kinetics_{kind}.bam",
        "precompute": precompute,
        "tissue": f"{precompute}/sorted_tissue_classification.tsv",
        "repeats": f"{root}/output/repetitive_dna/{species}/{sample}",
        "workdir": f"{root}/output/nonrepetitive_indels/{species}/{sample}",
    }


def busiest_window(results, contigs):
    """The scan window holding the most unique calls in the final results."""
    counts = Counter()
    with open(f"{results}/indels.tsv") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        contig, position, unique = (header.index(c) for c in ("contig", "position", "unique"))
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[unique] == "unique" and fields[contig] in contigs:
                counts[(fields[contig], (int(fields[position]) - 1) // SCAN_WINDOW)] += 1
    if not counts:
        sys.exit(f"No unique calls in {results}/indels.tsv; give a window.")
    (contig, index), n = counts.most_common(1)[0]
    log(f"window with the most unique calls: {n} calls")
    return contig, index * SCAN_WINDOW, min((index + 1) * SCAN_WINDOW, contigs[contig])


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


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("species")
    parser.add_argument("sample")
    parser.add_argument("window", nargs="?", help="contig:start-end, 0-based half-open")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--out", default="mirror",
                        help="the window goes to {out}/{species}/{sample}/{contig}_{start}_{end}/")
    parser.add_argument("--keep-kinetics", action="store_true", help="keep the per-base kinetics tags")
    args = parser.parse_args()

    with open(args.config) as handle:
        config = yaml.safe_load(handle)
    kind = config["assembly_type"]
    source = cluster_paths(config, args.species, args.sample)
    with open(f"{source['assembly']}.fai") as handle:
        contigs = {f[0]: int(f[1]) for f in (line.split("\t") for line in handle)
                   if int(f[1]) > config["min_contig_length"]}

    if args.window:
        contig, span = args.window.rsplit(":", 1)
        start, end = (int(x.replace(",", "")) for x in span.split("-"))
    else:
        contig, start, end = busiest_window(source["results"], contigs)
    if contig not in contigs or not 0 <= start < end <= contigs[contig]:
        sys.exit(f"{contig}:{start}-{end} is not inside a contig longer than min_contig_length")
    is_scan_job = start % SCAN_WINDOW == 0 and end == min(start + SCAN_WINDOW, contigs[contig])
    root = os.path.abspath(os.path.join(args.out, args.species, args.sample, f"{contig}_{start}_{end}"))
    target = mirror_paths(root, args.species, args.sample, kind)
    log(f"{args.species}/{args.sample} {contig}:{start}-{end}"
        f" ({'one scan job' if is_scan_job else 'not a scan job: calls compared, opportunities not'})")

    # The reads: every alignment overlapping the window, sorted as fetched.
    log("alignments")
    os.makedirs(os.path.dirname(target["alignment"]), exist_ok=True)
    names, records = set(), 0
    with pysam.AlignmentFile(source["alignment"], "rb") as bam, \
            pysam.AlignmentFile(target["alignment"], "wb", template=bam) as out:
        for read in bam.fetch(contig, start, end):
            if not args.keep_kinetics:
                for tag in KINETICS_TAGS:
                    if read.has_tag(tag):
                        read.set_tag(tag, None)
            out.write(read)
            names.add(read.query_name)
            records += 1
    pysam.index(target["alignment"])
    log(f"  {records:,} alignment records of {len(names):,} reads")

    # The contig at full length, so coordinates match; N away from the window.
    log("assembly")
    os.makedirs(os.path.dirname(target["assembly"]), exist_ok=True)
    length = contigs[contig]
    left, right = max(0, start - MARGIN), min(length, end + MARGIN)
    with pysam.FastaFile(source["assembly"]) as fasta:
        sequence = fasta.fetch(contig, left, right)
    with open(target["assembly"], "w") as out:
        out.write(f">{contig}\n")
        full_lines, rest = divmod(left, 60)
        out.writelines("N" * 60 + "\n" for _ in range(full_lines))
        body = "N" * rest + sequence + "N" * (length - right)
        out.writelines(body[i:i + 60] + "\n" for i in range(0, len(body), 60))
    pysam.faidx(target["assembly"])

    # Workflows 1 and 2: coverage, unresolved regions, tissue labels, repeat tracts.
    log("precomputed inputs and repeat tracts")
    os.makedirs(target["precompute"], exist_ok=True)
    shutil.copy(f"{source['precompute']}/avg_coverage.txt", f"{target['precompute']}/avg_coverage.txt")
    with open(f"{target['precompute']}/long_contigs.bed", "w") as out:
        out.write(f"{contig}\t0\t{length}\n")
    copy_overlapping(f"{source['precompute']}/unresolved_haplotypes/{contig}.bed",
                     f"{target['precompute']}/unresolved_haplotypes/{contig}.bed", contig, left, right, header=False)
    if source["tissue"]:
        with open(source["tissue"]) as handle, open(target["tissue"], "w") as out:
            out.write(handle.readline())
            out.writelines(line for line in handle if line.split("\t", 1)[0] in names)
    for model in ("M0", "M1", "M2"):
        copy_overlapping(f"{source['repeats']}/{contig}_{model}.bed",
                         f"{target['repeats']}/{contig}_{model}.bed", contig, left, right, header=True)

    # Workflow 3: the mask, and the haploid alignments of the window's reads.
    log("mask and haploid alignment records")
    copy_overlapping(f"{source['workdir']}/mask/{contig}.bed",
                     f"{target['workdir']}/mask/{contig}.bed", contig, left, right, header=False)
    os.makedirs(f"{target['workdir']}/haploid", exist_ok=True)
    with open(f"{source['workdir']}/haploid/locations.tsv") as handle, \
            open(f"{target['workdir']}/haploid/locations.tsv", "w") as out:
        out.writelines(line for line in handle if line.split("\t", 1)[0] in names)

    # The cluster's results for the window, to compare a local run with.
    log("the cluster's results for the window")
    compare = f"{root}/cluster"
    os.makedirs(compare, exist_ok=True)
    chunk = f"{contig}_{start}_{end}"
    if is_scan_job:
        for kind_of_table in ("calls", "opportunities"):
            path = f"{source['workdir']}/{kind_of_table}/chunks/{chunk}.tsv"
            if os.path.exists(path):
                shutil.copy(path, f"{compare}/{kind_of_table}.tsv")
    else:
        copy_calls_in_window(f"{source['workdir']}/calls/{contig}.tsv", f"{compare}/calls.tsv", contig, start, end)
    copy_calls_in_window(f"{source['workdir']}/merged/indels.tsv", f"{compare}/merged_indels.tsv", contig, start, end)
    copy_calls_in_window(f"{source['results']}/indels.tsv", f"{compare}/indels.tsv", contig, start, end)

    # A config.yaml for the mirror, so the workflows' path logic finds these files.
    with open(f"{root}/config.yaml", "w") as out:
        yaml.safe_dump({
            "samples": {args.species: [args.sample]},
            "tissue_samples": {args.species: [args.sample]} if source["tissue"] else {},
            "methylation_samples": {},
            "alignment_overrides": {},
            "assembly_type": kind,
            "min_contig_length": 0,
            "paths": {"fasta": f"{root}/fasta", "bam": f"{root}/bam", "output": f"{root}/output",
                      "results": f"{root}/results", "plots": f"{root}/plots"},
        }, out, sort_keys=False)

    relative = lambda path: os.path.relpath(path, root)
    tissue = relative(target["tissue"]) if source["tissue"] else "NA"
    workdir = relative(target["workdir"])
    with open(f"{root}/run_local.sh", "w") as out:
        out.write(f"""#!/usr/bin/env bash
# Reruns workflow 3 on {contig}:{start}-{end} and compares with the cluster.
# REPO is the checkout whose scripts/ to run; by default the one this folder
# sits in ({args.out}/{args.species}/{args.sample}/ under the repository root).
set -euo pipefail
cd "$(dirname "$0")"
REPO="${{REPO:-$(cd ../../.. && pwd)}}"
mkdir -p local

# Step 2: the scan.
python "$REPO/scripts/nonrepetitive_scan.py" {workdir}/mask/{contig}.bed {contig} \\
    {relative(target['alignment'])} {relative(target['assembly'])} \\
    {relative(target['precompute'])}/avg_coverage.txt {tissue} {start} {end} \\
    local/calls.tsv local/opportunities.tsv local/reads.fastq
# Step 5: the haploid check, against the cluster's haploid alignments.
python "$REPO/scripts/filtering_haploid_mappping.py" local/calls.tsv {workdir}/haploid/locations.tsv \\
    > local/merged_indels.tsv
# Step 7: final filtering.
python "$REPO/scripts/filtering_small_indels.py" local/merged_indels.tsv -o local/indels.tsv \\
    --haploid-check True --min-quality 93 --report local/filter_report.tsv
python "$REPO/scripts/filtering_opportunities.py" local/opportunities.tsv \\
    -o local/opportunities_q93.tsv --min-kmerq 93

# Compare with the cluster, ignoring row order.
same() {{ cmp -s <(tail -n +2 "$1" | sort) <(tail -n +2 "$2" | sort) \\
          && echo "same:      $1" || echo "DIFFERENT: $1 vs $2"; }}
same local/calls.tsv cluster/calls.tsv
[ -f cluster/opportunities.tsv ] && same local/opportunities.tsv cluster/opportunities.tsv
echo "calls after the haploid check: local $(($(wc -l < local/merged_indels.tsv) - 1)),"\\
     "cluster $(($(wc -l < cluster/merged_indels.tsv) - 1))"
echo "final calls:                   local $(($(wc -l < local/indels.tsv) - 1)),"\\
     "cluster $(($(wc -l < cluster/indels.tsv) - 1))"
""")
    os.chmod(f"{root}/run_local.sh", 0o755)

    sizes = {os.path.relpath(os.path.join(d, f), root): os.path.getsize(os.path.join(d, f))
             for d, _, files in os.walk(root) for f in files}
    with open(f"{root}/README.txt", "w") as out:
        out.write(f"Made by scripts/export_local_subset.py on {time.strftime('%Y-%m-%d %H:%M')}\n"
                  f"command: {' '.join(sys.argv)}\n"
                  f"window: {contig}:{start}-{end} ({'one scan job' if is_scan_job else 'not a scan job'})\n"
                  f"BAM: {source['alignment']}\n"
                  f"alignment records: {records:,} of {len(names):,} reads\n"
                  f"kinetics kept: {args.keep_kinetics}\n\nfiles (MB):\n")
        out.writelines(f"  {name}: {size / 1e6:.1f}\n" for name, size in sorted(sizes.items()))
    log(f"done: {sum(sizes.values()) / 1e6:.0f} MB in {root}")


if __name__ == "__main__":
    main()
