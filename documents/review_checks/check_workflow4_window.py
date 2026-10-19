"""
How workflow 4's scan treats real reads, on one window mirrored from the
cluster with scripts/export_local_subset.py. These are checks of software
behaviour; the mirror is for software testing only, so no rate is drawn from
it (../workflow4_review.md takes those from results/).

    python documents/review_checks/check_workflow4_window.py mirror/human/chk/h1tg000001l_0_10000000

It reads the window's BAM and M0 tracts, and, if the window has been replayed
(run_local.sh), the replay's M0 calls.
"""
import collections
import os
import sys

import pysam
import yaml
from intervaltree import IntervalTree

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
from utilities import canonical_repeat_indel  # noqa: E402

if len(sys.argv) != 2:
    sys.exit(__doc__)
root = sys.argv[1].rstrip("/")
contig, start, end = os.path.basename(root).rsplit("_", 2)
start, end = int(start), int(end)
with open(f"{root}/config.yaml") as handle:
    config = yaml.safe_load(handle)
(species, samples), = config["samples"].items()
sample, kind = samples[0], config["assembly_type"]
BAM = f"{root}/bam/{species}/{sample}/kinetics/{sample}_kinetics_{kind}.bam"
FASTA = f"{root}/fasta/{species}/{sample}/{kind}_assembly/{sample}_{kind}.fa"

tracts = IntervalTree()
starts = {}
with open(f"{root}/output/repetitive_dna/{species}/{sample}/{contig}_M0.bed") as handle:
    handle.readline()
    for line in handle:
        fields = line.split("\t")
        s, e, motif, k = int(fields[1]), int(fields[2]), fields[3], int(fields[4])
        tracts[s:e] = (motif, k)
        starts[s] = (motif, k, e)

reference = pysam.FastaFile(FASTA).fetch(contig, 0, end + 100_000).upper()


def indels(read):
    """(anchor, sequence) of each InDel in the read: the base before the gap, and the bases inserted or deleted."""
    query, ref = 0, read.reference_start
    sequence = read.query_sequence
    for op, length in read.cigartuples:
        if op in (0, 7, 8):          # match
            query, ref = query + length, ref + length
        elif op == 1:                # insertion
            yield ref - 1, sequence[query:query + length]
            query += length
        elif op == 2:                # deletion
            yield ref - 1, reference[ref:ref + length]
            ref += length
        elif op == 4:                # soft clip
            query += length
        elif op == 3:                # skipped reference
            ref += length


# Every whole-unit InDel of an M0 tract's motif, in every primary read, by where
# its anchor (the base before the gap, as pbmm2 places it) lies.
places = collections.Counter()
with pysam.AlignmentFile(BAM) as bam:
    for read in bam.fetch(contig, start, end):
        if read.is_secondary or read.is_supplementary:
            continue
        for anchor, sequence in indels(read):
            if not start <= anchor < end:
                continue
            before = starts.get(anchor + 1)
            inside = tracts[anchor]
            if before and canonical_repeat_indel(sequence) == before[0]:
                places["anchored on the base before the tract"] += 1
            elif inside:
                motif, k = next(iter(inside)).data
                if len(sequence) == k and canonical_repeat_indel(sequence) == motif:
                    places["anchored inside the tract"] += 1

print(f"{species}/{sample} {contig}:{start}-{end}")
print("One-unit InDels of an M0 tract's motif, in the reads (primary alignments):")
for place, n in places.most_common():
    print(f"  {place}: {n:,}")
print("Workflow 4 scans only anchors inside a tract, so it sees only the second kind.")

calls = f"{root}/local/repetitive_indels/M0/calls.tsv"
if os.path.exists(calls):
    offsets = collections.Counter()
    with open(calls) as handle:
        header = handle.readline().lstrip("#").rstrip("\n").split("\t")
        position, string = header.index("position"), header.index("indel_string")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            anchor = int(fields[position]) - 1
            for tract in tracts[anchor]:
                offsets["first base of the tract" if anchor == tract.begin else
                        "last unit of the tract" if anchor >= tract.end - tract.data[1] else "elsewhere"] += 1
    print("\nWorkflow 4's M0 calls in the replay, by where their anchor lies in the tract:")
    for place, n in offsets.most_common():
        print(f"  {place}: {n:,}")
