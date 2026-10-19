"""
The real-data tests of ../nonrepetitive_indel_review.md, on one window mirrored
from the cluster with scripts/export_local_subset.py.

    python documents/review_checks/check_mirror_window.py mirror/human/chk/h1tg000001l_0_10000000

Run it from the repository root in an environment with pysam and pandas
(indels-local, or .venv). It only reads the window's files: bam/, fasta/,
output/ and the cluster's results in cluster/. It takes a few minutes.
"""
import bisect
import collections
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
import pysam
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
from filtering_haploid_mappping import check_call  # noqa: E402
from old_haploid_walk import find_indel_position_in_read  # noqa: E402
from utilities import reverse_complement, subread_array  # noqa: E402

warnings.simplefilter("ignore")
pysam.set_verbosity(0)          # the copied index can be older than the BAM; harmless

if len(sys.argv) != 2:
    sys.exit(__doc__)
root = sys.argv[1].rstrip("/")
contig, start, end = os.path.basename(root).rsplit("_", 2)
start, end = int(start), int(end)
with open(f"{root}/config.yaml") as handle:
    config = yaml.safe_load(handle)
species, samples = next(iter(config["samples"].items()))
sample, kind = samples[0], config["assembly_type"]
BAM = f"{root}/bam/{species}/{sample}/kinetics/{sample}_kinetics_{kind}.bam"
FASTA = f"{root}/fasta/{species}/{sample}/{kind}_assembly/{sample}_{kind}.fa"
WORK = f"{root}/output/nonrepetitive_indels/{species}/{sample}"
PRE = f"{root}/output/precomute/{species}/{sample}"
# Workflow 3's tables from the cluster; windows exported before workflow 4 was
# mirrored too keep them directly under cluster/.
CLUSTER = f"{root}/cluster/nonrepetitive_indels"
if not os.path.isdir(CLUSTER):
    CLUSTER = f"{root}/cluster"

bam = pysam.AlignmentFile(BAM)
reference = pysam.FastaFile(FASTA).fetch(contig, 0, end + 100_000).upper()
calls = pd.read_csv(f"{CLUSTER}/calls.tsv", sep="\t", low_memory=False).rename(columns={"#contig": "contig"})
calls["size"] = calls.indel_string.str.len()
mask_bed = pd.read_csv(f"{WORK}/mask/{contig}.bed", sep="\t", header=None, names=["c", "s", "e"])
mask_starts, mask_ends = mask_bed.s.tolist(), mask_bed.e.tolist()
print(f"{species}/{sample} {contig}:{start}-{end}: {len(calls):,} calls before the haploid check\n")


def section(title):
    print(f"\n=== {title}")


def masked(position):
    i = bisect.bisect_right(mask_starts, position) - 1
    return i >= 0 and mask_ends[i] > position


def decode(sa):
    """Run-length pairs (length, coverage) -> per-base subread coverage."""
    return np.repeat(np.array(sa[1::2]), np.array(sa[0::2]))


def reference_position(read, query_position):
    q, r = 0, read.reference_start
    for op, n in read.cigartuples:
        if op in (0, 7, 8):
            if q <= query_position < q + n:
                return r + query_position - q
            q += n
            r += n
        elif op in (1, 4):
            if q <= query_position < q + n:
                return None
            q += n
        elif op in (2, 3):
            r += n
    return None


def call_record(call):
    """The alignment record that carries a call."""
    anchor, read_pos = call.position - 1, call.start - 1
    return next((r for r in bam.fetch(call.contig, anchor, anchor + 1)
                 if r.query_name == call.readID and r.query_length == call.read_length
                 and reference_position(r, read_pos) == anchor), None)


records = [call_record(c) for c in calls.itertuples(index=False)]

# ---------------------------------------------------------------------------
section("S1: the subread tags (sm, sx, sa) and the read's orientation")
# ---------------------------------------------------------------------------
per_read = {"+": [], "-": []}
for r in bam.fetch(contig, start, end):
    if r.is_secondary or r.is_supplementary:
        continue
    aligned = np.array(r.query_qualities, float)
    native = np.array(r.get_forward_qualities(), float)
    sx = np.array(r.get_tag("sx"), float)
    support = np.array(r.get_tag("sm"), float) / np.maximum(decode(r.get_tag("sa")), 1)
    if sx.std() == 0 or aligned.std() == 0:
        continue
    per_read["-" if r.is_reverse else "+"].append(
        [np.corrcoef(aligned, support)[0, 1], np.corrcoef(native, support)[0, 1],
         np.corrcoef(aligned, sx)[0, 1], np.corrcoef(native, sx)[0, 1]])
print("Median correlation per read between base quality and the subread tags, with the tags taken")
print("as the scan indexes them (aligned orientation) and in the read's own orientation:")
print(pd.DataFrame({s: np.median(np.array(v), axis=0) for s, v in per_read.items()},
                   index=["support, as scanned", "support, own orientation",
                          "sx, as scanned", "sx, own orientation"]).T.assign(
    reads=[len(per_read["+"]), len(per_read["-"])]).round(3).to_string())

rows = []
for c, rec in zip(calls.itertuples(index=False), records):
    if rec is None:
        rows.append({})
        continue
    sm, sx, sa = rec.get_tag("sm"), rec.get_tag("sx"), rec.get_tag("sa")
    length, i = rec.query_length, c.start - 1
    own = (lambda j: length - 1 - j) if rec.is_reverse else (lambda j: j)
    scanned = (sm[i + 1], sx[i + 1], subread_array(sa, i), sm[i], sx[i])
    correct = (sm[own(i + 1)], sx[own(i + 1)], subread_array(sa, own(i)), sm[own(i)], sx[own(i)])
    rows.append({"reproduces": scanned == (c.sm, c.sx, c.sa, c.sm_pre, c.sx_pre),
                 **{f"wrong_{k}": a != b for k, a, b in zip(("sm", "sx", "sa", "sm_pre", "sx_pre"), scanned, correct)},
                 "sx_scanned": scanned[1], "sx_correct": correct[1],
                 "full_scanned": sm[i + 1] == subread_array(sa, i + 1),
                 "full_correct": correct[0] == subread_array(sa, own(i + 1)),
                 "used": i + 1, "correct_index": own(i + 1),
                 "recorded": f"{c.sm}/{c.sa} {c.sx} {c.sx_pre}",
                 "correct": f"{correct[0]}/{correct[2]} {correct[1]} {correct[4]}",
                 "flag": rec.flag, "mapq": rec.mapping_quality})
tags = pd.concat([calls, pd.DataFrame(rows, index=calls.index)], axis=1)
print(f"\nRecomputing the tags the scan's way reproduces the recorded sm, sx, sa, sm_pre, sx_pre for "
      f"{tags.reproduces.mean():.1%} of calls")
print("\nAt the call sites (all calls, any kmerQ):")
print(tags.groupby("diploid_direction").agg(
    calls=("sx_correct", "size"),
    sx_1_or_more_as_scanned=("sx_scanned", lambda s: (s >= 1).mean()),
    sx_1_or_more_correct=("sx_correct", lambda s: (s >= 1).mean()),
    full_support_as_scanned=("full_scanned", "mean"),
    full_support_correct=("full_correct", "mean")).round(3).to_string())
reverse = tags[tags.diploid_direction == "-"]
print("\nReverse-strand calls: share of rows whose recorded value is wrong: "
      + ", ".join(f"{k} {reverse[f'wrong_{k}'].mean():.0%}" for k in ("sm", "sx", "sa", "sm_pre", "sx_pre")))
q93 = tags[tags.kmerQ >= 93]
print(f"\nkmerQ >= 93: {len(q93)} calls, {(q93.unique == 'unique').sum()} unique")
print(q93.groupby(["unique", "diploid_direction"]).agg(
    calls=("sx_correct", "size"), full_support_as_scanned=("full_scanned", "mean"),
    full_support_correct=("full_correct", "mean")).round(2).to_string())
print("\nUnique kmerQ >= 93 calls, recorded against correct (sm/sa sx sx_pre):")
print(q93[q93.unique == "unique"].sort_values(["diploid_direction", "position"])[
    ["diploid_direction", "indel_type", "indel_string", "position", "recorded", "correct",
     "used", "correct_index", "read_length"]].to_string(index=False))

# ---------------------------------------------------------------------------
section("H1-H4: the haploid check")
# ---------------------------------------------------------------------------
merged = pd.read_csv(f"{CLUSTER}/merged_indels.tsv", sep="\t", low_memory=False)
first = {}
with open(f"{WORK}/haploid/locations.tsv") as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        first.setdefault(fields[0], fields)          # the record the pipeline uses


def events(cigar):
    """(op, length, read position) of every I and D, counting I, M, =, X and S as read bases."""
    q, out = 0, []
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=XS":
            q += n
        elif op == "I":
            out.append(("I", n, q))
            q += n
        elif op == "D":
            out.append(("D", n, q))
    return out


def walk_with_h1_fixed(cigar, position, op_wanted, size):
    """The pipeline's walk, with inserted bases counted (H1) and nothing else changed."""
    q = 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=XS":
            q += n
        elif op in "ID":
            if op == op_wanted and q >= position:
                return n == size
            if op == "I":
                q += n
    return False


def slide(read, kind, gap, sequence):
    """The leftmost and rightmost equivalent placements of an InDel, in the read's aligned orientation."""
    n, length = len(sequence), len(read)
    low = high = gap
    if kind == "D":
        s = sequence
        while low - 1 >= 0 and read[low - 1] == s[-1]:
            s = s[-1] + s[:-1]
            low -= 1
        s = sequence
        while high < length and read[high] == s[0]:
            s = s[1:] + s[0]
            high += 1
    else:
        while low - 1 >= 0 and low - 1 + n < length and read[low - 1] == read[low - 1 + n]:
            low -= 1
        while high + n < length and read[high] == read[high + n]:
            high += 1
    return low, high


rows = []
for c in merged.itertuples(index=False):
    rec = first.get(c.readID)
    if rec is None or int(rec[1]) not in (0, 16):
        continue
    cigar = rec[5]
    flipped = ("-" if int(rec[1]) == 16 else "+") != c.diploid_direction
    kind, n, length = ("I" if c.indel_type == "insertion" else "D"), len(c.indel_string), c.read_length
    position_used = length - c.start - 1 if flipped else c.start
    now, _, _, _ = find_indel_position_in_read(cigar, position_used, c.indel_type, c.indel_string)
    gap = c.start                                    # the read base after the anchor
    mapped = (length - gap if kind == "D" else length - gap - n) if flipped else gap
    same = [q for op, ln, q in events(cigar) if op == kind and ln == n]
    read = next((r.query_sequence for r in bam.fetch(c.contig, c.position - 1, c.position)
                 if r.query_name == c.readID and r.query_length == length), None)
    if read is None:
        continue
    low, high = slide(read, kind, gap, c.indel_string if kind == "D" else read[gap:gap + n])
    if flipped:
        low, high = (length - high, length - low) if kind == "D" else (length - high - n, length - low - n)
    # The pipeline's check as fixed since the review takes the read as sequenced.
    forward = read if c.diploid_direction == "+" else reverse_complement(read)
    fixed, _, _, _ = check_call(forward, int(rec[1]), cigar, c.start, c.indel_type, c.indel_string, c.diploid_direction)
    rows.append({"flipped": flipped, "type": c.indel_type, "size": n, "recorded": c.haploid_check,
                 "now": now, "H1 fixed": walk_with_h1_fixed(cigar, position_used, kind, n),
                 "H1+H2 fixed": mapped in same, "H1+H2+H3 fixed": any(low <= q <= high for q in same),
                 "pipeline fixed": fixed, "contig_hap": rec[2]})
haploid = pd.DataFrame(rows)
print(f"The walk as reviewed (old_haploid_walk.py) reproduces the recorded haploid_check for "
      f"{(haploid.now == haploid.recorded).mean():.1%} of {len(haploid):,} calls")
print(f"InDels per read in the haploid CIGARs (median): {np.median([len(events(v[5])) for v in first.values()]):.0f}")
print("\nShare of calls passing the haploid check, now and with the bugs fixed one after another:")
print(haploid.groupby(["flipped", "type", "size"]).agg(
    calls=("now", "size"), now=("now", "mean"), H1_fixed=("H1 fixed", "mean"),
    H1_H2_fixed=("H1+H2 fixed", "mean"), H1_H2_H3_fixed=("H1+H2+H3 fixed", "mean"),
    pipeline_fixed=("pipeline fixed", "mean")).round(3).to_string())
print(f"The pipeline's check as fixed agrees with H1+H2+H3 fixed for "
      f"{(haploid['pipeline fixed'] == haploid['H1+H2+H3 fixed']).mean():.1%} of calls")
wrong = haploid[haploid.now & ~haploid["H1+H2+H3 fixed"]]
print(f"\nH4: calls passing now with no matching event where the InDel is: {len(wrong)} of {haploid.now.sum():,}")
print(f"Overall pass rate: now {haploid.now.mean():.1%}, with H1-H3 fixed {haploid['H1+H2+H3 fixed'].mean():.1%}")
print("Haploid contigs of the calls' reads:", haploid.groupby("flipped").contig_hap.value_counts().to_dict())

# ---------------------------------------------------------------------------
section("L1: InDels from masked repeats anchored outside the mask")
# ---------------------------------------------------------------------------
window_mask = mask_bed[(mask_bed.e > start) & (mask_bed.s < end)]
masked_bp = (window_mask.e.clip(upper=end) - window_mask.s.clip(lower=start)).sum()
unmasked_bp = (end - start) - masked_bp
chance = sum(1 for s in window_mask.s if start < s <= end and not masked(s - 1)) / unmasked_bp
print(f"The mask covers {masked_bp / 1e6:.2f} Mb of the window; by chance {chance:.2%} of unmasked positions "
      "sit right before a mask interval")
calls["next_masked"] = [masked(p) for p in calls.position]     # 1-based anchor = 0-based base after it
print("\nShare of calls whose first affected base (the one after the anchor) is masked:")
print(calls.groupby(["size", "unique"]).agg(calls=("next_masked", "size"), masked=("next_masked", "mean")).round(3).to_string())
print(calls[calls.kmerQ >= 93].groupby("size").agg(
    calls_kmerQ_93=("next_masked", "size"), masked=("next_masked", "mean")).round(3).to_string())
placements = collections.Counter()
for k, r in enumerate(bam.fetch(contig, start + 200_000, end - 200_000)):
    if k >= 3000:
        break
    seq, q, p = r.query_sequence, 0, r.reference_start
    for op, n in r.cigartuples:
        if op in (0, 7, 8):
            q += n
            p += n
        elif op == 4:
            q += n
        elif op == 2:
            if 0 < p and p + n < len(reference):
                left, right = reference[p - 1] == reference[p + n - 1], reference[p] == reference[p + n]
                placements["both" if left and right else "left" if left else "right" if right else "none"] += 1
            p += n
        elif op == 1:
            if 0 < q and q + n < len(seq):
                left, right = seq[q - 1] == seq[q + n - 1], seq[q] == seq[q + n]
                placements["both" if left and right else "left" if left else "right" if right else "none"] += 1
            q += n
print("\nInDels in 3,000 diploid alignments, by the way they could still slide "
      "(none could slide left if the aligner places gaps leftmost):", dict(placements))

# ---------------------------------------------------------------------------
section("U1: unresolved-haplotype regions")
# ---------------------------------------------------------------------------
sites = collections.defaultdict(list)
for r in bam.fetch(contig, start, end):
    seq, qual, q, p = r.query_sequence, r.query_qualities, 0, r.reference_start
    for op, n in r.cigartuples:
        if op in (0, 7, 8):
            if op != 7:
                for j in range(n):
                    pos = p + j
                    if start <= pos < end and qual[q + j] >= 40 and seq[q + j] != reference[pos]:
                        sites[(pos, seq[q + j])].append((r.reference_start, r.reference_end))
            q += n
            p += n
        elif op in (1, 4):
            q += n
        elif op in (2, 3):
            p += n
hits = sorted((pos, min(s for s, _ in v), max(e for _, e in v)) for (pos, _), v in sites.items() if len(v) > 1)


def regions(fixed):
    out, region_start, region_end = [], 0, 0
    for _, s, e in hits:
        if region_start == 0 and region_end == 0:
            region_start, region_end = s, e
        elif fixed and s <= region_end:
            region_start, region_end = min(region_start, s), max(region_end, e)
        elif not fixed and s < region_end:
            region_end = e                        # as coded: the end can move back, the start never moves
        elif fixed or s > region_end:
            out.append((region_start, region_end))
            region_start, region_end = s, e
    if region_end:
        out.append((region_start, region_end))
    return out


def merged_span(intervals):
    joined = []
    for s, e in sorted(intervals):
        if joined and s <= joined[-1][1]:
            joined[-1][1] = max(joined[-1][1], e)
        else:
            joined.append([s, e])
    return joined, sum(min(e, end) - max(s, start) for s, e in joined if e > start and s < end)


as_coded, as_coded_bp = merged_span(regions(fixed=False))
fixed, fixed_bp = merged_span(regions(fixed=True))
cluster_bed = pd.read_csv(f"{PRE}/unresolved_haplotypes/{contig}.bed", sep="\t", header=None)
_, cluster_bp = merged_span(list(zip(cluster_bed[1], cluster_bed[2])))


def inside(intervals, position):
    i = bisect.bisect_right([s for s, _ in intervals], position) - 1
    return i >= 0 and intervals[i][1] > position


newly = calls[[inside(fixed, p - 1) and not inside(as_coded, p - 1) for p in calls.position]]
print(f"Unresolved bp in the window: cluster {cluster_bp:,}; rebuilt as coded {as_coded_bp:,}; "
      f"with the bounds fixed {fixed_bp:,}")
print(f"Calls the fixed regions would mask: {len(newly)}, of them kmerQ >= 93: {(newly.kmerQ >= 93).sum()}")

# ---------------------------------------------------------------------------
section("M1: exact tandem repeats the mask misses (k = 2-6, at least 3 copies)")
# ---------------------------------------------------------------------------
length = end - start
mask_array = np.zeros(length, bool)
for s, e in zip(mask_bed.s, mask_bed.e):
    mask_array[max(s - start, 0):max(min(e - start, length), 0)] = True
bases = np.frombuffer(reference[start:end].encode(), np.uint8)
tract = np.zeros(length, bool)
for k in range(2, 7):
    equal = np.concatenate([bases[:-k] == bases[k:], np.zeros(k, bool)])
    edges = np.flatnonzero(np.diff(np.concatenate([[0], equal.astype(np.int8), [0]])))
    for a, b in zip(edges[0::2], edges[1::2]):
        if b - a >= 2 * k:
            unit = reference[start + a:start + a + k]
            primitive = all(unit != unit[d:] + unit[:d] for d in range(1, k) if k % d == 0)
            if len(set(unit)) > 1 and primitive:
                tract[a:b + k] = True
holes = tract & ~mask_array
touching = holes[np.clip(calls.position.values - start, 0, length - 1)]
print(f"Exact repeat bases: {tract.sum():,}; not masked: {holes.sum():,} ({holes.sum() / tract.sum():.1%})")
print(f"Calls whose first affected base is in an unmasked exact repeat: {touching.sum()} "
      f"(kmerQ >= 93: {touching[calls.kmerQ.values >= 93].sum()})")

# ---------------------------------------------------------------------------
section("D1, D4: the 1 bp opportunity rule and the kmerQ windows, on the reads")
# ---------------------------------------------------------------------------
three = four = positions = differs = 0
for k, r in enumerate(bam.fetch(contig, start + 1_000_000, end - 1_000_000)):
    if k >= 2000:
        break
    quality = np.array(r.query_qualities)
    sequence = np.frombuffer(r.query_sequence.encode(), np.uint8)
    top = quality == 93
    three += (top[:-2] & top[1:-1] & top[2:]).sum()
    four += (top[:-3] & top[1:-2] & top[2:-1] & top[3:]).sum()
    positions += len(quality) - 3
    differs += ((sequence[1:-1] != sequence[:-2]) & (sequence[1:-1] != sequence[2:])).sum()
print(f"D1: read bases that differ from both neighbours (the 1 bp opportunity rule): {differs / positions:.1%}; "
      f"1 bp insertion rates are about {positions / differs:.1f}x too high")
print(f"D4: all-Q93 windows: three bases {three / positions:.1%}, four bases {four / positions:.1%}; "
      f"1 bp deletion rates are low by about {1 - four / three:.0%}")

# ---------------------------------------------------------------------------
section("D3, F1, N1")
# ---------------------------------------------------------------------------
opportunities = pd.read_csv(f"{CLUSTER}/opportunities.tsv", sep="\t", keep_default_na=False)
read_ends = sum(1 for r in bam.fetch(contig, start, end)
                if not (r.is_secondary or r.is_supplementary) and start <= r.reference_end - 1 < end)
print(f"D3: at most {3 * read_ends:,} truncated k-mers against {opportunities['count'].sum():,} opportunities "
      f"({3 * read_ends / opportunities['count'].sum():.1e})")
print("F1: alignment flags of the records carrying calls:", tags.flag.value_counts().to_dict())
q93_calls = tags[tags.kmerQ >= 93]
print("N1: diploid mapQ of the kmerQ >= 93 calls:",
      pd.cut(q93_calls.mapq, [-1, 0, 1, 30, 59, 60]).value_counts(sort=False).to_dict())
one_bp = opportunities[(opportunities.kmerQ == 93) & (opportunities.context.str.len() == 1)]
bins = [-1, 1, 59, 60]
unique_one_bp = q93_calls[(q93_calls.unique == "unique") & (q93_calls["size"] == 1)]
print(pd.DataFrame({"opportunities": one_bp.groupby(pd.cut(one_bp.mapQ, bins), observed=False)["count"].sum(),
                    "unique 1 bp calls": unique_one_bp.groupby(pd.cut(unique_one_bp.mapq, bins), observed=False).size()})
      .to_string())

# ---------------------------------------------------------------------------
section("The yardstick: Q93 read positions passing the subread filter, read correctly")
# ---------------------------------------------------------------------------
rng = np.random.default_rng(1)
sampled = passing = 0
for k, r in enumerate(bam.fetch(contig, start + 1_000_000, end - 1_000_000)):
    if k >= 1500 or r.is_secondary or r.is_supplementary:
        continue
    length_read, quality = r.query_length, np.array(r.query_qualities)
    sm, sx, sa = np.array(r.get_tag("sm")), np.array(r.get_tag("sx")), decode(r.get_tag("sa"))
    own = (lambda j: length_read - 1 - j) if r.is_reverse else (lambda j: j)
    pairs = r.get_aligned_pairs(matches_only=True)
    for qi, ri in (pairs[j] for j in rng.choice(len(pairs), 50, replace=False)):
        if not 2 <= qi <= length_read - 3 or masked(ri) or min(quality[qi - 1:qi + 2]) < 93:
            continue
        sampled += 1
        anchor, after = own(qi - 1), own(qi)
        passing += sx[anchor] == 0 and sx[after] == 0 and sm[after] == sa[after]
print(f"{sampled:,} Q93 read positions sampled; {passing / sampled:.1%} pass "
      "(no subread disagreeing at the anchor or the base after it, every subread matching the latter)")
