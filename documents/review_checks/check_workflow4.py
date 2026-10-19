"""
The genome-wide numbers behind ../workflow4_review.md, from the tables in
results/ (workflow 4's final calls, opportunities and filter reports, and
workflow 2's tract summaries).

    .venv/bin/python documents/review_checks/check_workflow4.py

Medians are over samples. A sample counts for a statistic only if it has calls.
"""
import glob
import os

import numpy as np
import pandas as pd

from common import NONREP, REPEAT, read_calls

MODELS = ["M0", "M1", "M2"]


def section(title):
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")


def label(path, depth):
    return "/".join(path.split(os.sep)[-depth - 1:-1])


# ---------------------------------------------------------------------------
section("1. The haploid check: share of calls it removes (filter_report.tsv)")
# ---------------------------------------------------------------------------
def removed_by_haploid_check(pattern, depth):
    shares = {}
    for path in glob.glob(pattern):
        report = pd.read_csv(path, sep="\t").set_index("filter")
        before = report.loc["haploid_check", "rows_before"]
        if before:
            shares[label(path, depth)] = 1 - report.loc["haploid_check", "rows_after"] / before
    return pd.Series(shares)


outside = removed_by_haploid_check(f"{NONREP}/*/*/filter_report.tsv", 2)
print(f"workflow 3:    median {outside.median():.1%}  (range {outside.min():.1%}-{outside.max():.1%}, {len(outside)} samples)")
for model in MODELS:
    inside = removed_by_haploid_check(f"{REPEAT}/*/*/{model}/filter_report.tsv", 3)
    print(f"workflow 4 {model}: median {inside.median():.1%}  (range {inside.min():.1%}-{inside.max():.1%}, {len(inside)} samples)")


# ---------------------------------------------------------------------------
section("2. Calls on flipped reads (haploid strand differs) among the final calls")
# ---------------------------------------------------------------------------
# Whether a read is flipped does not depend on the InDel it carries, so
# without a strand-dependent loss every class would show the same share.
# 1 bp InDels outside repeats cannot slide (homopolymer rule), so they are
# the reference.
rows = []
for path in glob.glob(f"{NONREP}/*/*/indels.tsv"):
    calls = read_calls(path, [])
    one = calls[calls["size"] == 1]
    if len(one):
        rows.append(("outside, 1 bp", label(path, 2), one.flip.mean(), len(one)))
for model in MODELS:
    for path in glob.glob(f"{REPEAT}/*/*/{model}/indels.tsv"):
        calls = read_calls(path, ["k"])
        for k, group in calls.groupby("k"):
            rows.append((f"{model} k={k}", label(path, 3), group.flip.mean(), len(group)))
flips = pd.DataFrame(rows, columns=["class", "sample", "flipped", "calls"])
flips = flips[flips.calls >= 20]
print(flips.groupby("class").flipped.agg(["median", "min", "max", "count"]).round(3).to_string())


# ---------------------------------------------------------------------------
section("3. One-unit InDels per Q93 read base, by model and period")
# ---------------------------------------------------------------------------
# Calls whose InDel is exactly one copy of the tract's motif (any rotation, any
# strand), over the opportunities of the same motif and period.
from utilities import canonical_repeat_indel  # noqa: E402

rates = []
for model in MODELS:
    for path in glob.glob(f"{REPEAT}/*/*/{model}/indels.tsv"):
        sample = label(path, 2 + 1).rsplit("/", 1)[0]   # species/sample, without the model
        if not os.path.exists(path.replace("indels.tsv", "opportunities.tsv")):
            continue
        calls = pd.read_csv(path, sep="\t", usecols=["indel_string", "canonical_motif", "k", "unique"],
                            low_memory=False)
        calls = calls[calls.unique == "unique"]
        one_unit = calls[calls.indel_string.str.len() == calls.k]
        one_unit = one_unit[np.array([canonical_repeat_indel(s) == m for s, m in zip(one_unit.indel_string,
                                                                                    one_unit.canonical_motif)],
                                     dtype=bool)]
        opportunities = pd.read_csv(path.replace("indels.tsv", "opportunities.tsv"), sep="\t",
                                    usecols=["k", "count"])
        per_k_calls = one_unit.groupby("k").size()
        per_k_opps = opportunities.groupby("k")["count"].sum()
        for k in per_k_opps.index:
            if per_k_opps[k] > 0:
                rates.append((model, k, sample, per_k_calls.get(k, 0) / per_k_opps[k]))
rates = pd.DataFrame(rates, columns=["model", "k", "sample", "rate"])
table = rates.groupby(["k", "model"]).rate.median().unstack()
print(table.to_string(float_format=lambda x: f"{x:.2e}"))
print("\nM1 / M0, median over samples, per k:")
paired = rates.pivot_table(index=["sample", "k"], columns="model", values="rate")
paired = paired[(paired.M0 > 0) & (paired.M1 > 0)]
print((paired.M1 / paired.M0).groupby(level="k").median().round(1).to_string())


# ---------------------------------------------------------------------------
section("4. M1 calls by tract purity, against the tract sequence in each purity class")
# ---------------------------------------------------------------------------
# The opportunity table has no purity, so the tract summaries stand in: the
# bases of tract sequence per purity class (copies x k x tracts). Within a
# sample, calls per tract base in exactly pure tracts are compared with tracts
# 90% to < 100% pure. Depth cancels out of the ratio.
def purity_class(purity):
    return np.where(purity >= 100, "pure", np.where(purity >= 90, "90-99.9", "<90"))


ratios = []
pooled = {}   # (copies bin, purity class) -> [calls, tract bases], summed over samples
COPY_BINS = [0, 15, 20, 30]
for path in glob.glob(f"{REPEAT}/*/*/M1/indels.tsv"):
    species, sample = label(path, 3).split("/")[:2]
    summary_path = f"results/repetitive_dna/{species}/{sample}/M1_summary.tsv"
    if not os.path.exists(summary_path):
        continue
    calls = pd.read_csv(path, sep="\t", usecols=["indel_string", "k", "copies", "trf_purity", "unique"],
                        low_memory=False)
    calls = calls[(calls.unique == "unique") & (calls.k == 2) & (calls.indel_string.str.len() == 2)]
    summary = pd.read_csv(summary_path, sep="\t")
    summary = summary[summary.k == 2]
    bases = (summary.copies * summary.k * summary["count"]).groupby(purity_class(summary.trf_purity)).sum()
    counts = calls.groupby(purity_class(calls.trf_purity)).size()
    if counts.get("pure", 0) and counts.get("90-99.9", 0):
        ratios.append((counts["pure"] / bases["pure"]) / (counts["90-99.9"] / bases["90-99.9"]))
    # The same, per bin of copy number, pooled: pure calls are too few per sample.
    for low, high in zip(COPY_BINS, COPY_BINS[1:]):
        tracts = summary[(summary.copies > low) & (summary.copies <= high)]
        called = calls[(calls.copies > low) & (calls.copies <= high)]
        tract_bases = (tracts.copies * tracts.k * tracts["count"]).groupby(purity_class(tracts.trf_purity)).sum()
        call_counts = called.groupby(purity_class(called.trf_purity)).size()
        for purity in ("pure", "90-99.9"):
            totals = pooled.setdefault((f"{low}-{high} copies", purity), [0, 0])
            totals[0] += call_counts.get(purity, 0)
            totals[1] += tract_bases.get(purity, 0)
print(f"k = 2, one-unit calls per tract base, pure / 90-99.9% pure: median {np.median(ratios):.2f} "
      f"(range {min(ratios):.2f}-{max(ratios):.2f}, {len(ratios)} samples)")
print("Within bins of tract length, pooled over samples:")
for low, high in zip(COPY_BINS, COPY_BINS[1:]):
    name = f"{low}-{high} copies"
    (pure_calls, pure_bases), (impure_calls, impure_bases) = pooled[(name, "pure")], pooled[(name, "90-99.9")]
    print(f"  {name:13s} pure {pure_calls:5d} calls / {pure_bases:12,.0f} bp, 90-99.9% {impure_calls:5d} calls / "
          f"{impure_bases:12,.0f} bp: ratio {(pure_calls / pure_bases) / (impure_calls / impure_bases):.2f}")


# ---------------------------------------------------------------------------
section("5. The read-end cut: calls only (filter_report.tsv), and read lengths")
# ---------------------------------------------------------------------------
for model in MODELS:
    shares, lengths = [], []
    for path in glob.glob(f"{REPEAT}/*/*/{model}/filter_report.tsv"):
        report = pd.read_csv(path, sep="\t").set_index("filter")
        before = report.loc["min_read_end_distance", "rows_before"]
        if before:
            shares.append(report.loc["min_read_end_distance", "rows_removed"] / before)
            calls = pd.read_csv(path.replace("filter_report.tsv", "indels.tsv"), sep="\t", usecols=["read_length"])
            lengths.append(calls.read_length.median())
    print(f"{model}: calls removed median {np.median(shares):.2%} (max {max(shares):.2%}); "
          f"median read length {np.median(lengths):,.0f}, so about {2 * 100 / np.median(lengths):.1%} "
          f"of read bases lie within 100 bp of an end")


# ---------------------------------------------------------------------------
section("6. Calls of more than one unit (their kmerQ window is not the opportunities')")
# ---------------------------------------------------------------------------
for model in MODELS:
    shares = []
    for path in glob.glob(f"{REPEAT}/*/*/{model}/indels.tsv"):
        calls = pd.read_csv(path, sep="\t", usecols=["indel_string", "k"])
        if len(calls):
            shares.append((calls.indel_string.str.len() > calls.k).mean())
    print(f"{model}: median {np.median(shares):.1%} of final calls are 2 or more units")


# ---------------------------------------------------------------------------
section("7. S1 in the repeat scan: full subread support by strand, unique M1 calls")
# ---------------------------------------------------------------------------
# sm == sa and sx == 0 at the tested base. The tags are only read correctly on
# forward-strand reads.
rows = []
for path in glob.glob(f"{REPEAT}/*/*/M1/indels.tsv"):
    calls = pd.read_csv(path, sep="\t", usecols=["unique", "sm", "sx", "sa", "diploid_direction"], low_memory=False)
    calls = calls[(calls.unique == "unique")].dropna(subset=["sm", "sx", "sa"])
    if len(calls) < 50:
        continue
    full = (calls.sm == calls.sa) & (calls.sx == 0)
    rows.append((full[calls.diploid_direction == "+"].mean(), full[calls.diploid_direction == "-"].mean()))
rows = np.array(rows)
print(f"full support, median over {len(rows)} samples: forward {np.median(rows[:, 0]):.1%}, "
      f"reverse {np.median(rows[:, 1]):.1%}")
