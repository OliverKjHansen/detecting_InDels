"""
Finding L1 of ../nonrepetitive_indel_review.md: InDels in repeat tracts that
the aligner places at the base before the tract, outside workflow 3's mask.

    .venv/bin/python documents/review_checks/check_repeat_leak.py
"""
import os

import pandas as pd

from common import REPEAT, comparison_view, labels, load_calls, opportunities, show
from utilities import canonical_repeat_indel

calls = comparison_view(load_calls(["indel_type", "unique", "distance_to_indel_in_hap", "ToreadEnd",
                                    "canonical_indel"]))

# ---------------------------------------------------------------------------
# 1. The 2 bp rate per context outside repeats, all samples pooled.
# ---------------------------------------------------------------------------
pooled = pd.concat([o for o in map(opportunities, labels()) if o is not None])
pooled = pooled[pooled["size"] == 2].groupby("context")["count"].sum()
counts = calls[calls["size"] == 2].groupby(["canonical_indel", "indel_type"]).size().unstack(fill_value=0)
contexts = pd.DataFrame({"opportunities": pooled}).join(counts).fillna(0)
contexts["deletion rate"] = contexts.deletion / contexts.opportunities
contexts["insertion rate"] = contexts.insertion / contexts.opportunities
print("1. 2 bp InDel rate per context outside repeats, all samples pooled")
show(contexts.sort_values("deletion rate", ascending=False))

# ---------------------------------------------------------------------------
# 2. Outside repeats against inside dinucleotide tracts, per sample: the rate
#    of 2 bp InDels outside, and of one-unit InDels (a copy of the tract's
#    motif) inside M0 (exact) and M1 (impure) tracts with k = 2. The repeat
#    calls are filtered the same way as the calls outside.
# ---------------------------------------------------------------------------
def inside(label, model):
    folder = f"{REPEAT}/{label}/{model}"
    if not os.path.exists(f"{folder}/opportunities.tsv"):
        return float("nan"), float("nan")
    bases = pd.read_csv(f"{folder}/opportunities.tsv", sep="\t", usecols=["k", "count"])
    bases = bases.loc[bases.k == 2, "count"].sum()
    repeat_calls = pd.read_csv(f"{folder}/indels.tsv", sep="\t", low_memory=False,
                               usecols=["indel_string", "canonical_motif", "k", "unique",
                                        "distance_to_indel_in_hap", "ToreadEnd"])
    repeat_calls = comparison_view(repeat_calls)
    repeat_calls = repeat_calls[(repeat_calls.k == 2) & (repeat_calls.indel_string.str.len() == 2)]
    one_unit = sum(canonical_repeat_indel(s) == m
                   for s, m in zip(repeat_calls.indel_string, repeat_calls.canonical_motif))
    return (one_unit / bases if bases else float("nan")), bases


rows = []
for label in labels():
    table = opportunities(label)
    if table is None:
        continue
    outside_bases = table.loc[table["size"] == 2, "count"].sum()
    outside_calls = (calls[calls.label == label]["size"] == 2).sum()
    if not outside_bases or not outside_calls:
        continue
    m0_rate, m0_bases = inside(label, "M0")
    m1_rate, _ = inside(label, "M1")
    rows.append(dict(sample=label, outside=outside_calls / outside_bases, M0=m0_rate, M1=m1_rate,
                     outside_calls=outside_calls, m0_bases_vs_outside=m0_bases / outside_bases))
side = pd.DataFrame(rows).set_index("sample")
side["M0 / outside"] = side.M0 / side.outside
side["M1 / outside"] = side.M1 / side.outside
print("2. 2 bp InDels per read base: outside repeats, and one-unit InDels inside k = 2 tracts")
show(side[["outside", "M0", "M1", "M0 / outside", "M1 / outside"]])
print(f"   median M0 / outside: {side['M0 / outside'].median():.2f}; "
      f"median M1 / outside: {side['M1 / outside'].median():.1f}; "
      f"samples with M0 below outside: {(side['M0 / outside'] < 1).sum()} of {side['M0 / outside'].notna().sum()}\n")

# ---------------------------------------------------------------------------
# 3. How many of the outside 2 bp calls the events missing from M0 could be,
#    if exact k = 2 tracts mutated at a fraction f of the M1 per-base rate.
# ---------------------------------------------------------------------------
print("3. Share of outside 2 bp calls explained by events missing from M0 k = 2 tracts,")
print("   if those tracts mutated at a fraction f of the M1 rate (median over samples)")
print(f"   M0 k = 2 read bases relative to outside: {side.m0_bases_vs_outside.median():.3f}")
for f in (0.03, 0.1, 0.3):
    # missing events per outside read base, over outside calls per outside read base
    missing = (f * side.M1 - side.M0).clip(lower=0) * side.m0_bases_vs_outside
    print(f"   f = {f}: {(missing / side.outside).median():.0%}")
print()

# ---------------------------------------------------------------------------
# 4. How dense the mask is: tracts per model, and read bases inside each
#    model's tracts relative to read positions outside the mask.
# ---------------------------------------------------------------------------
print("4. Tracts per model, and read bases in them relative to outside (Q93)")
rows = []
for label in ("human/da1", "chimp/carl", "rhesus_macaque/R22013"):
    outside = opportunities(label)
    outside = outside.loc[outside["size"] == 3, "count"].sum()
    row = {"sample": label}
    for model in ("M0", "M1", "M2"):
        summary = pd.read_csv(f"results/repetitive_dna/{label}/{model}_summary.tsv", sep="\t")
        bases = pd.read_csv(f"{REPEAT}/{label}/{model}/opportunities.tsv", sep="\t", usecols=["count"])
        row[f"{model} tracts (M)"] = summary["count"].sum() / 1e6
        row[f"{model} bases / outside"] = bases["count"].sum() / outside
    rows.append(row)
show(pd.DataFrame(rows).set_index("sample"))
