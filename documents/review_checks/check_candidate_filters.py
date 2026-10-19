"""
../candidate_filters.md and finding S1 of ../nonrepetitive_indel_review.md:
which properties of a call mark artefacts, measured against the real-sequence
yardstick (see yardstick_calls() in common.py).

    .venv/bin/python documents/review_checks/check_candidate_filters.py
"""
import math

import numpy as np
import pandas as pd

from common import NONREP, filtered_rates, load_calls, relative_rate, show, yardstick_calls

calls, opportunities = yardstick_calls(["ec", "sm", "sx", "sa", "sx_pre", "map_mapq", "softclip_haploid",
                                        "other_ops", "read_length"])
for column in ("ec", "sm", "sx", "sa", "sx_pre", "map_mapq"):
    calls[column] = pd.to_numeric(calls[column], errors="coerce")
calls["support"] = calls.sm / calls.sa
calls["haploid mapQ"] = np.where(calls.map_mapq >= 60, "60", "<60")
calls["haploid soft clip"] = calls.softclip_haploid.astype(str)
calls["diploid strand"] = calls.diploid_direction
calls["haplotype"] = calls.contig.str[:2]
calls["H/P/N operations"] = calls.other_ops.fillna("none")
revio = calls[calls.instrument == "Revio"]


def relative(feature, bins=None, labels=None, data=calls):
    print(f"Rate of unique calls by {feature}, relative to the sample's average")
    show(relative_rate(data, feature, bins, labels), 2)


# ---------------------------------------------------------------------------
# 1. Properties of the call and its read not covered in check_rate_floor.py.
# ---------------------------------------------------------------------------
relative("ec", [0, 8, 12, 16, 20, 30, 1000], ["<8", "8-12", "12-16", "16-20", "20-30", ">=30"])
relative("sx_pre", [0, 1, 2, 3, 1000], ["0", "1", "2", ">=3"])
for feature in ("haploid mapQ", "haploid soft clip", "diploid strand", "haplotype", "H/P/N operations"):
    relative(feature)

unique_calls = calls[calls.set == "unique"]
print("Share of unique calls a filter would remove, median over samples:")
show(pd.DataFrame({
    "read < 12 kb": unique_calls.groupby("label").read_length.apply(lambda s: (s < 12000).mean()),
    "< 1 kb from a read end": unique_calls.groupby("label").ToreadEnd.apply(lambda s: (s < 1000).mean()),
    "haploid mapQ < 60": unique_calls.groupby("label").map_mapq.apply(lambda s: (s < 60).mean()),
    "haploid soft clip": unique_calls.groupby("label")["haploid soft clip"].apply(lambda s: (s == "True").mean()),
}).median().to_frame("median"), 2)

# ---------------------------------------------------------------------------
# 2. S1: the subread tags on forward- and reverse-strand reads (Revio samples).
#    The controls set the baseline; if the tags were read at the right base on
#    both strands, the unique calls would differ from it equally on both.
# ---------------------------------------------------------------------------
print("S1. Share of calls with every subread agreeing (sm == sa) at the base after the anchor")
show(revio.dropna(subset=["support"]).groupby(["diploid strand", "set"]).support
     .apply(lambda s: (s == 1).mean()).unstack(), 3)
print("S1. Share of calls with two or more disagreeing subreads (sx >= 2)")
show(revio.dropna(subset=["sx"]).groupby(["diploid strand", "set"]).sx
     .apply(lambda s: (s >= 2).mean()).unstack(), 3)
for strand in "+-":
    relative("support", [0, .8, .9, .95, .9999, 1.01], ["<0.8", "0.8-0.9", "0.9-0.95", "0.95-<1", "1"],
             data=revio[revio["diploid strand"] == strand].assign(**{"support": lambda d: d.support}))
    print(f"   (above: {strand} strand only)\n")

# ---------------------------------------------------------------------------
# 3. What a working subread filter could do: forward-strand calls only, where
#    the tags are read at the right base, sensitivity-corrected per class.
# ---------------------------------------------------------------------------
forward = lambda d: d[d["diploid strand"] == "+"]
agree = lambda d: forward(d)[(forward(d).sx == 0) & (forward(d).sm == forward(d).sa) & (forward(d).sx_pre == 0)]
filters = {
    "none": lambda d: d,
    "+ strand": forward,
    "+ strand, subreads agree": agree,
    "+ strand, subreads agree, >=12 kb, >=1 kb from end":
        lambda d: agree(d)[(agree(d).read_length >= 12000) & (agree(d).ToreadEnd >= 1000)],
}
per_filter = {name: filtered_rates(revio, opportunities, keep) for name, keep in filters.items()}
print("Rate per class, median over Revio samples, under each filter")
show(pd.DataFrame({name: t.groupby("class").rate.median() for name, t in per_filter.items()}), 2)
summed = pd.DataFrame({name: t.groupby("label").rate.sum() for name, t in per_filter.items()})
print("Summed over the six classes, median over samples:")
show(summed.median().to_frame("rate"), 2)
print("Share of controls each filter keeps (median over samples and classes):")
show(pd.DataFrame({name: [t.kept.median()] for name, t in per_filter.items()}, index=["kept"]).T, 2)
strictest = list(filters)[-1]
print(f"Lowest samples under the strictest filter: {summed[strictest].nsmallest(3).to_dict()}\n")

# ---------------------------------------------------------------------------
# 4. Clustering: another unique call nearby, against the share expected if the
#    calls were spread at random over a diploid assembly of about 6 Gb.
# ---------------------------------------------------------------------------
everything = load_calls(["contig", "position", "unique", "readID", "read_length"])
rows = []
for label, u in calls[calls.set == "unique"].groupby("label"):
    if len(u) < 50:
        continue
    u = u.sort_values(["contig", "position"])
    gap_before = u.groupby("contig").position.diff()
    gap_after = -u.groupby("contig").position.diff(-1)
    nearest = pd.concat([gap_before, gap_after], axis=1).min(axis=1)
    rows.append({"sample": label, "unique calls": len(u),
                 "another within 1 kb": (nearest <= 1000).mean(),
                 "another within 100 bp": (nearest <= 100).mean(),
                 "expected within 1 kb": 1 - math.exp(-2 * 1000 * len(u) / 6e9)})
clusters = pd.DataFrame(rows).set_index("sample")
print("Unique calls with another unique call nearby, median over samples")
show(clusters.median().to_frame("median"), 3)

# ---------------------------------------------------------------------------
# 5. Contigs with more unique calls than chance allows (the notebook's test:
#    one-sided Poisson against the sample's rate, p < 0.001 / contigs), and
#    the rate once they are dropped with their opportunities.
# ---------------------------------------------------------------------------
def poisson_at_least(k, mu):
    """P(X >= k) for X ~ Poisson(mu)."""
    if k <= 0:
        return 1.0
    term = math.exp(-mu + k * math.log(mu) - math.lgamma(k + 1)) if mu > 0 else 0.0
    total, i = term, k
    while term > total * 1e-15 and i < k + 100_000:
        i += 1
        term *= mu / i
        total += term
    return min(1.0, total)


rows = []
for label, u in calls[calls.set == "unique"].groupby("label"):
    per_contig = pd.read_csv(f"{NONREP}/{label}/opportunities.tsv", sep="\t",
                             usecols=["contig", "count"]).groupby("contig")["count"].sum()
    t = pd.DataFrame({"positions": per_contig, "calls": u.groupby("contig").size()}).fillna(0)
    rate = t.calls.sum() / t.positions.sum()
    flagged = np.array([poisson_at_least(int(k), rate * p) for k, p in zip(t.calls, t.positions)]) < 0.001 / len(t)
    rows.append({"sample": label, "contigs flagged": flagged.sum(),
                 "share of calls": t.calls[flagged].sum() / t.calls.sum(),
                 "share of positions": t.positions[flagged].sum() / t.positions.sum(),
                 "rate after dropping": (t.calls[~flagged].sum() / t.positions[~flagged].sum()) / rate})
contigs = pd.DataFrame(rows).set_index("sample")
print("Contigs with an excess of unique calls, median over samples")
show(contigs.median().to_frame("median"), 3)

# ---------------------------------------------------------------------------
# 6. Supplementary alignments: a read whose calls show two read lengths had a
#    hard-clipped second alignment.
# ---------------------------------------------------------------------------
lengths = everything.groupby(["label", "readID"]).read_length.nunique()
print(f"Reads whose calls show more than one read_length: {(lengths > 1).sum()} of {len(lengths):,}")
