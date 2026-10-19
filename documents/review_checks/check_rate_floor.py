"""
../rate_floor.md: how far filters on read-level features lower the rate
outside repeats.

    .venv/bin/python documents/review_checks/check_rate_floor.py

The opportunity table cannot be split by read length, passes or subread
support, so unique calls are compared with control calls that are real
sequence in the read (see yardstick_calls() in common.py).

The subread tags are read at the wrong base on reverse-strand reads (S1 in the
review), so the subread filter here works on about half the calls only;
check_candidate_filters.py estimates what it does when it works.
"""
import pandas as pd

from common import filtered_rates, relative_rate, show, yardstick_calls

SPERM = {"human/aa_sperm", "human/chk", "human/da1", "human/la", "human/la4", "human/ph", "human/tr", "human/tsaed"}

calls, opportunities = yardstick_calls(["np", "sm", "sx", "sa", "read_length"])
for column in ("np", "sm", "sx", "sa"):
    calls[column] = pd.to_numeric(calls[column], errors="coerce")
calls["support"] = calls.sm / calls.sa
revio = set(calls.loc[calls.instrument == "Revio", "label"])
print(f"{(calls.set == 'unique').sum():,} unique calls and {(calls.set == 'control').sum():,} controls "
      f"in {calls.label.nunique()} samples ({len(revio)} Revio)\n")


def relative(feature, bins, labels):
    print(f"Rate of unique calls by {feature}, relative to the sample's average")
    show(relative_rate(calls, feature, bins, labels), 2)


# 1. Read-level features.
relative("read_length", [0, 4000, 8000, 12000, 16000, 20000, 1e9],
         ["<4 kb", "4-8 kb", "8-12 kb", "12-16 kb", "16-20 kb", ">=20 kb"])
relative("np", [0, 10, 15, 20, 30, 1000], ["<10", "10-14", "15-19", "20-29", ">=30"])
relative("support", [0, .8, .9, .95, .9999, 1.01], ["<0.8", "0.8-0.9", "0.9-0.95", "0.95-<1", "1"])
relative("sx", [0, 1, 2, 3, 1000], ["0", "1", "2", ">=3"])

# 2. SMRT Cells: relative rate of each cell within its sample.
t = calls.groupby(["label", "cell", "set"]).size().unstack(fill_value=0)
t = t[t.control >= 100]
share = t / t.groupby("label").transform("sum")
ratio = (share.unique / share.control).groupby("label").agg(["size", "min", "max"])
ratio = ratio[ratio["size"] >= 2].rename(columns={"size": "cells", "min": "lowest", "max": "highest"})
print("Rate of unique calls per SMRT Cell, relative to the sample's average")
show(ratio, 2)

# 3. Instrument: the 1 bp rate (both types) of every sample with unique calls
#    (human/tsaed has none at all, see N5 in the review).
unique = calls[calls.set == "unique"]
one_bp = (unique[unique["size"] == 1].groupby("label").size()
          / pd.Series({k: v[1] for k, v in opportunities.items()})).dropna()
instrument = calls.groupby("label").instrument.agg(lambda s: s.mode()[0])
by_instrument = one_bp.groupby(instrument).agg(["size", "median"])
print("1 bp rate (unique, both types) by instrument")
show(by_instrument)
print(f"   Sequel II / Revio medians: {by_instrument.loc['Sequel II', 'median'] / by_instrument.loc['Revio', 'median']:.1f}x\n")

# 4. Which base the subread tags test: the base after the anchor, which is the
#    first inserted base of an insertion but a base flanking a deletion.
full = calls.dropna(subset=["support"]).groupby(["indel_type", "size", "set"]).support.apply(lambda s: (s == 1).mean())
print("Share of calls with every subread agreeing at the base after the anchor")
show(full.unstack(), 3)


# 5. The cleanest subset, Revio samples, per class, sensitivity-corrected with
#    the controls (see filtered_rates() in common.py).
def clean(d):
    return d[(d.read_length >= 12000) & (d.ToreadEnd >= 1000) & (d.sx == 0) & (d.sm == d.sa)]


r = calls[calls.label.isin(revio)]
before = filtered_rates(r, opportunities, lambda d: d)
after = filtered_rates(r, opportunities, clean)
floor = before[["label", "class"]].assign(before=before.rate, after=after.rate,
                                          kept=after.kept, **{"clean calls": after.calls})
print("Share of controls the clean filters keep, over samples and classes: "
      f"median {floor.kept.median():.0%}, range {floor.kept.min():.0%} to {floor.kept.max():.0%}; "
      f"lowest in {', '.join(floor.groupby('label').kept.median().nsmallest(2).index)}\n")
print("Rate per class, median over Revio samples, before and after the clean filters")
show(floor.groupby("class")[["before", "after"]].median())
summed = floor.groupby("label")[["before", "after"]].sum()
print(f"Summed over the six classes, median over samples: before {summed.before.median():.2g}, "
      f"after {summed.after.median():.2g} ({summed.before.median() / summed.after.median():.1f}x lower); "
      f"lowest sample after: {summed.after.min():.2g} ({summed.after.idxmin()})")
two_del = floor[floor["class"] == "2 bp deletion"].set_index("label").after
print(f"2 bp deletions as a share of what is left, median over samples: {(two_del / summed.after).median():.0%}\n")
sperm = floor[(floor["class"] == "1 bp deletion") & floor.label.isin(SPERM)].set_index("label")
print("Sperm samples, 1 bp deletions after the clean filters")
show(sperm[["clean calls", "after"]])
