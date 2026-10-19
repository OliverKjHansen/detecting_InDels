"""
Findings H1-H4 of ../nonrepetitive_indel_review.md: how
scripts/filtering_haploid_mapping.py (then filtering_haploid_mappping.py) matched a call to its read's alignment
to the haploid assembly when the review was written, and, in section 1, how
it does since the fix. Sections 2-5 read results/, which the walk as reviewed
produced.

    .venv/bin/python documents/review_checks/check_haploid_check.py
"""
import pandas as pd

from common import load_calls, show
from filtering_haploid_mapping import check_call
from old_haploid_walk import find_indel_position_in_read

# ---------------------------------------------------------------------------
# 1. The CIGAR walk on reads whose haploid alignment carries the same InDel.
#
# A 40 bp read. In the diploid alignment (forward strand) the InDel follows
# read base 17, so the call's 1-based start is 17. The pipeline handed the
# walk `start` on the same strand and `read_length - start - 1` on a flipped
# one. Every CIGAR below holds exactly the called InDel, so the right answer
# is always (True, 0). The walk as reviewed (old_haploid_walk.py) is shown
# next to the pipeline's check as fixed since, which takes the read itself;
# nothing in this read can slide.
# ---------------------------------------------------------------------------
L, START = 40, 17
FLIPPED = L - START - 1
READ = "GATCCAGTGCATGACTCAGTCGATCAAGCTGGTACACGTA"
cases = [
    # what is tested,                                  haploid CIGAR,        flipped, type,        sequence
    ("same strand, 1 bp deletion",                      "17M1D23M",           False,   "deletion",  "G"),
    ("same strand, 1 bp deletion, 2 bp insertion before", "5M2I10M1D23M",     False,   "deletion",  "G"),
    ("  as above, plus a 1 bp deletion further on",     "5M2I10M1D10M1D13M",  False,   "deletion",  "G"),
    ("flipped strand, 1 bp deletion",                   "23M1D17M",           True,    "deletion",  "G"),
    ("  as above, 1 inserted base before it",           "10M1I12M1D17M",      True,    "deletion",  "G"),
    ("  as above, 2 inserted bases before it",          "10M2I11M1D17M",      True,    "deletion",  "G"),
    ("flipped strand, 1 bp insertion",                  "22M1I17M",           True,    "insertion", "A"),
    ("flipped strand, 2 bp insertion",                  "21M2I17M",           True,    "insertion", "AG"),
]
rows = []
for what, cigar, flipped, kind, sequence in cases:
    check, distance, _, _ = find_indel_position_in_read(cigar, FLIPPED if flipped else START, kind, sequence)
    fixed, fixed_distance, _, _ = check_call(READ, 16 if flipped else 0, cigar, START, kind, sequence, "+")
    rows.append((what, cigar, check, distance, fixed, fixed_distance))
print("1. The haploid check when the haploid alignment has the same InDel, as reviewed and as fixed")
print("   (correct answer: haploid_check True, distance 0)\n")
show(pd.DataFrame(rows, columns=["case", "haploid CIGAR", "haploid_check", "distance",
                                 "fixed: haploid_check", "fixed: distance"]).set_index("case"))

# ---------------------------------------------------------------------------
# 2-4. The final calls in results/ (haploid_check True and kmerQ >= 93).
# ---------------------------------------------------------------------------
calls = load_calls(["contig", "indel_type", "unique", "distance_to_indel_in_hap", "start", "read_length"])
print(f"{len(calls):,} calls; {calls.flip.mean():.1%} on reads whose strand differs between "
      "the diploid and the haploid alignment\n")

# 2. Flipping is a property of the read, so among correctly checked calls its
#    share should not depend on the InDel's size.
share = calls.groupby(["label", "size"]).flip.mean().unstack()
print("2. Share of calls on flipped reads, range over samples, by InDel size")
show(share.describe().loc[["min", "50%", "max"]].rename(index={"50%": "median"}), 2)
pooled = calls.groupby(["indel_type", calls["size"].map({1: "1 bp", 2: "2-3 bp", 3: "2-3 bp"})]).flip.mean()
print("   All samples pooled:")
show(pooled.unstack(), 3)

# 3. Where the matched haploid event is, by strand flip, type and size.
d = calls.distance_to_indel_in_hap
calls["dist"] = pd.cut(d, [-1, 0, 1, 10, float("inf")], labels=["0", "1", "2-10", ">10"])
table = calls.groupby(["flip", "indel_type", "size"]).dist.value_counts(normalize=True).unstack()
table["calls"] = calls.groupby(["flip", "indel_type", "size"]).size()
table["median distance"] = calls.groupby(["flip", "indel_type", "size"]).distance_to_indel_in_hap.median()
print("3. Distance to the matched haploid event (share of calls)")
show(table, 3)

dels = calls[calls.indel_type == "deletion"]
print("   Deletions kept by the notebook's max_hap_distance = 0: same strand "
      f"{(dels[~dels.flip].distance_to_indel_in_hap == 0).mean():.1%}, flipped "
      f"{(dels[dels.flip].distance_to_indel_in_hap == 0).mean():.1%}\n")

# 4. The walk skips the true event once an insertion precedes it in the read,
#    and those are more likely the further along the read the call is.
same = calls[~calls.flip].copy()
same["along the read"] = pd.cut(same.start / same.read_length, [0, .1, .25, .5, .75, .9, 1],
                                include_lowest=True)
trend = (same.groupby(["indel_type", "size", "along the read"], observed=True)
             .distance_to_indel_in_hap.apply(lambda x: (x == 0).mean()).unstack())
print("4. Same-strand calls matched at distance 0, by position of the call along the read")
show(trend, 2)

# 5. H2 and H3 genome-wide. Whether a read is flipped depends on its contig:
#    a contig whose haploid counterpart runs the other way flips every read on
#    it. If the check loses InDels on flipped reads, the rate on those contigs
#    falls, and the opportunities are per contig, so the comparison is exact.
#    Notebook view: unique calls matched at distance 0.
from common import NONREP  # noqa: E402

rows = []
for label in ("human/chk", "human/la", "human/ph", "human/tr", "chimp/yoran",
              "rhesus_macaque/R22013", "gorilla/kipenzi"):
    sample_calls = calls[calls.label == label]
    flipped_contig = sample_calls.groupby("contig").flip.mean() > 0.5
    table = pd.read_csv(f"{NONREP}/{label}/opportunities.tsv", sep="\t", usecols=["contig", "context", "count"],
                        keep_default_na=False)
    table["size"] = table.context.str.len()
    table["flipped"] = table.contig.map(flipped_contig).fillna(False).astype(bool)
    shown = sample_calls[(sample_calls.unique == "unique") & (sample_calls.distance_to_indel_in_hap == 0)].copy()
    shown["flipped"] = shown.contig.map(flipped_contig).astype(bool)
    rate = (shown.groupby(["flipped", "size"]).size()
            / table.groupby(["flipped", "size"])["count"].sum()).unstack()
    rows.append({"sample": label, **{f"{k} bp": rate.loc[True, k] / rate.loc[False, k] for k in (1, 2, 3)},
                 "flipped share of positions": table.loc[table.flipped, "count"].sum() / table["count"].sum()})
print("5. Rate on contigs whose haploid counterpart runs the other way, relative to the other contigs")
show(pd.DataFrame(rows).set_index("sample"), 2)
