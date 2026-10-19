"""
The denominator section and findings D1, D2 and N2 of
../nonrepetitive_indel_review.md: is the opportunity table the size it should
be?

    .venv/bin/python documents/review_checks/check_denominator.py
"""
import glob
import os

import pandas as pd

from common import NONREP, labels, opportunities_by_size, show

# CpG dinucleotides per base in primate genomes, about 28 million in 3.1 Gb.
CPG_PER_BASE = 0.009

# ---------------------------------------------------------------------------
# 1. Q93 opportunities against the bases sequenced, for the samples whose
#    read count is known (classification_reads.tsv counts every read in the
#    sample's methylation likelihood table). Bases are estimated two ways:
#    reads times the median length of the reads that carry a call (those run
#    a little long, so this overestimates the bases), and CpGs on all reads
#    over the genome's CpG density.
# ---------------------------------------------------------------------------
rows = []
for path in sorted(set(glob.glob("results/*/*/*/classification_reads.tsv"))):
    label = "/".join(path.split(os.sep)[-3:-1])
    if any(r["sample"] == label for r in rows):
        continue
    counts = pd.read_csv(path, sep="\t")
    by_size = opportunities_by_size(label)
    lengths = pd.read_csv(f"{NONREP}/{label}/indels.tsv", sep="\t", usecols=["readID", "read_length"])
    median_length = lengths.drop_duplicates("readID").read_length.median()
    bases_from_length = counts.reads.sum() * median_length
    bases_from_cpg = counts.CpGs.sum() / CPG_PER_BASE
    rows.append({"sample": label, "reads (M)": counts.reads.sum() / 1e6,
                 "median length": median_length,
                 "Q93 opportunities, 3 bp (G)": by_size[3] / 1e9,
                 "share of bases (length)": by_size[3] / bases_from_length,
                 "share of bases (CpG)": by_size[3] / bases_from_cpg})
print("1. Q93 opportunities of size 3 (nearly every read position) as a share of the bases sequenced")
show(pd.DataFrame(rows).set_index("sample"))

# ---------------------------------------------------------------------------
# 2-4. Every sample.
# ---------------------------------------------------------------------------
rows = []
for label in labels():
    path = f"{NONREP}/{label}/opportunities.tsv"
    if not os.path.exists(path):
        continue
    table = pd.read_csv(path, sep="\t", usecols=["context", "mapQ", "count"], keep_default_na=False)
    size = table.context.str.len()
    total = table["count"].sum()
    rows.append({"sample": label,
                 "1 bp / 3 bp opportunities": (table.loc[size == 1, "count"].sum()
                                               / table.loc[size == 3, "count"].sum()) if total else float("nan"),
                 "mapQ <= 1 share": table.loc[table.mapQ <= 1, "count"].sum() / total if total else float("nan"),
                 "Q93 opportunities": total})
every = pd.DataFrame(rows).set_index("sample")

print("2. 1 bp opportunities relative to size-3 ones (insertion sites are almost every position)")
show(every["1 bp / 3 bp opportunities"].describe()[["min", "50%", "max"]].to_frame().T)

print("3. Share of Q93 opportunities at diploid mapQ 0 or 1 (the calls carry no diploid mapQ)")
show(every["mapQ <= 1 share"].describe()[["min", "50%", "max"]].to_frame().T)

print("4. Samples with no Q93 opportunity at all, and what final filtering did to their calls")
for label in every.index[every["Q93 opportunities"] == 0]:
    report = pd.read_csv(f"{NONREP}/{label}/filter_report.tsv", sep="\t").set_index("filter")
    kmerq = report.loc["min_kmerQ"]
    print(f"   {label}: kmerQ >= 93 kept {kmerq.rows_after:,} of {kmerq.rows_before:,} calls")
