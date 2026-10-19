"""
Finding E1 of ../nonrepetitive_indel_review.md: the rate per sample, and how
much of the call set the kmerQ >= 93 filter removes.

    .venv/bin/python documents/review_checks/check_error_floor.py
"""
import pandas as pd

from common import NONREP, comparison_view, labels, load_calls, opportunities_by_size, show

calls = comparison_view(load_calls(["unique", "distance_to_indel_in_hap", "ToreadEnd"]))

rows = []
for label in labels():
    by_size = opportunities_by_size(label)
    if by_size is None or not by_size.sum():
        continue
    report = pd.read_csv(f"{NONREP}/{label}/filter_report.tsv", sep="\t").set_index("filter")
    n = calls[calls.label == label]["size"].value_counts()
    rows.append({"sample": label,
                 "kept by kmerQ >= 93": report.loc["min_kmerQ", "rows_after"] / report.loc["min_kmerQ", "rows_before"],
                 **{f"{k} bp rate": n.get(k, 0) / by_size[k] for k in (1, 2, 3)}})
rates = pd.DataFrame(rows).set_index("sample").sort_values("1 bp rate")

print("InDels per read position outside repeats: unique, haploid distance 0, >= 100 bp from a read end")
show(rates)
with_calls = rates[rates["1 bp rate"] > 0]
print("Range over samples with calls:")
show(with_calls.describe().loc[["min", "50%", "max"]])
spread = with_calls.max() / with_calls.min()
print("max / min: " + ", ".join(f"{c} {spread[c]:.1f}x" for c in ("1 bp rate", "2 bp rate", "3 bp rate")))
human = with_calls[with_calls.index.str.startswith("human/")]["1 bp rate"]
print(f"max / min of the 1 bp rate within human: {human.max() / human.min():.1f}x")
print(f"samples whose 2 bp rate is above their 1 bp rate: "
      f"{(with_calls['2 bp rate'] > with_calls['1 bp rate']).sum()} of {len(with_calls)}")
