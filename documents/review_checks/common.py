"""
Loading code shared by the checks behind ../nonrepetitive_indel_review.md.

Run every check from the repository root with the project's virtualenv, e.g.

    .venv/bin/python documents/review_checks/check_haploid_check.py

The checks read only results/, the small final tables copied from the cluster.
"""
import glob
import os
import sys
import warnings

import numpy as np
import pandas as pd

# pandas warns when it stacks a sample whose table has an all-NA column.
warnings.simplefilter("ignore", FutureWarning)

NONREP = "results/nonrepetitive_indels"
REPEAT = "results/repetitive_indels"

# scripts/ holds the pipeline's own functions; the checks call them directly.
sys.path.insert(0, "scripts")

if not os.path.isdir(NONREP):
    sys.exit(f"No {NONREP}/ here. Run the checks from the repository root.")


def labels(root=NONREP):
    """species/sample for every sample directory under root."""
    return sorted(os.path.relpath(d, root) for d in glob.glob(f"{root}/*/*") if os.path.isdir(d))


def read_calls(path, columns):
    """One indels.tsv: the given columns plus size, flip and label."""
    wanted = set(columns) | {"indel_string", "diploid_direction", "bitflag_hap"}
    calls = pd.read_csv(path, sep="\t", low_memory=False, usecols=lambda c: c in wanted)
    calls["size"] = calls.indel_string.str.len()
    # Only reads aligned as flag 0 or 16 survive the haploid check, so 16 is
    # the reverse strand and everything else the forward one.
    haploid_direction = np.where(calls.bitflag_hap == 16, "-", "+")
    calls["flip"] = haploid_direction != calls.diploid_direction
    # .../{species}/{sample}/indels.tsv
    calls["label"] = "/".join(path.split(os.sep)[-3:-1])
    return calls


def load_calls(columns, root=NONREP):
    """indels.tsv of every sample, stacked."""
    paths = sorted(glob.glob(f"{root}/*/*/indels.tsv"))
    return pd.concat([read_calls(p, columns) for p in paths], ignore_index=True)


def comparison_view(calls):
    """
    The calls compared throughout the review: unique, matched at haploid
    distance 0 (the notebook's defaults), and at least 100 bp from a read end
    (the filter workflow 4 applies), so calls inside and outside repeats are
    filtered alike.
    """
    return calls[(calls.unique == "unique")
                 & (calls.distance_to_indel_in_hap == 0)
                 & (calls.ToreadEnd >= 100)]


def opportunities(label, root=NONREP):
    """A sample's Q93 opportunities summed per context, or None if absent."""
    path = f"{root}/{label}/opportunities.tsv"
    if not os.path.exists(path):
        return None
    # keep_default_na=False: a context is a k-mer, never missing.
    table = pd.read_csv(path, sep="\t", usecols=["context", "count"], keep_default_na=False)
    table = table.groupby("context", as_index=False)["count"].sum()
    table["size"] = table.context.str.len()
    return table


def opportunities_by_size(label, root=NONREP):
    """A sample's Q93 opportunities per InDel size: a Series indexed 1, 2, 3."""
    table = opportunities(label, root)
    if table is None:
        return None
    return table.groupby("size")["count"].sum().reindex([1, 2, 3], fill_value=0)


def show(frame, digits=3):
    """Print a table with numbers to a few significant digits."""
    print(frame.to_string(float_format=lambda x: f"{x:.{digits}g}"))
    print()


# ---------------------------------------------------------------------------
# The real-sequence yardstick, used by check_rate_floor.py and
# check_candidate_filters.py.
#
# The opportunity table cannot be split by read length, passes, subread support
# or any other property of a read, so unique calls are compared with control
# calls that are real sequence in the read: non-unique calls whose InDel sits at
# the same site in at least three reads (germline differences or assembly
# errors). Real sequence is there whatever the read's quality, so the controls
# are spread over reads the way Q93 positions are, and stand in for the
# opportunities.
# ---------------------------------------------------------------------------

CLASSES = [(kind, size) for size in (1, 2, 3) for kind in ("deletion", "insertion")]


def yardstick_calls(extra_columns=()):
    """
    The unique calls and the controls of every sample with Q93 opportunities.
    Both sets are matched at haploid distance 0 and lie at least 100 bp from a
    read end. Returns the calls, with `set` ("unique" or "control"), `cell`
    (SMRT Cell) and `instrument` added, and each sample's opportunities per
    InDel size.
    """
    columns = {"contig", "position", "indel_type", "unique", "distance_to_indel_in_hap",
               "ToreadEnd", "readID", *extra_columns}
    parts, by_size = [], {}
    for path in sorted(glob.glob(f"{NONREP}/*/*/indels.tsv")):
        label = "/".join(path.split(os.sep)[-3:-1])
        sizes = opportunities_by_size(label)
        calls = read_calls(path, columns)
        if sizes is None or not sizes.sum() or not len(calls):
            continue
        by_size[label] = sizes
        calls = calls[calls.distance_to_indel_in_hap == 0].copy()
        site = ["contig", "position", "indel_type", "indel_string"]
        calls["reads_at_site"] = calls.groupby(site).readID.transform("size")
        parts.append(calls)
    calls = pd.concat(parts, ignore_index=True)
    calls = calls[calls.ToreadEnd >= 100].copy()
    calls["set"] = np.select([calls.unique == "unique",
                              (calls.unique == "non-unique") & (calls.reads_at_site >= 3)],
                             ["unique", "control"], "")
    calls = calls[calls.set != ""].copy()
    calls["cell"] = calls.readID.str.split("/").str[0]
    calls["instrument"] = np.where(calls.cell.str.startswith("m84"), "Revio", "Sequel II")
    return calls, by_size


def relative_rate(calls, feature, bins=None, labels=None):
    """
    Rate of unique calls in each group of reads relative to the sample's
    average: within a sample, the share of its unique calls in the group over
    the share of its controls in it. Median and quartiles over samples with at
    least 200 controls, counting a group only where it holds 30 controls.
    """
    d = calls.dropna(subset=[feature]).copy()
    d["bin"] = pd.cut(d[feature], bins, labels=labels, right=False) if bins is not None else d[feature]
    t = d.groupby(["label", "bin", "set"], observed=True).size().unstack(fill_value=0)
    t = t[t.groupby("label").control.transform("sum") >= 200]
    share = t / t.groupby("label").transform("sum")
    ratio = (share.unique / share.control).where(t.control >= 30).unstack()
    return pd.DataFrame({"median": ratio.median(), "q25": ratio.quantile(.25), "q75": ratio.quantile(.75),
                         "share of unique calls": share.unique.unstack().median(),
                         "samples": ratio.notna().sum()})


def filtered_rates(calls, by_size, keep):
    """
    Rate per sample and class (size and type) after a filter on the calls,
    divided by the share of the sample's controls the same filter keeps, which
    estimates the share of real mutations it keeps. A class's own controls are
    the yardstick where the sample has at least 100 of them, all of its controls
    otherwise. A class without calls counts as a rate of 0.
    """
    rows = []
    for label, sample in calls.groupby("label"):
        controls = sample[sample.set == "control"]
        unique = keep(sample[sample.set == "unique"])
        for kind, size in CLASSES:
            own = controls[(controls.indel_type == kind) & (controls["size"] == size)]
            yardstick = own if len(own) >= 100 else controls
            kept = len(keep(yardstick)) / len(yardstick) if len(yardstick) else np.nan
            n = int(((unique.indel_type == kind) & (unique["size"] == size)).sum())
            positions = by_size[label][size]
            rows.append({"label": label, "class": f"{size} bp {kind}", "kept": kept, "calls": n,
                         "rate": n / (positions * kept) if positions and kept else np.nan})
    return pd.DataFrame(rows)
