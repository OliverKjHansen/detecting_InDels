import os
import sys

from parse_trf import canonicalise

# Usage: python recanonicalise_motifs.py [--dry-run] <file> [<file> ...]
#
# Rewrites the canonical_motif column of repeat BEDs and {model}_summary.tsv
# files so a motif and its reverse complement share one label, matching the
# current canonicalise() in parse_trf.py.
#
# Files written before that change can be converted without re-running mreps or
# TRF: the stored motif already identifies the rotation class, and the collapsed
# form follows from it.
#
#   BED      (contig, start, end, canonical_motif, ...): column 4 is rewritten,
#            one line per tract, nothing merged.
#   summary  (canonical_motif, k, copies, ..., count): column 1 is rewritten and
#            rows that become identical are merged, adding their counts.
#
# Each file is rewritten in place, via a temporary file, so an interrupted run
# cannot leave a half-written table behind.

BED_HEADER = "contig"
SUMMARY_HEADER = "canonical_motif"
COUNT_HEADER = "count"


def rewrite_bed(rows):
    """Collapse the motif column of a repeat BED, keeping one row per tract."""
    for row in rows:
        row[3] = canonicalise(row[3])
    return rows


def rewrite_summary(rows):
    """Collapse the motif column of a summary and add up the merged counts."""
    merged = {}
    for row in rows:
        row[0] = canonicalise(row[0])
        key = tuple(row[:-1])
        merged[key] = merged.get(key, 0) + int(row[-1])

    # Same order as the workflow writes: period, then motif, then copies.
    ordered = sorted(merged.items(), key=lambda kv: (int(kv[0][1]), kv[0][0], float(kv[0][2])))
    return [list(key) + [str(count)] for key, count in ordered]


def rewrite(path, dry_run=False):
    with open(path) as handle:
        header = handle.readline()
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]

    columns = [c.strip().lstrip("#") for c in header.rstrip("\n").split("\t")]

    # Both the repeat BED and the InDel opportunity table start with "contig",
    # so the motif column has to be identified by name, not by position: in the
    # opportunity table column 4 is "tissue", which must not be touched.
    if len(columns) > 3 and columns[0] == BED_HEADER and columns[3] == SUMMARY_HEADER:
        motif_column, new_rows = 3, rewrite_bed(rows)
    elif columns[0] == SUMMARY_HEADER and columns[-1] == COUNT_HEADER:
        motif_column, new_rows = 0, rewrite_summary(rows)
    else:
        print(f"{path}: not a repeat BED or repeat summary (columns: "
              f"{', '.join(columns[:4])}...), skipped", file=sys.stderr)
        return

    motifs_after = len({row[motif_column] for row in new_rows})
    print(f"{path}: {len(rows)} rows -> {len(new_rows)} rows, "
          f"{motifs_after} motifs" + (" (dry run, nothing written)" if dry_run else ""))

    if dry_run:
        return

    temporary = path + ".tmp"
    with open(temporary, "w") as handle:
        handle.write(header)
        for row in new_rows:
            handle.write("\t".join(row) + "\n")
    os.replace(temporary, path)


def main():
    arguments = sys.argv[1:]
    dry_run = "--dry-run" in arguments
    paths = [a for a in arguments if a != "--dry-run"]

    if not paths:
        sys.exit("Usage: recanonicalise_motifs.py [--dry-run] <file> [<file> ...]")

    for path in paths:
        rewrite(path, dry_run)


if __name__ == "__main__":
    main()
