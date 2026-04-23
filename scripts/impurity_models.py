import sys
import bisect
import pysam
from collections import defaultdict

# Usage:
#   python impurity_models.py <trs_raw.tsv> <contig> <k> <assembly.fa> <output_prefix>
#
# Reads the per-contig, per-period raw TSV produced by parse_mreps_canonical.py and applies
# two biologically-grounded impurity models (M0, M1) that decide whether adjacent repeat
# tracts sharing the same canonical motif should be merged.
#
# Output: one TSV per model at {output_prefix}_M{0,1}_merged.tsv
# Columns: contig, start, end, canonical_motif, k, copies, num_merged, model, species, sample
#
# Models in the current framework:
#
#   M0 — exact baseline (this script): mreps pure/exact repeats, never merged.
#         One tract per mreps hit. This is the strictest possible annotation.
#
#   M1 — TRF standard (separate pipeline): Tandem Repeat Finder with standard
#         parameters (2 7 7 80 10 50 k), run per period k. Handled by the
#         trf_k{k} GWF targets and a dedicated parse step (not this script).
#
# Preserved future model — M_fission (not yet active in the workflow):
#   Merge two adjacent same-motif mreps tracts when d = 1, under either:
#     - Condition 1: point mutation (substitution) — gap base is a mismatch
#     - Condition 2: 1-bp indel — gap base matches motif but causes phase offset
#   See merge_tracts() below for the full implementation.

models = [0]  # only M0 is produced from this script; M1 = TRF (separate pipeline)
header = ("contig", "start", "end", "canonical_motif", "k",
          "copies", "num_merged", "model", "species", "sample")


def compute_mismatches(gap_seq, gap_start_pos, canonical_motif, k):
    """
    Count positions in gap_seq that deviate from the expected motif continuation.

    The expected character at gap position j is canonical_motif[(gap_start_pos + j) % k].
    gap_start_pos = end(Ri) % k, i.e. the motif-cycle offset at the first gap base.
    """
    mismatches = 0
    for j, base in enumerate(gap_seq):
        expected = canonical_motif[(gap_start_pos + j) % k]
        if base.upper() != expected.upper():
            mismatches += 1
    return mismatches


def gap_is_clear(gap_start, gap_end, all_starts_sorted):
    """
    Return True if no tract of any motif starts in the half-open interval
    [gap_start, gap_end).  Uses binary search on the globally sorted start list.

    Tracts of the current motif cannot appear in this interval during greedy
    left-to-right iteration (they would have been encountered as `nxt` first),
    so this check catches only other-motif tracts that lie between the two
    candidates and would make a merge biologically unjustified.
    """
    lo = bisect.bisect_left(all_starts_sorted, gap_start)
    return lo >= len(all_starts_sorted) or all_starts_sorted[lo] >= gap_end


def merge_tracts(tracts, k, get_gap_fn, canonical_motif, model, all_starts_sorted):
    """
    Greedy left-to-right merge of a start-sorted list of tracts under one impurity model.

    Parameters
    ----------
    tracts : list of dict
        Each dict must contain at least: start (int), end (int), copies (float).
    k : int
        Period length.
    get_gap_fn : callable(start, end) -> str
        Returns the reference sequence in the half-open interval [start, end).
    canonical_motif : str
        The canonical motif string for this group (length k).
    model : int
        0 or 1.
    all_starts_sorted : list of int
        Sorted list of tract start positions across ALL motifs on this contig.
        Used to ensure the two candidate tracts are genomically adjacent (no
        other repeat tract of any motif starts between them).

    Returns
    -------
    list of dict  — merged tracts, each with an added 'num_merged' key.
    Overlapping tracts (d < 0) are never merged and emitted as-is.
    """
    if not tracts:
        return []

    merged = []
    cur = dict(tracts[0])
    cur["num_merged"] = 1

    for nxt in tracts[1:]:
        d = nxt["start"] - cur["end"]

        can_merge = False
        if model == 0:
            # M0: exact baseline — never merge across any gap.
            pass

        elif model == "M_fission":  # future model — not active in workflow
            # 1-bp fission/fusion: merge when d == 0 (abutting) or d == 1 (one-base gap).
            if d == 0:
                # Abutting tracts: always merge — no gap to be occupied.
                can_merge = True
            elif d == 1:
                # Before evaluating the gap sequence, confirm no other repeat
                # tract starts in the gap.  If another motif sits between these
                # two tracts they are not genomically adjacent and must not merge.
                if not gap_is_clear(cur["end"], nxt["start"], all_starts_sorted):
                    pass  # blocked — a different repeat occupies the gap
                else:
                    gap_seq = get_gap_fn(cur["end"], nxt["start"])
                    mismatches = compute_mismatches(gap_seq, cur["end"] % k, canonical_motif, k)
                    if mismatches == 1:
                        # Condition 1: point mutation — gap base is a single substitution.
                        # e.g. ATATATATCTATATATATA: gap 'C' replaces expected 'A'.
                        can_merge = True
                    elif mismatches == 0:
                        # Condition 2: 1-bp indel — gap base is a correct motif continuation
                        # but the following tract is phase-shifted (insertion or deletion).
                        # e.g. ATATATATTATATAT: extra 'T' causes phase offset.
                        can_merge = True
                    # d == 1 guarantees len(gap_seq) == 1, so mismatches ∈ {0, 1}; no other case.

        if can_merge:
            cur["end"] = nxt["end"]
            cur["copies"] = float(cur["copies"]) + float(nxt["copies"])
            cur["num_merged"] += 1
        else:
            merged.append(cur)
            cur = dict(nxt)
            cur["num_merged"] = 1

    merged.append(cur)
    return merged


def read_trs_raw(path):
    """
    Read a *_k{k}_trs_raw.tsv file.
    Returns (tracts_by_motif, species, sample):
        tracts_by_motif : dict  canonical_motif -> list of tract dicts
        species, sample : str
    """
    tracts_by_motif = defaultdict(list)
    species = sample = None

    with open(path) as f:
        header = f.readline().strip().split("\t")
        idx = {col: i for i, col in enumerate(header)}

        for line in f:
            fields = line.strip().split("\t")
            if not fields or len(fields) < len(header):
                continue
            tract = {col: fields[idx[col]] for col in header}
            tract["start"] = int(tract["start"])
            tract["end"] = int(tract["end"])
            tract["copies"] = float(tract["copies"])
            species = tract["species"]
            sample = tract["sample"]
            tracts_by_motif[tract["canonical_motif"]].append(tract)

    return tracts_by_motif, species, sample


def main():
    trs_raw    = sys.argv[1]
    contig     = sys.argv[2]
    k          = int(sys.argv[3])
    assembly   = sys.argv[4]
    out_prefix = sys.argv[5]

    tracts_by_motif, species, sample = read_trs_raw(trs_raw)
    fasta = pysam.FastaFile(assembly)

    out_files = {}
    for m in models:
        path = f"{out_prefix}_M{m}_merged.tsv"
        out_files[m] = open(path, "w")
        out_files[m].write("\t".join(header) + "\n")

    # Build a globally sorted list of all tract start positions across all motifs.
    # Used by merge_tracts to verify genomic adjacency before any merge.
    all_starts_sorted = sorted(
        t["start"] for motif_tracts in tracts_by_motif.values() for t in motif_tracts
    )

    def get_gap_fn(start, end, _contig=contig, _fasta=fasta):
        if end <= start:
            return ""
        return _fasta.fetch(_contig, start, end).upper()

    # Collect all merged tracts per model, then sort by start position before writing.
    results = {m: [] for m in models}
    for canonical_motif, tracts in tracts_by_motif.items():
        tracts.sort(key=lambda t: t["start"])
        for m in models:
            for t in merge_tracts(tracts, k, get_gap_fn, canonical_motif, m, all_starts_sorted):
                results[m].append((t["contig"], t["start"], t["end"],
                                   canonical_motif, k, t["copies"],
                                   t["num_merged"], m, species, sample))

    for m in models:
        for row in sorted(results[m], key=lambda r: r[1]):  # sort by start
            out_files[m].write("\t".join(str(x) for x in row) + "\n")

    fasta.close()
    for m in models:
        out_files[m].close()


if __name__ == "__main__":
    main()
