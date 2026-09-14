import sys

# Usage: parse_trf.py <trf_ngs_txt> <k> <species> <sample> <model>
#
# Parses TRF -ngs output for a single period k and emits the unified repeat
# annotation TSV, tagged with the TRF model it came from (M1, M2).
#
# Two impurity measures are reported:
#
#   trf_purity    — TRF's own SW-alignment-based % matches (column 6, 1-indexed).
#                   Computed internally by TRF using Smith-Waterman; can absorb
#                   indels, so it reflects biological impurity more accurately for
#                   tracts with insertions/deletions relative to the consensus.
#
#   model_purity  — Our Hamming-based measure: for each position j in the
#                   repeat_sequence (TRF column 15), compare to motif[j % period]
#                   (TRF's own consensus motif, column 14, tiled from offset 0).
#                   model_purity = matches / tract_length * 100
#                   Strictly positional — no gap concept.  Consistent and
#                   reproducible without re-running the aligner.
#
# Coordinate convention: TRF 1-based inclusive → 0-based half-open
#   start = start_1based - 1,  end = end_1based (already exclusive)
#
# Output TSV columns (unified format shared with M0 from impurity_models.py):
#   contig, start, end, canonical_motif, k, copies,
#   trf_purity, model_purity, model, species, sample


def canonicalise(motif):
    """Return the lexicographically smallest rotation of the motif (rotations only, no RC)."""
    motif = motif.upper()
    k = len(motif)
    return min(motif[i:] + motif[:i] for i in range(k))


def hamming_pct_matches(repeat_seq, motif):
    """
    Compute % base matches between repeat_seq and motif tiled from offset 0.
    Returns model_purity as a float in [0.0, 100.0].

    Each position j is compared to motif[j % len(motif)].  The motif should be
    the TRF consensus (not the canonical rotation) to stay in phase with the
    reported repeat sequence.

    Examples
    --------
    Perfect AT repeat (TRF row 1: trf_purity=100, pct_indels=0):
        repeat_seq = "ATATATATATATATATATATATATATATATATATATATATATATATATATATAT"   (53 bp)
        motif      = "AT"
        tiled      = "ATATATATATATATATATATATATATATATATATATATATATATATATATATAT"
        matches    = 53 / 53 → model_purity = 100.00  (trf_purity = 100.0)

    Impure AT repeat (TRF row 2: trf_purity=83, pct_indels=9):
        repeat_seq = "ATATATATAATATATATATATATTTTATATATATACATATAT"   (41 bp)
        motif      = "AT"
        tiled      = "ATATATATATATATATATATATATATATATATATATATATAT"
        pos  8: A==A ✓    pos  9: A!=T ✗    pos 22: T!=A ✗    pos 23: T==T ✓
        pos 24: T!=A ✗    pos 25: T==T ✓    pos 26: T!=A ✗    pos 34: C!=A ✗  ...
        (TRF's trf_purity=83 uses SW alignment and absorbs 9% indels; our
        Hamming is purely positional and will give a different value.)
        ≈ 34 matches out of 41 → model_purity ≈ 82.93  (trf_purity = 83.0)

    Impure TG repeat (TRF row 3: trf_purity=81, pct_indels=0):
        repeat_seq = "TGTGTGTCTGTGTATGTGTATCTGAGTGTGTGTGTGTGTGTGTGTGTGTGT"   (51 bp)
        motif      = "TG"
        tiled      = "TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT"
        pos  7: C!=G ✗    pos 10: A!=T ✗    pos 12: A!=G ✗    pos 14: A!=T ✗
        pos 18: C!=T ✗    pos 22: A!=T ✗    pos 24: G!=T? ...
        ≈ 41 matches out of 51 → model_purity ≈ 80.39  (trf_purity = 81.0)
    """
    n = len(repeat_seq)
    if n == 0:
        return 100.0
    matches = sum(
        1 for j, b in enumerate(repeat_seq)
        if b.upper() == motif[j % len(motif)].upper()
    )
    return matches / n * 100.0


def main():
    input_file = sys.argv[1]
    k          = int(sys.argv[2])
    species    = sys.argv[3]
    sample     = sys.argv[4]
    model      = sys.argv[5]

    out_header = (
        "contig", "start", "end", "canonical_motif", "k", "copies",
        "trf_purity", "model_purity", "model", "species", "sample"
    )
    print("\t".join(out_header))

    current_contig = None
    with open(input_file) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith("@"):
                current_contig = line[1:].strip()
                continue
            if not line.strip():
                continue

            fields = line.split()
            if len(fields) < 15:
                continue

            period = int(fields[2])
            if period != k:
                # awk pre-filter in the GWF target should guarantee this;
                # kept as a safety check.
                continue

            start_1based = int(fields[0])
            end_1based   = int(fields[1])
            copies       = float(fields[3])
            trf_purity   = float(fields[5])  # TRF's SW-alignment % matches
            motif_seq    = fields[13]         # TRF consensus motif (used for tiling)
            repeat_seq   = fields[14]         # actual genomic repeat sequence

            # 0-based half-open coordinates
            start = start_1based - 1
            end   = end_1based

            canonical     = canonicalise(motif_seq)
            model_purity  = hamming_pct_matches(repeat_seq, motif_seq)

            print(current_contig, start, end, canonical, k, copies,
                  f"{trf_purity:.2f}", f"{model_purity:.2f}", model,
                  species, sample, sep="\t")


if __name__ == "__main__":
    main()
