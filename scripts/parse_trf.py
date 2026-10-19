import sys

# Usage: parse_trf.py <trf_ngs_txt> <k> <species> <sample> <model>
#                     [<contig> <offset> [<core_start> <core_end>]]
#
# Parses TRF -ngs output for a single period k and emits the unified repeat
# annotation TSV, tagged with the TRF model it came from (M1, M2).
#
# k can also be a range of periods, MIN-MAX, as for the tracts with periods
# above the repeat models' (the excluded group). Then a tract whose consensus
# motif is whole copies of a motif shorter than MIN is skipped: TRF can report
# a short-period repeat again at a multiple of its period, (AC)n as ACACACAC,
# and that tract belongs to the short period. For a single k, as M1 is split,
# nothing is skipped; overlap resolution keeps the smallest period instead.
#
# The optional arguments are for runs on a slice of a contig, where TRF
# sees a window rather than the whole sequence:
#   contig      name to report, since the FASTA header is "contig:start-end"
#   offset      0-based start of the slice, added to every coordinate
#   core_start  half-open range of contig positions this window owns; a tract
#   core_end    is kept only if its start falls inside it, so windows that
#               overlap do not report the same tract twice. Without them every
#               tract is kept.
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
# Output TSV columns (unified format shared with M0 and M2 from
# parse_mreps_canonical.py):
#   contig, start, end, canonical_motif, k, copies,
#   trf_purity, model_purity, model, species, sample


def canonicalise(motif):
    """
    Return the lexicographically smallest rotation of the motif or of its
    reverse complement.

    A repeat read from the other strand is the same locus, so AG and CT (and
    AC and GT) are one motif. Periods start at 2, so the older rotation-only
    form, which kept the A and T homopolymers apart, no longer applies.
    """
    motif = motif.upper()
    complement = motif.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return min(strand[i:] + strand[:i]
               for strand in (motif, complement)
               for i in range(len(strand)))


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


def has_shorter_period(motif, shorter_than):
    """True if the motif is whole copies of a motif shorter than `shorter_than` bp."""
    return any(len(motif) % d == 0 and motif == motif[:d] * (len(motif) // d)
               for d in range(1, min(shorter_than, len(motif))))


def main():
    input_file = sys.argv[1]
    # One period k, or a range MIN-MAX.
    low, dash, high = sys.argv[2].partition("-")
    is_range = bool(dash)
    low, high  = int(low), int(high or low)
    species    = sys.argv[3]
    sample     = sys.argv[4]
    model      = sys.argv[5]

    # Window arguments; without them the whole contig was scanned.
    contig_name = sys.argv[6] if len(sys.argv) > 6 else None
    offset      = int(sys.argv[7]) if len(sys.argv) > 7 else 0
    core_start  = int(sys.argv[8]) if len(sys.argv) > 8 else None
    core_end    = int(sys.argv[9]) if len(sys.argv) > 9 else None

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
            if not low <= period <= high:
                # For a single k the awk pre-filter in the GWF target should
                # guarantee this; kept as a safety check.
                continue

            start_1based = int(fields[0])
            end_1based   = int(fields[1])
            copies       = float(fields[3])
            trf_purity   = float(fields[5])  # TRF's SW-alignment % matches
            motif_seq    = fields[13]         # TRF consensus motif (used for tiling)
            repeat_seq   = fields[14]         # actual genomic repeat sequence

            if is_range and has_shorter_period(motif_seq.upper(), low):
                continue

            # 0-based half-open coordinates, shifted to contig coordinates
            start = start_1based - 1 + offset
            end   = end_1based + offset

            # Tracts starting outside this window's core belong to a neighbour.
            if core_start is not None and not (core_start <= start < core_end):
                continue

            canonical     = canonicalise(motif_seq)
            model_purity  = hamming_pct_matches(repeat_seq, motif_seq)

            print(contig_name or current_contig, start, end, canonical, period, copies,
                  f"{trf_purity:.2f}", f"{model_purity:.2f}", model,
                  species, sample, sep="\t")


if __name__ == "__main__":
    main()
