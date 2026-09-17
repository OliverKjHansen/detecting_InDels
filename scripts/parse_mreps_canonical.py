import re
import sys
from collections import Counter

from parse_trf import canonicalise, hamming_pct_matches

# Usage: parse_mreps_canonical.py <mreps_txt> <contig> <k> <species> <sample> <model>
#
# Parses mreps output for a single period k and emits the unified 11-column
# repeat BED shared with parse_trf.py:
#   contig, start, end, canonical_motif, k, copies,
#   trf_purity, model_purity, model, species, sample
#
# Works for exact runs (M0, no -res) and runs with -res (M2): mreps prints the
# same table for both. Each repeat is one line:
#       19  ->        26 :     8   <2>   [4.00]   0.167     AC AT AC AC
#   from -> to : size <period> [exponent] err-rate copies
#
# Purity:
#   trf_purity   — mreps' own measure, (1 - err-rate) * 100. err-rate is
#                  mismatches / (length - period), where a single odd base
#                  between two matching copies counts once. 100.00 when exact.
#   model_purity — the Hamming measure parse_trf.py uses: the repeat compared
#                  base by base with its consensus motif tiled from the start.
#
# Motif: mreps prints the repeat split into copies, and in an impure repeat any
# copy can carry the mismatch. The motif is therefore the consensus of the full
# copies: the most common base at each position, ties going to the first copy.
#
# Coordinates: mreps is 1-based inclusive (201 -> 220 is 20 bp). The BED is
# 0-based half-open, so start = from - 1 and end = to.
#
# Canonical motif: lexicographic minimum over the rotations of the motif and of
# its reverse complement, so a repeat and the same repeat read from the other
# strand share one motif (AC and GT both become AC).

REPEAT_LINE = re.compile(
    r"^\s*(\d+)\s*->\s*(\d+)\s*:"   # from -> to :
    r"\s*(\d+)"                     # size
    r"\s*<(\d+)>"                   # <period>
    r"\s*\[([\d.]+)\]"              # [exponent]
    r"\s*([\d.]+)"                  # err-rate
    r"\s*(.*)$"                     # copies
)


def consensus(copies, k):
    """Most common base at each of the k positions across the full copies."""
    full = [copy.upper() for copy in copies if len(copy) == k]
    motif = ""
    for j in range(k):
        counts = Counter(copy[j] for copy in full)
        motif += max(counts, key=lambda base: (counts[base], base == full[0][j]))
    return motif


def main():
    input_file = sys.argv[1]
    contig     = sys.argv[2]
    k          = int(sys.argv[3])
    species    = sys.argv[4]
    sample     = sys.argv[5]
    model      = sys.argv[6]

    print("contig", "start", "end", "canonical_motif", "k", "copies",
          "trf_purity", "model_purity", "model", "species", "sample", sep="\t")

    with open(input_file) as handle:
        for line in handle:
            match = REPEAT_LINE.match(line)
            if not match:
                continue

            from_pos, to_pos, _size, period, exponent, err_rate, sequence = match.groups()
            copies = sequence.split()

            # -minperiod/-maxperiod should guarantee the period; a repeat
            # without its sequence (-noprint) cannot be given a motif.
            if int(period) != k or not any(len(copy) == k for copy in copies):
                continue

            start = int(from_pos) - 1
            end   = int(to_pos)

            motif = consensus(copies, k)
            mreps_purity = (1 - float(err_rate)) * 100
            model_purity = hamming_pct_matches("".join(copies), motif)

            print(contig, start, end, canonicalise(motif), k, exponent,
                  f"{mreps_purity:.2f}", f"{model_purity:.2f}", model,
                  species, sample, sep="\t")


if __name__ == "__main__":
    main()
