import sys

# Usage: parse_mreps_canonical.py <mreps_txt> <contig> <k> <species> <sample>
#
# Parses mreps -fasta output for a single period k.
# Emits TSV: contig, start, end, canonical_motif, k, copies, phase, raw_motif, species, sample
#
# Coordinate convention (matching kolpakov.py):
#   mreps output is 1-based inclusive → shift start to 0-based, leave end as-is (exclusive).
#
# Canonical motif: lexicographic minimum over all rotations of the motif itself.
# Reverse complement is NOT included: AAAAA and TTTTT are distinct genomic tracts and
# must not be collapsed into the same motif group at this stage. Collapsing across strands
# would incorrectly merge tracts and break the adjacency-graph / impurity-model logic.
# Phase: start % k (0-based start), representing where in the motif cycle the tract begins.

def canonicalise(motif):
    """Return the lexicographically smallest rotation of the motif (rotations only, no RC)."""
    motif = motif.upper()
    k = len(motif)
    return min(motif[i:] + motif[:i] for i in range(k))


input_file = sys.argv[1]
contig      = sys.argv[2]
k           = int(sys.argv[3])
species     = sys.argv[4]
sample      = sys.argv[5]

with open(input_file) as f:
    lines = f.readlines()

# mreps -fasta output: skip the 15-line header and last 3 footer lines
# (same convention as kolpakov.py)
body = lines[15:-3]

print("contig", "start", "end", "canonical_motif", "k", "copies", "phase",
      "raw_motif", "species", "sample", sep="\t")

for line in body:
    line = line.strip()
    if not line:
        continue

    fields = line.split("\t")
    if len(fields) < 4:
        continue

    try:
        start_raw, end_raw = fields[0].replace(" ", "").replace(":", "").split("->")
        start_1based = int(start_raw)
        end_1based   = int(end_raw)

        copies   = fields[3].strip(" [] ")
        raw_motif = fields[-1].split()[0]  # first token after the sequence body
    except (ValueError, IndexError):
        continue

    # Safety net for malformed lines: -exp 2 guarantees ≥ 2 copies (exponent = length/period).
    try:
        if float(copies) < 2.0:
            continue
    except ValueError:
        continue

    # 0-based half-open coordinates #check this!!!!!!
    start = start_1based - 1
    end   = end_1based        # already exclusive

    canonical = canonicalise(raw_motif)
    phase      = start % k

    print(contig, start, end, canonical, k, copies, phase,
          raw_motif, species, sample, sep="\t")
