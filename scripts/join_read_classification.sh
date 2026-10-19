#!/usr/bin/env bash
#
# Join per-read cell-type likelihoods onto a table of InDel calls.
#
# Usage:
#   join_read_classification.sh <calls.tsv> <likelihoods.tsv> <out.tsv>
#                               <germline_prefix> <somatic_types> <bin_edges> <bin_labels>
#
# The last four come from the `classification` block of config.yaml; the lists
# are comma-separated.
#
# The likelihood table has one row per read: readsID, chr, start, one column per
# cell type holding a log-likelihood (higher = better fit), then CpGs_onRead and
# warning. It covers every read of the sample and is far too big to load whole,
# so it is streamed once and only the reads present in <calls.tsv> are kept.
#
# For each read:
#   max_germline        highest likelihood among the germline columns, i.e. the
#                       ones whose name starts with germline_prefix (Sperm.)
#   max_somatic         highest likelihood among somatic_types
#   classification_diff max_germline - max_somatic; positive means the read
#                       looks germline
#   classification_bin  the difference cut into bins: the label of the first
#                       bin edge it is below, or the last label
#   CpGs_onRead         carried through: the evidence behind the likelihoods.
#                       A read with few CpGs scores near 0 for lack of
#                       information, not because it is genuinely ambiguous.
#
# Every row of <calls.tsv> is kept. Calls whose read is absent from the
# likelihood table get NA in the five new columns.

set -euo pipefail

if [ "$#" -ne 7 ]; then
    echo "Usage: $(basename "$0") <calls.tsv> <likelihoods.tsv> <out.tsv> <germline_prefix> <somatic_types> <bin_edges> <bin_labels>" >&2
    exit 1
fi

calls=$1
likelihoods=$2
out=$3
germline_prefix=$4
somatic_types=$5
bin_edges=$6
bin_labels=$7

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# 1. The reads this call table needs. The header may or may not start with '#'.
awk -F'\t' '
    FNR == 1 {
        for (i = 1; i <= NF; i++) { name = $i; sub(/^#/, "", name); if (name == "readID") col = i }
        if (!col) { print "no readID column in " FILENAME > "/dev/stderr"; exit 1 }
        next
    }
    { print $col }
' "$calls" | sort -u > "$work/wanted_reads.txt"

# 2. Scores per read, reading the big table once and keeping only those reads.
awk -F'\t' -v OFS='\t' -v prefix="$germline_prefix" -v somatic="$somatic_types" \
    -v edges="$bin_edges" -v labels="$bin_labels" '
    # The label of the first edge d is below, or the last label.
    function bin(d,    i) {
        for (i = 1; i <= n_edges; i++) if (d < edge[i]) return label[i]
        return label[n_edges + 1]
    }

    BEGIN {
        n_edges = split(edges, edge, ",")
        for (i = 1; i <= n_edges; i++) edge[i] += 0
        if (split(labels, label, ",") != n_edges + 1) {
            print "need one more bin label than bin edges" > "/dev/stderr"; exit 1
        }
    }

    NR == FNR { wanted[$1]; next }

    FNR == 1 {
        split(somatic, list, ",")
        for (i in list) is_somatic[list[i]]
        for (i = 1; i <= NF; i++) {
            if (index($i, prefix) == 1) germline_col[i]
            if ($i in is_somatic)       somatic_col[i]
            if ($i == "CpGs_onRead")    cpg_col = i
        }
        if (length(germline_col) == 0) { print "no columns start with " prefix > "/dev/stderr"; exit 1 }
        if (length(somatic_col) == 0)  { print "none of " somatic " found"    > "/dev/stderr"; exit 1 }
        print "readID", "max_germline", "max_somatic", "classification_diff", "classification_bin", "CpGs_onRead"
        next
    }

    !($1 in wanted) { next }

    {
        first = 1
        for (i in germline_col) { v = $i + 0; if (first || v > g) { g = v; first = 0 } }
        first = 1
        for (i in somatic_col)  { v = $i + 0; if (first || v > s) { s = v; first = 0 } }
        d = g - s
        print $1, g, s, d, bin(d), (cpg_col ? $cpg_col : "NA")
    }
' "$work/wanted_reads.txt" "$likelihoods" > "$work/read_scores.tsv"

# 3. Left join onto the calls: every call row is kept, unmatched reads get NA.
awk -F'\t' -v OFS='\t' '
    NR == FNR {
        if (FNR > 1) { germline[$1] = $2; somatic[$1] = $3; diff[$1] = $4; bin[$1] = $5; cpg[$1] = $6 }
        next
    }
    FNR == 1 {
        for (i = 1; i <= NF; i++) { name = $i; sub(/^#/, "", name); if (name == "readID") col = i }
        print $0, "max_germline", "max_somatic", "classification_diff", "classification_bin", "CpGs_onRead"
        next
    }
    {
        read = $col
        if (read in diff) print $0, germline[read], somatic[read], diff[read], bin[read], cpg[read]
        else              print $0, "NA", "NA", "NA", "NA", "NA"
    }
' "$work/read_scores.tsv" "$calls" > "$out.tmp"

mv "$out.tmp" "$out"

matched=$(($(wc -l < "$work/read_scores.tsv") - 1))
wanted=$(wc -l < "$work/wanted_reads.txt")
echo "$(basename "$calls"): $matched of $wanted reads classified" >&2
