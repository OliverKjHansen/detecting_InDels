#!/usr/bin/env bash
#
# Count the sample's reads per germline/somatic bin.
#
# Usage:
#   classification_read_counts.sh <likelihoods.tsv> <out.tsv>
#                                 <germline_prefix> <somatic_types> <bin_edges> <bin_labels>
#
# The last four come from the `classification` block of config.yaml; the lists
# are comma-separated.
#
# The calls carry a classification bin, from join_read_classification.sh, but
# the opportunity tables do not: they count read positions per contig, context,
# tissue and mapQ, with read identity aggregated away. So the denominator a rate
# per bin needs,
#
#     P_b = the read positions contributed by the reads in bin b
#
# cannot be recovered from them. This script counts, over every read of the
# sample, how many reads fall in each bin and how many CpGs those reads carry.
# The bins can then be given a share of the sample's read positions P:
#
#     P_b ~= P * reads_b / reads_total      flat
#     P_b ~= P * CpGs_b  / CpGs_total       length-weighted
#
# The flat split assumes the positions per read do not depend on the bin. They
# do a little: a longer read carries more CpGs, is classified more confidently
# and lands in a more extreme bin, so the extreme bins hold longer reads than
# average and the flat split understates their share. CpGs per read stand in for
# read length, which is why both columns are written.
#
# The scores per read are those of join_read_classification.sh: the highest
# likelihood among the germline columns (the ones whose name starts with
# germline_prefix) minus the highest among somatic_types, cut into the same bins.
# Both scripts take the bins from config.yaml; tests/test_classification_bins.py
# checks that they bin alike.
#
# Output: one row per bin, in bin order.
#
#     classification_bin   reads   CpGs

set -euo pipefail

if [ "$#" -ne 6 ]; then
    echo "Usage: $(basename "$0") <likelihoods.tsv> <out.tsv> <germline_prefix> <somatic_types> <bin_edges> <bin_labels>" >&2
    exit 1
fi

likelihoods=$1
out=$2
germline_prefix=$3
somatic_types=$4
bin_edges=$5
bin_labels=$6

if [ ! -r "$likelihoods" ]; then
    echo "cannot read $likelihoods" >&2
    exit 1
fi

# One row per read, streamed: the table covers the whole sample and is far too
# big to load. Nothing is kept per read, only the totals per bin.
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
        next
    }

    {
        first = 1
        for (i in germline_col) { v = $i + 0; if (first || v > g) { g = v; first = 0 } }
        first = 1
        for (i in somatic_col)  { v = $i + 0; if (first || v > s) { s = v; first = 0 } }
        b = bin(g - s)
        reads[b]++
        if (cpg_col) cpgs[b] += $cpg_col + 0
    }

    END {
        for (b in reads) print b, reads[b], cpgs[b] + 0
    }
' "$likelihoods" | sort -k1,1n > "$out.body"

{
    printf 'classification_bin\treads\tCpGs\n'
    cat "$out.body"
} > "$out.tmp"
rm -f "$out.body"
mv "$out.tmp" "$out"

echo "$(basename "$out"): $(($(wc -l < "$out") - 1)) bins over" \
     "$(awk -F'\t' 'NR > 1 { n += $2 } END { print n }' "$out") reads" >&2
