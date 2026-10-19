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
# germline_prefix) minus the highest among somatic_types, cut into the same bins
# by the same program, classify_reads.awk.
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
# big to load. Reads are scored and binned by classify_reads.awk, as in
# join_read_classification.sh; only the totals per bin are kept.
awk -F'\t' -v OFS='\t' -v prefix="$germline_prefix" -v somatic="$somatic_types" \
    -v edges="$bin_edges" -v labels="$bin_labels" \
    -f "$(dirname "$0")/classify_reads.awk" "$likelihoods" \
| awk -F'\t' -v OFS='\t' '
    NR == 1 { next }
    { reads[$5]++; if ($6 != "NA") cpgs[$5] += $6 }
    END { for (b in reads) print b, reads[b], cpgs[b] + 0 }
' | sort -k1,1n > "$out.body"

{
    printf 'classification_bin\treads\tCpGs\n'
    cat "$out.body"
} > "$out.tmp"
rm -f "$out.body"
mv "$out.tmp" "$out"

echo "$(basename "$out"): $(($(wc -l < "$out") - 1)) bins over" \
     "$(awk -F'\t' 'NR > 1 { n += $2 } END { print n }' "$out") reads" >&2
