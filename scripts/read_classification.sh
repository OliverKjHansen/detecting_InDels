#!/usr/bin/env bash
#
# Classify every read of a sample as germline-like or somatic-like.
#
# Usage:
#   read_classification.sh <likelihoods.tsv> <out.tsv>
#                          <germline_prefix> <somatic_types> <bin_edges> <bin_labels>
#
# The last four come from the `classification` block of config.yaml; the lists
# are comma-separated. Writes one row per read (readID, max_germline,
# max_somatic, classification_diff, classification_bin, CpGs_onRead), scored
# and binned by classify_reads.awk. The scans read each read's bin from it, so
# calls and opportunities carry the same bins.

set -euo pipefail

if [ "$#" -ne 6 ]; then
    echo "Usage: $(basename "$0") <likelihoods.tsv> <out.tsv> <germline_prefix> <somatic_types> <bin_edges> <bin_labels>" >&2
    exit 1
fi

likelihoods=$1
out=$2

if [ ! -r "$likelihoods" ]; then
    echo "cannot read $likelihoods" >&2
    exit 1
fi

awk -F'\t' -v OFS='\t' -v prefix="$3" -v somatic="$4" -v edges="$5" -v labels="$6" \
    -f "$(dirname "$0")/classify_reads.awk" "$likelihoods" > "$out.tmp"
mv "$out.tmp" "$out"
echo "$(basename "$out"): $(($(wc -l < "$out") - 1)) reads classified" >&2
