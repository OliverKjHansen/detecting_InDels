# Classify reads as germline-like or somatic-like from a table of per-read
# cell-type log-likelihoods: one row per read, its ID first, one column per
# cell type (higher = better fit), and CpGs_onRead. Shared by
# read_classification.sh, join_read_classification.sh and
# classification_read_counts.sh, so they score and bin reads the same way.
#
#   awk -F'\t' -v OFS='\t' -f classify_reads.awk \
#       -v prefix=Sperm. -v somatic=Blood.Monocytes,... \
#       -v edges=-10,...,10 -v labels=-12,...,12 [-v wanted=reads.txt] likelihoods.tsv
#
# Writes a header and, per read:
#   readID, max_germline (highest likelihood among the columns starting with
#   prefix), max_somatic (highest among the somatic columns),
#   classification_diff (their difference; positive looks germline),
#   classification_bin (the label of the first edge the difference is below,
#   or the last label), CpGs_onRead (NA without the column).
# With wanted, a file of read IDs, only those reads.

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
    if (wanted != "")
        while ((getline id < wanted) > 0) keep[id]
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
    print "readID", "max_germline", "max_somatic", "classification_diff", "classification_bin", "CpGs_onRead"
    next
}

wanted != "" && !($1 in keep) { next }

{
    first = 1
    for (i in germline_col) { v = $i + 0; if (first || v > g) { g = v; first = 0 } }
    first = 1
    for (i in somatic_col)  { v = $i + 0; if (first || v > s) { s = v; first = 0 } }
    d = g - s
    print $1, g, s, d, bin(d), (cpg_col ? $cpg_col : "NA")
}
