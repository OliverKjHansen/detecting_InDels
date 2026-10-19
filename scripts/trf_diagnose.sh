#!/bin/bash
#SBATCH --account=mutationalscanning
#SBATCH --cpus-per-task=1
#SBATCH --mem=16g
#SBATCH --time=12:00:00
#SBATCH --job-name=trf_diagnose
#SBATCH --output=trf_diagnose_%j.log

# What makes TRF slow in one region of a contig.
#
# Splits the region into chunks and reports, per chunk:
#   compress_ratio  how well it compresses: ordinary sequence ~0.28-0.30,
#                   satellite arrays ~0.1
#   N_share         share of N bases
#   top_6mers       the three most common 6-mers and their shares; in a
#                   satellite one or two of them dominate
#   trf_seconds     TRF runtime with the parameters of workflows/repetitive_dna.py,
#                   or TIMEOUT when it runs over the time limit
#   tracts_k2to6    tracts of period 2-6 TRF reports, and the longest (bp)
#
# Usage: sbatch scripts/trf_diagnose.sh [assembly.fa] [contig] [start] [end] [chunk_bp] [time_limit]
# The defaults cover the window of M1_trf_h1tg000011l_0_J12013_diploid with
# trf_window 4 Mb and trf_overlap 1 Mb. start is 0-based, as in the workflow.
# The table is written to trf_diagnose_{contig}_{start}-{end}.tsv and the log.

set -euo pipefail

assembly=${1:-$HOME/mutationalscanning/DerivedData/fasta/HiFi/long_tailed_macaque/J12013/diploid_assembly/J12013_diploid.fa}
contig=${2:-h1tg000011l}
start=${3:-0}
end=${4:-5000000}
chunk=${5:-500000}
limit=${6:-30m}
trf=${TRF:-$HOME/mutationalscanning/Workspaces/oliver/microsatellites/development/TRF/build/src/trf}
trf_params="2 7 7 80 10 50 6"

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
report=trf_diagnose_${contig}_${start}-${end}.tsv

printf "start\tend\tcompress_ratio\tN_share\ttop_6mers\ttrf_seconds\ttracts_k2to6\tlongest_bp\n" | tee "$report"

for ((s = start; s < end; s += chunk)); do
    e=$(( s + chunk < end ? s + chunk : end ))
    samtools faidx "$assembly" "$contig:$((s + 1))-$e" > "$work/chunk.fa"

    stats=$(python3 - "$work/chunk.fa" <<'EOF'
import sys, zlib
from collections import Counter
seq = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">")).upper()
ratio = len(zlib.compress(seq.encode(), 9)) / len(seq)
kmers = Counter(seq[i:i + 6] for i in range(len(seq) - 5))
kmers.pop("NNNNNN", None)
total = sum(kmers.values()) or 1
top = ",".join(f"{k}:{n / total:.3f}" for k, n in kmers.most_common(3))
print(f"{ratio:.3f}\t{seq.count('N') / len(seq):.3f}\t{top}")
EOF
)

    t0=$SECONDS
    status=0
    timeout "$limit" "$trf" "$work/chunk.fa" $trf_params -ngs -h > "$work/trf.txt" 2> /dev/null || status=$?
    if [ "$status" -eq 124 ]; then
        secs="TIMEOUT(>$limit)"; tracts=NA; longest=NA
    elif [ "$status" -ne 0 ]; then
        secs="FAILED(exit $status)"; tracts=NA; longest=NA
    else
        secs=$(( SECONDS - t0 ))
        read -r tracts longest < <(awk '!/^@/ && NF >= 15 && $3 >= 2 {n++; l = $2 - $1 + 1; if (l > m) m = l}
                                       END {print n + 0, m + 0}' "$work/trf.txt")
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$s" "$e" "$stats" "$secs" "$tracts" "$longest" | tee -a "$report"
done
