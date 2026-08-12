# InDel Call Plan

Two parallel pipelines for calling InDels from self-aligned PacBio HiFi long reads.
Both use pysam pileup over a diploid assembly BAM, validated with haploid remapping.
The key distinction is *where* in the genome they scan.

---

## Shared inputs (per sample)

| Input | Path |
|---|---|
| Diploid assembly | `~/DerivedData/fasta/HiFi/{species}/{sample}/diploid_assembly/{sample}_diploid.fa` |
| Haploid assembly | `~/DerivedData/fasta/HiFi/{species}/{sample}/haploid_assembly/{sample}_haploid.fa` |
| Self-aligned BAM (diploid) | `~/DerivedData/bam/HiFi/{species}/{sample}/kinetics/{sample}_kinetics_diploid.bam` |
| Average coverage | `microsatellites/ratemicrosat/coverage/{sample}_avg_coverage.txt` |
| Tissue classification | `microsatellites/marmosetanalysis/{species}/{sample}_sorted_tissue_classification.tsv` |
| TR BED M0 (exact/pure, mreps) | `output/tandem_repeats/{species}/{sample}/{contig}_M0.bed` ← from `repetitiveDNA.py` |
| TR BED M1 (impurity-tolerant, TRF) | `output/tandem_repeats/{species}/{sample}/{contig}_M1.bed` ← from `repetitiveDNA.py` |

---

## Pipeline 1 — Non-repetitive InDel Call

**Workflow file**: `indel_mutation_call.py` ✅ (largely complete)

Scans positions **outside** the filter BED (microsatellites + unresolved haplotypes).
Calls 1–3 bp insertions and deletions.

### Steps

**Step 1 — Build filter BED** *(per contig)*
- Merge microsatellite intervals (`{contig}_trs_tidy.tsv`, cols 1-3) and unresolved haplotype regions (`{contig}_unresolved_haplotypes.bed`) using `bedtools merge`.
- Output: `output/small_indels/{species}/{sample}/{contig}_filters.tsv`
- Purpose: defines regions to *exclude* from scanning.

**Step 2 — Detect 1–3 bp InDels** *(per contig)*
- Script: `scripts/small_indels.py` → `utilities.small_indels()`
- Scans pileup; skips positions inside filter BED (IntervalTree exclusion).
- Filters: coverage 0.5×–2× average, homopolymer check, InDel size 1–3 bp, kmerQ ≥ 93.
- Output: `output/small_indels/{species}/{sample}/{contig}_small_indels.tsv`
- Columns: `#contig, start, end, ref, alt, indel_type, unique, indel_string, tissue, read_length, np, ec, sm, sx, sa, sm_pre, sx_pre, ToreadEnd, kmerQ, position, readID, diploid_direction`

**Step 3 — Extract mutation FASTQ** *(per contig, parallel with Step 4)*
- Script: `scripts/build_mutation_fastq.py`
- Extracts reads supporting detected InDels into FASTQ for haploid remapping.
- Output: `output/mutation_fastq/{species}/{sample}/{contig}_remapping.fastq`

**Step 4 — Count opportunities** *(per contig, parallel with Step 3)*
- Script: `scripts/opportunities.py` → `utilities.indel_opportunities()`
- Counts non-InDel positions per (context, tissue, mapQ, kmerQ), using the same filter BED for exclusion.
- Output: `output/opportunities/{species}/{sample}/{contig}_opportunities.tsv`
- Columns: `contig, context, tissue, mapQ, kmerQ, count`

**Step 5 — Haploid remapping** *(per sample, depends on all Step 3 contigs)*
- Concatenate all per-contig FASTQs → single FASTQ.
- Remap to haploid assembly: `pbmm2 align --preset HIFI --sort`
- Extract alignment info: `samtools view | awk` → sorted TSV.
- Output: `output/haploid_mapping/{species}_{sample}_haploid_locations.tsv`
- Columns: `readID, bitflag, contig_hap, pos_1based_hap, mapq, cigar_hap`

**Step 6 — Haploid concordance filter** *(per contig, depends on Steps 2 and 5)*
- Script: `scripts/filtering_haploid_mappping.py`
- Checks each InDel against the haploid mapping; retains only InDels confirmed in haploid with matching type and size.
- Appends `species` and `sample` columns.
- Output: `output/small_indels/{species}/{sample}/{contig}_small_indels_haploid_filtered.tsv`

**Step 7 — Merge and plot** *(per sample, depends on all Step 6 contigs and all Step 4 contigs)*
- `awk` merge of per-contig filtered TSVs and opportunities TSVs (header from first file only).
- R diagnostic plots: `scripts/plot_indels.R`
- Output: `output/final/{sample}_haploid_combined.tsv`, `output/final/{sample}_opportunities_combined.tsv`, `output/plots/{sample}_indel_plots.pdf`

**Step 8 — Final filtering** *(per sample, depends on Step 7)*
- Scripts: `scripts/filtering_small_indels.py` and `scripts/filtering_opportunities.py`
- Applies: kmerQ ≥ 93, haploid concordance check, min read-end distance 100 bp.
- Produces step-by-step filter report.
- Output: `output/final/{sample}_combined_filtered.tsv`, `output/final/{sample}_opportunities_combined_filtered.tsv`, `output/final/{sample}_filter_report.tsv`

### Open items (Pipeline 1)
- "Only one quality InDel per read" filter (flagged in README, not yet implemented).
- Distance-to-read-end filter (parameter exists in `filtering_small_indels.py`, not yet wired).

---

## Pipeline 2 — Tandem Repeat InDel Call

**Workflow file**: `repetitive_indel_call.py` 🆕

Scans positions **inside** M0 or M1 TR BED regions.
Calls InDels whose size is a multiple of the repeat period (any size).
Runs twice per sample: once for M0 (exact/pure repeats) and once for M1 (impurity-tolerant).

**Prerequisite**: M0 and M1 BED files must exist from `repetitiveDNA.py`.
Run: `gwf -f repetitiveDNA.py run` before starting this pipeline.

### Steps

**Step 1 — Detect period-multiple InDels** *(per contig per model)*
- Script: `scripts/tandem_repeat_indels.py` → `utilities.tandem_repeat_indels()` 🆕
- Loads M0 or M1 BED into IntervalTree (with metadata payload per interval).
- Scans pileup; only processes positions **inside** a TR interval (inverted from Pipeline 1).
- Filters:
  - Coverage 0.5×–2× average.
  - InDel size filter: `abs(indel) % k == 0` where k is the period of the overlapping TR interval. Homopolymer filter is NOT applied (k=1 repeats ARE homopolymers — they are the signal).
  - When a position overlaps multiple TR intervals, the matching interval with the largest k is used.
  - kmerQ is recorded but not used as a hard filter here (applied at final filtering step).
- Extra output columns vs. Pipeline 1: `canonical_motif, k, copies, trf_purity, model_purity, model`
- Output: `output/tandem_repeat_indels/{species}/{sample}/{contig}_{model}_tr_indels.tsv`
- Columns: `#contig, start, end, ref, alt, indel_type, unique, indel_string, tissue, read_length, np, ec, sm, sx, sa, sm_pre, sx_pre, ToreadEnd, kmerQ, position, readID, diploid_direction, canonical_motif, k, copies, trf_purity, model_purity, model`

**Step 2 — Extract mutation FASTQ** *(per contig per model, parallel with Step 3)*
- Reuse `scripts/build_mutation_fastq.py` unchanged (reads `position` and `readID` columns).
- Output: `output/tandem_repeat_mutation_fastq/{species}/{sample}/{contig}_{model}_remapping.fastq`

**Step 3 — Count TR opportunities** *(per contig per model, parallel with Step 2)*
- Script: `scripts/tandem_repeat_opportunities.py` → `utilities.tandem_repeat_opportunities()` 🆕
- Loads M0 or M1 BED into IntervalTree (same as Step 1).
- Scans pileup inside TR intervals; counts non-InDel reads per position.
- For each position, uses the overlapping TR interval with the smallest k (most frequent period).
- kmerQ window uses period k as the InDel length proxy.
- Output: `output/tandem_repeat_opportunities/{species}/{sample}/{contig}_{model}_tr_opportunities.tsv`
- Columns: `contig, canonical_motif, k, tissue, mapQ, kmerQ, count`

**Step 4 — Haploid remapping** *(per sample per model, depends on all Step 2 contigs)*
- Reuse same `pbmm2 align --preset HIFI --sort` command as Pipeline 1.
- Output: `output/tandem_repeat_haploid_mapping/{species}_{sample}_{model}_haploid_locations.tsv`

**Step 5 — Haploid concordance filter** *(per contig per model, depends on Steps 1 and 4)*
- Reuse `scripts/filtering_haploid_mappping.py` unchanged.
- Output: `output/tandem_repeat_indels/{species}/{sample}/{contig}_{model}_tr_indels_haploid_filtered.tsv`

**Step 6 — Merge and plot** *(per sample per model, depends on all Step 5 contigs and all Step 3 contigs)*
- `awk` merge of per-contig filtered TSVs and opportunities TSVs.
- R plots reusing `scripts/plot_indels.R` (period-stratified plots to be added separately).
- Output: `output/final/{sample}_{model}_tr_haploid_combined.tsv`, `output/final/{sample}_{model}_tr_opportunities_combined.tsv`, `output/plots/{sample}_{model}_tr_indel_plots.pdf`

**Step 7 — Final filtering** *(per sample per model, depends on Step 6)*
- Reuse `scripts/filtering_small_indels.py` and `scripts/filtering_opportunities.py`.
- TR-specific filters can be added via the existing column-based filter framework: e.g., minimum `model_purity`, minimum `copies`.
- Output: `output/final/{sample}_{model}_tr_combined_filtered.tsv`, `output/final/{sample}_{model}_tr_opportunities_combined_filtered.tsv`, `output/final/{sample}_{model}_tr_filter_report.tsv`

### Open items (Pipeline 2)
- Period-stratified diagnostic plots (currently reuses Pipeline 1 R script; a TR-specific R script may be needed).
- Decision on minimum `model_purity` threshold for final filtering.
- Decide whether k=1 (homopolymers) InDels should be separated from k>1 (microsatellites) in downstream analysis.

---

## New files to create

| File | Purpose |
|---|---|
| `repetitive_indel_call.py` | GWF workflow for Pipeline 2 |
| `scripts/tandem_repeat_indels.py` | CLI wrapper for `utilities.tandem_repeat_indels()` |
| `scripts/tandem_repeat_opportunities.py` | CLI wrapper for `utilities.tandem_repeat_opportunities()` |
| `utilities.py` (additions) | `load_tr_bed_with_tree()`, `tandem_repeat_indels()`, `tandem_repeat_opportunities()` |

## Scripts reused unchanged (Pipeline 2)

| Script | Reuse in step |
|---|---|
| `scripts/build_mutation_fastq.py` | Step 2 |
| `scripts/filtering_haploid_mappping.py` | Step 5 |
| `scripts/filtering_small_indels.py` | Step 7 |
| `scripts/filtering_opportunities.py` | Step 7 |
| `scripts/plot_indels.R` | Step 6 |

---

## Running order

```
# 1. Build M0/M1 BED files (prerequisite for Pipeline 2)
gwf -f repetitiveDNA.py run

# 2. Run non-repetitive InDel pipeline
gwf -f Indel_mutation_call.py run

# 3. Run tandem repeat InDel pipeline (can run in parallel with step 2)
gwf -f repetitive_indel_call.py run
```
