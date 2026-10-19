# Project Overview

## What this project does

`detecting_InDels` is a bioinformatics analysis pipeline that detects **small insertions and deletions (InDels) in single PacBio HiFi long reads**, in order to estimate somatic/germline InDel mutation rates across primate species.

The central idea is **self-alignment**. Reads are mapped back to a de novo assembly built from those same reads. A read that disagrees with its own assembly at a given position carries a candidate mutation rather than a genuine reference difference. Because a single sequencing error looks identical to a real mutation at this level, most of the code is dedicated to *filtering out artefacts* rather than to calling variants.

Species and samples currently configured (`indel_mutation_call.py:9-17`, `repetitive_indel_call.py:4-8`, `repetitiveDNA.py:4-7`):

| Species | Samples |
|---|---|
| chimp | yoran |
| macaque | jako |
| marmoset | marmoset, lambiek |

Human, gorilla, baboon and gibbon entries exist but are commented out.

## Execution model

Everything is orchestrated by [GWF](https://gwf.app/), a Python workflow manager that submits jobs to a Slurm cluster. A "workflow file" is a plain Python script that builds a dependency graph of `gwf.target(...)` declarations; each target names its `inputs`, its `outputs`, and a bash recipe. GWF re-runs a target only when its inputs are newer than its outputs.

Default resources are declared once per workflow (`account='mutationalscanning'`), with heavier per-target overrides for the merge/remap stages.

Two structural rules shape the whole graph:

- **Per-contig fan-out, per-sample fan-in.** Detection runs independently per contig so it parallelises across the cluster; results are concatenated per sample with `awk 'FNR==1 && NR==1 || FNR>1'` (keep header from the first file only).
- **Only contigs longer than 1 Mb are processed** (`if int(length) > 10**6`), read from the assembly's `.fa.fai` index. Short contigs are skipped entirely.

There is no database, no service and no network I/O. Every stage communicates by writing TSV/BED/FASTQ files to `output/`.

## Repository layout

```
├── repetitiveDNA.py            # Workflow A — annotate tandem repeats (M0/M1 BED)
├── indel_mutation_call.py      # Workflow B — InDels OUTSIDE repeats
├── repetitive_indel_call.py    # Workflow C — InDels INSIDE repeats
├── scripts/
│   ├── utilities.py            # 1190 lines — all pileup-scanning analysis logic
│   ├── small_indels.py         # thin argv wrappers around utilities functions
│   ├── opportunities.py
│   ├── tandem_repeat_indels.py
│   ├── tandem_repeat_opportunities.py
│   ├── large_indels.py
│   ├── small_indels_summary.py
│   ├── build_mutation_fastq.py         # extract mutation-supporting reads → FASTQ
│   ├── filtering_haploid_mappping.py   # haploid concordance check
│   ├── filtering_small_indels.py       # final column-based filtering + report
│   ├── filtering_opportunities.py      # same, for the denominator
│   ├── parse_mreps_canonical.py        # mreps output → unified TSV
│   ├── parse_trf.py                    # TRF output → M1 BED
│   ├── impurity_models.py              # M0 annotation
│   └── plot_indels.R                   # ggplot2 diagnostics
├── InDel_call_plan.md          # design doc for workflows B and C
├── tandem_repeat_plan.md       # design doc for workflow A
└── todo.md
```

The architecture is deliberately two-layered: **`utilities.py` holds the algorithms**, and each `scripts/*.py` is a ~25-line wrapper that reads positional `sys.argv` arguments and calls one function. Wrappers write results to stdout; the GWF recipe redirects stdout to the target's output file. Diagnostics go to stderr.

> Note: the workflows do not execute the scripts from this repository. All three set
> `scripts = expanduser("~/mutationalscanning/Workspaces/oliver/endjoining/scripts")`,
> so the code that actually runs on the cluster is a separate checkout at that path.

## Shared inputs

| Input | Path pattern |
|---|---|
| Diploid assembly | `~/mutationalscanning/DerivedData/fasta/HiFi/{species}/{sample}/diploid_assembly/{sample}_diploid.fa` |
| Assembly index | same path + `.fai` — drives the contig loop |
| Haploid assembly | `.../haploid_assembly/{sample}_haploid.fa` |
| Self-aligned BAM | `~/mutationalscanning/DerivedData/bam/HiFi/{species}/{sample}/kinetics/{sample}_kinetics_diploid.bam` |
| Average coverage | `.../ratemicrosat/coverage/{sample}_avg_coverage.txt` — a single float |
| Tissue classification | `.../marmosetanalysis/{species}/{sample}_sorted_tissue_classification.tsv` — readID → tissue |

The BAM is a *kinetics* BAM, so it carries PacBio consensus tags (`np`, `ec`, `sm`, `sx`, `sa`) that are recorded alongside every call for later quality modelling.

---

## Workflow A — `repetitiveDNA.py`: tandem repeat annotation

Produces the repeat annotation that the other two workflows depend on. Two independent repeat models are built per contig, for periods **k = 1…6**:

- **M0 — `mreps`, exact/pure repeats.** `mreps -minperiod $k -maxperiod $k -exp 2 -allowsmall`. Output parsed by `parse_mreps_canonical.py`, then annotated by `impurity_models.py`. Purity is 100% by definition.
- **M1 — `TRF`, impurity-tolerant repeats.** `trf {contig}.fa 2 7 7 80 10 50 $k -ngs -h`, post-filtered with awk to period-`k` rows, parsed by `parse_trf.py`.

`parse_trf.py` records two distinct purity measures (`parse_trf.py:8-20`): `trf_purity` (TRF's own Smith-Waterman % matches, which absorbs indels) and `model_purity` (a positional Hamming comparison against the consensus motif tiled from offset 0). Motifs are canonicalised to their lexicographically smallest **rotation** — deliberately *without* reverse complement, so `A` and `T` homopolymers stay distinct.

The k-loop runs inside bash in a `mktemp -d` scratch dir with `trap 'rm -rf "$WORK"' EXIT`, so only two files per contig survive:

```
output/tandem_repeats/{species}/{sample}/{contig}_M0.bed
output/tandem_repeats/{species}/{sample}/{contig}_M1.bed
```

Both are then concatenated per sample. BED columns are unified across models:
`contig, start, end, canonical_motif, k, copies, trf_purity, model_purity, model, species, sample` (0-based half-open).

---

## Workflow B — `indel_mutation_call.py`: non-repetitive InDels

Scans positions **outside** repeats and calls 1–3 bp InDels.

**Step 1 — build the exclusion BED.** Merge microsatellite intervals and unresolved-haplotype regions with `bedtools merge` → `{contig}_filters.tsv`.

**Step 2 — detect InDels** (`scripts/small_indels.py` → `utilities.small_indels()`). Iterates `bam.pileup()` with `stepper="nofilter"`, `min_base_quality=0`, `min_mapping_quality=0` — pysam's own filtering is disabled so that all filtering is explicit in the code. Per reference position:

- skip unless `avg_coverage/2 ≤ nsegments ≤ avg_coverage*2`;
- skip if the position falls in the exclusion BED (`IntervalTree` lookup);
- mark the position `unique` or `non-unique` depending on how many reads carry an indel there;
- per read, keep only `1 ≤ |indel| ≤ 3`;
- reject homopolymer contexts;
- record PacBio tags, `ToreadEnd`, and a 4-mer minimum base quality (`kmerQ`) spanning the InDel and its flanks.

**Step 3 — extract supporting reads** (`build_mutation_fastq.py`) into per-contig FASTQ.

**Step 4 — count opportunities** (`utilities.indel_opportunities()`), the denominator for a rate: non-InDel read positions tallied by `(context, tissue, mapQ, kmerQ)`.

**Step 5 — haploid remapping.** Concatenate FASTQs and remap to the *haploid* assembly with `pbmm2 align --preset HIFI --sort`, then flatten to `readID, bitflag, contig, pos, mapq, cigar` with `samtools view | awk | sort`.

**Step 6 — haploid concordance** (`filtering_haploid_mappping.py`). Joins calls to the haploid alignment on `readID`, drops reads that map to multiple haploid locations, flips the read coordinate when the diploid and haploid strands disagree (`read_length - read_position - 1`), then walks the haploid CIGAR to check whether an InDel of the **same type and same size** appears at the expected read position. Adds `haploid_check`, `distance_to_indel_in_hap`, `softclip_haploid`, `other_ops`, plus `species`/`sample`.

**Step 7 — merge** per-contig results into `output/final/`.

**Step 8 — final filtering** (`filtering_small_indels.py`, `filtering_opportunities.py`) with `--haploid-check True --min-quality 93`.

---

## Workflow C — `repetitive_indel_call.py`: tandem repeat InDels

The mirror image of Workflow B: it scans **inside** M0/M1 intervals and runs the whole chain twice per sample, once per model.

The two substantive differences (`utilities.tandem_repeat_indels()`):

- **Size rule instead of a size cap.** Any InDel is kept if `abs(indel) % k == 0`, i.e. whole-motif expansions and contractions of any length, rather than 1–3 bp.
- **No homopolymer filter**, because at k=1 homopolymer slippage *is* the signal being measured.

Six extra columns are carried through: `canonical_motif, k, copies, trf_purity, model_purity, model`. Steps 2, 5, 7 and 8 reuse Workflow B's scripts unchanged. This workflow additionally calls `plot_indels.R` to produce a diagnostics PDF.

---

## Filtering strategy

Filters are applied at three different stages, which is worth keeping in mind when reading the code:

| Filter | Where applied |
|---|---|
| Coverage 0.5×–2× sample average | during pileup, `utilities.py` |
| Repeat / unresolved-region masking | during pileup, via `IntervalTree` |
| InDel size (1–3 bp, or period-multiple) | during pileup |
| Homopolymer context | during pileup (Workflow B only) |
| kmerQ ≥ 93 | at final filtering, `--min-quality 93` |
| Haploid concordance | at final filtering, `--haploid-check True` |
| Distance to read end | at final filtering, `--min-read-end-distance` (Workflow C only) |

`filtering_small_indels.py` and `filtering_opportunities.py` share a generic design: a `filter_indels(df, **filters)` helper supporting `eq/ne/gt/gte/lt/lte/in/notin`, plus `apply_filters_with_report()` which applies steps sequentially and records rows removed at each stage into a `--report` TSV. Adding a new filter means adding one `add_argument` and one `filter_steps.append(...)` line — no new filtering logic.

## Output layout

```
output/
├── tandem_repeats/{species}/{sample}/{contig}_M{0,1}.bed
├── small_indels/{species}/{sample}/{contig}_small_indels[_haploid_filtered].tsv
├── opportunities/{species}/{sample}/{contig}_opportunities.tsv
├── mutation_fastq/{species}/{sample}/{contig}_remapping.fastq
├── haploid_mapping/{species}_{sample}_haploid_locations.tsv
├── tandem_repeat_*/…                       # Workflow C equivalents
├── final/{sample}[_{model}_tr]_combined_filtered.tsv
│        {sample}_opportunities_combined_filtered.tsv
│        {sample}_filter_report.tsv
└── plots/{sample}_indel_plots.pdf
```

The two headline products are a **filtered InDel table** (numerator) and a **filtered opportunity table** (denominator); a mutation rate is the ratio of the two, stratified by context, tissue, and — for Workflow C — repeat motif and period.

## Dependencies

- **Python 3** with `pysam`, `pandas`, `intervaltree`, and `gwf`
- **R** with `ggplot2`, `dplyr`, `tidyr`
- **External binaries**: `mreps`, `trf`, `pbmm2`, `samtools`, `bedtools`, `faSplit`, plus heavy use of `awk`

There is no dependency manifest in the repository — no `requirements.txt`, `pyproject.toml`, `environment.yml`, `Makefile` or `Dockerfile`. Versions are unpinned and the list above is reconstructed from imports and from `README.md`.

## Running order

```bash
gwf -f repetitiveDNA.py run          # 1. build M0/M1 BEDs (prerequisite for workflow C)
gwf -f indel_mutation_call.py run    # 2. non-repetitive InDels
gwf -f repetitive_indel_call.py run  # 3. tandem repeat InDels (can run parallel with 2)
```
