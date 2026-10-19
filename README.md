# detecting_InDels

Finds insertions and deletions (InDels) carried by single PacBio HiFi reads, and
turns them into InDel mutation rates inside and outside tandem repeats, across
primate species.

Each individual's HiFi reads are aligned back to the diploid assembly built from
those same reads. Where a single read disagrees with its own assembly, the read
may carry a mutation that arose in the cell it came from. A sequencing error in
one read looks exactly like a mutation in one read, so most of this pipeline is
filtering. Every set of calls comes with an **opportunity table**, counted in the
same pass over the reads under the same rules. A rate is calls divided by
opportunities, for any stratum: InDel size and sequence, repeat motif and period,
or tissue.

## Contents

- [Overview](#overview)
- [Getting started](#getting-started)
- [Workflow 1: precompute](#workflow-1-precompute)
- [Workflow 2: repeat annotation](#workflow-2-repeat-annotation)
- [Workflow 3: InDels outside repeats](#workflow-3-indels-outside-repeats)
- [Workflow 4: InDels inside repeats](#workflow-4-indels-inside-repeats)
- [Filtering in detail](#filtering-in-detail)
  - [Where each filter acts](#where-each-filter-acts)
  - [The anchor base](#the-anchor-base) ·
    [Coverage](#coverage) ·
    [Mask and tracts](#mask-and-tracts) ·
    [InDel size](#indel-size) ·
    [Homopolymer rule](#homopolymer-rule)
  - [kmerQ](#kmerq) ·
    [Haploid concordance](#haploid-concordance) ·
    [Distance to the read end](#distance-to-the-read-end) ·
    [The unique flag](#the-unique-flag)
  - [Final filtering and the filter report](#final-filtering-and-the-filter-report) ·
    [What the opportunities go through](#what-the-opportunities-go-through) ·
    [Filters in the notebooks](#filters-in-the-notebooks)
- [Tissue and cell-type labels](#tissue-and-cell-type-labels)
- [Outputs](#outputs)
- [Analysis notebooks](#analysis-notebooks)
- [Tests](#tests)
- [Repository layout](#repository-layout)

---

## Overview

### How it works

1. **Align.** Outside this repository, each sample's HiFi reads are aligned to
   its own diploid assembly. The BAM keeps PacBio's per-read tags.
2. **Scan.** The InDel callers walk that BAM one position at a time. A read with
   an insertion or deletion relative to the assembly is a candidate call.
3. **Count opportunities.** In the same pass, every read position where an InDel
   could have been seen, under the same filters, is counted. This is the
   denominator of the rate.
4. **Confirm.** Reads carrying a call are aligned again, to the haploid assembly.
   A call is kept only if the same InDel shows up there.
5. **Filter.** Thresholds on base quality and on position in the read remove
   what is most likely sequencing error.

Tandem repeats mutate by replication slippage, at rates far above the rest of the
genome, so they are measured separately. The pipeline annotates repeats first,
masks them out for one caller, and scans only them with a second caller that has
its own size rule.

### The four workflows

| # | Workflow | What it does | Needs | Writes |
|---|---|---|---|---|
| 1 | [precompute.py](workflows/precompute.py) | Per-contig FASTAs, mean coverage, unresolved-haplotype regions, read → tissue table | the inputs | `output/precomute/` |
| 2 | [repetitive_dna.py](workflows/repetitive_dna.py) | Tandem repeat tracts with periods 2–6, under three models (M0, M1, M2) | 1 | `output/repetitive_dna/`, `results/repetitive_dna/` |
| 3 | [indel_mutation_call.py](workflows/indel_mutation_call.py) | 1–3 bp InDels outside every repeat tract and every unresolved region | 1, 2 | `results/nonrepetitive_indels/` |
| 4 | [repetitive_indel_call.py](workflows/repetitive_indel_call.py) | InDels inside the tracts of each model, sized a whole number of repeat units | 1, 2 | `results/repetitive_indels/` |

```
precompute.py ──> repetitive_dna.py ──┬──> indel_mutation_call.py     InDels outside repeats
                                      └──> repetitive_indel_call.py   InDels inside repeats
```

### What the InDel callers share

Workflows 3 and 4 have the same shape:

```
10 Mb window of the BAM
 ├─ calls ────────────┐
 ├─ reads with calls ─┴─> remap to haploid assembly ─> haploid check ─┐
 └─ opportunities ────────────────────────────────────────────────────┴─> final filters
```

- **One pass per window.** Each contig is cut into 10 Mb windows. One job reads
  a window once and writes the calls, the opportunity counts and a FASTQ of the
  reads carrying a call. Every position belongs to exactly one window, so no call
  or opportunity is counted twice. The windows are then joined per contig: calls
  concatenated, opportunity counts summed, FASTQ records kept once per read.
- **Every filter is explicit.** The pileup runs with pysam's own read and base
  filters turned off (`stepper="nofilter"`, `min_base_quality=0`,
  `min_mapping_quality=0`), so everything that removes data is in this
  repository's code.
- **Numerator and denominator match.** Calls and opportunities come from the same
  pileup and pass the same coverage, region and context rules. The same
  base-quality threshold is applied to both at the end.
- **Long contigs only.** Every workflow processes only the contigs longer than
  `min_contig_length` (1 Mb) in the assembly's `.fai`.
- **Crash-safe outputs.** Scan and repeat-annotation jobs write to a `.tmp` file
  and rename it when done, so a job that dies never leaves a finished-looking
  file for GWF to skip.

---

## Getting started

### Requirements

| Kind | What |
|---|---|
| Workflow engine | [GWF](https://gwf.app/) ≥ 1.7, submitting to Slurm |
| Python | pysam, pandas, intervaltree, pyyaml; pytest for the tests |
| Binaries | samtools, bedtools, pbmm2, faSplit, TRF 4.10, mreps |
| R | ggplot2, dplyr, tidyr for the diagnostic plots; also readr and the IRkernel Jupyter kernel for the notebooks |

On the cluster, [environment.yml](environment.yml) provides everything except
mreps:

```bash
conda env create -f environment.yml
conda activate detecting-indels
```

mreps is not in bioconda. Build it from source and point `tools.mreps` in
config.yaml at the binary (more on the [mreps webpage](http://mreps.univ-mlv.fr/)):

```bash
git clone https://github.com/gregorykucherov/mreps.git
cd mreps
git checkout 5e0ec29    # mreps 2.6, the version used here
make
./mreps -res 1 -fasta test/sequence.fasta   # quick check
```

TRF must support `-ngs` and `-l`, as TRF 4.10 does. `tools.trf` points at the
binary, so a build from [source](https://github.com/Benson-Genomics-Lab/TRF)
works as well as the bioconda package:

```bash
git clone https://github.com/Benson-Genomics-Lab/TRF.git
cd TRF
git checkout 355c1f9    # the version used here: 13 commits after tag v4.10.0-rc.2
mkdir build && cd build
../configure && make    # the binary is build/src/trf
```

The local source copies in `TRF/` and `mreps/` are at those two commits. The
cluster binaries in config.yaml are built the same way, but which commits they
were built from hasn't been checked; `mreps -version` and `trf -v` on the
cluster would confirm it.

For local work (tests and notebooks), [requirements.txt](requirements.txt) is
enough:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

### Configuration

Every workflow reads [config.yaml](config.yaml) from the directory gwf is started
in.

| Key | Meaning |
|---|---|
| `samples` | species → sample names. Every workflow runs every sample listed here. |
| `tissue_samples` | Samples with a read → tissue label file. Reads of the other samples get tissue `NA`. |
| `methylation_samples` | Samples with per-read cell-type likelihoods. Both InDel workflows add germline/somatic columns for these. |
| `alignment_overrides` | species → sample → BAM, for samples whose alignment is not at the standard path. |
| `assembly_type` | `diploid`: the assembly the reads are aligned to. |
| `min_contig_length` | Contigs of this length or shorter are skipped. |
| `paths.fasta`, `paths.bam` | Roots of the assemblies and alignments (see [Inputs](#inputs)). |
| `paths.read_labels` | Root of the read → tissue label files. |
| `paths.classification` | Root of the methylation likelihood tables. |
| `paths.output`, `paths.results`, `paths.plots` | Where intermediates, final tables and PDFs go, relative to the repository root. |
| `tools.mreps`, `tools.trf` | The two repeat finders. |
| `layout` | Where every file lives: one path template per file, filled from `paths`, `assembly_type`, species, sample and, per file, contig, model, window start and end, or period k. The workflows build every path from these, through [workflows/pipeline.py](workflows/pipeline.py). |
| `precompute`, `repeats`, `scan`, `final_filters`, `classification` | The workflows' parameters, below. The scripts take them as arguments. |

No workflow reads `sperm_samples` or `sample_age`; the notebooks do.

Parameters:

| Parameter | Value | Key |
|---|---|---|
| Unresolved-haplotype scan window | 40 Mb | `precompute.unresolved_window_size` |
| Unresolved-haplotype evidence | ≥ 2 reads with the same base at quality ≥ 40 | `precompute.unresolved_min_reads`, `precompute.unresolved_min_base_quality` |
| Repeat models | M0, M1, M2 | `repeats.models` |
| Repeat periods | 2, 3, 4, 5, 6 | `repeats.periods` |
| mreps minimum exponent (M0, M2) | 3 copies | `repeats.mreps_min_exponent` |
| mreps resolution (M2) | 1 | `repeats.m2_resolution` |
| TRF scoring (M1) | `2 7 7 80 10 50`, maximum period 6 | `repeats.trf_params` |
| TRF window and overlap (M1) | 5 Mb, plus 1 Mb on each side | `repeats.trf_window_size`, `repeats.trf_overlap` |
| InDel scan window | 10 Mb | `scan.window_size` |
| Coverage window | 0.5×–2× the sample mean | `scan.coverage_band` |
| InDel sizes outside repeats (workflow 3) | 1–3 bp | `scan.indel_sizes` |
| Minimum kmerQ | 93 | `final_filters.min_kmerq` |
| Minimum distance to the read end (workflow 4) | 100 bp | `final_filters.min_read_end_distance` |
| Germline and somatic cell types; the 14 score bins | `Sperm.*` against six somatic types; edges −10 to 10 | `classification.germline_prefix`, `.somatic_types`, `.bin_edges`, `.bin_labels` |
| Slurm account and per-target cores, memory, walltime | `mutationalscanning` | every workflow |

### Inputs

Per sample, with `{fasta}`, `{bam}`, `{read_labels}` and `{classification}` the
roots under `paths`. All of these are made outside this repository.

| Input | Path | Used by |
|---|---|---|
| Diploid assembly, with `.fai` | `{fasta}/{species}/{sample}/diploid_assembly/{sample}_diploid.fa` | all |
| Haploid assembly | `{fasta}/{species}/{sample}/haploid_assembly/{sample}_haploid.fa` | 3, 4 |
| Reads aligned to the diploid assembly, indexed | `{bam}/{species}/{sample}/kinetics/{sample}_kinetics_diploid.bam`, or the `alignment_overrides` entry | 1, 3, 4 |
| Read → tissue labels | `{read_labels}/{species}/{sample}/{sample}_read_labels.txt` | 1, `tissue_samples` only |
| Per-read cell-type likelihoods | `{classification}/{species}/{sample}/{sample}/05_ReadsClassification/{sample}_All_Reads_Likelihood_merged.tsv` | 3, 4, `methylation_samples` only |

The callers record PacBio's per-read tags when the BAM has them: `np` and `ec`,
and the per-base subread tags `sa`, `sm` and `sx`. A missing tag is written as
`NA`. Nothing is filtered on these tags.

### Running the pipeline

GWF submits to Slurm once its backend is set (`gwf config set backend slurm`).
Then, from the repository root, in this order:

```bash
gwf -f workflows/precompute.py run
gwf -f workflows/repetitive_dna.py run
gwf -f workflows/indel_mutation_call.py run
gwf -f workflows/repetitive_indel_call.py run
```

- Start each workflow once the one before it has finished
  (`gwf -f workflows/<name>.py status`). GWF tracks dependencies within one
  workflow file, not between files. Workflows 3 and 4 can run at the same time.
- Run from the repository root. The workflows open `config.yaml` there and
  resolve `output/`, `results/`, `plots/` and `scripts/` from it.
- The target names in the step tables below are what `gwf status` and
  `gwf cancel` take. Wildcards select groups, for example
  `gwf -f workflows/repetitive_indel_call.py status 'TRScan_*'`.

To add a sample, list it under `samples` (and under `tissue_samples`,
`methylation_samples` or `alignment_overrides` if they apply), check that its
inputs are at the paths above, and run the four workflows again. GWF only runs
targets whose outputs are missing or out of date.

---

## Workflow 1: precompute

[workflows/precompute.py](workflows/precompute.py) prepares what the other
workflows share. Outputs go to `output/precomute/{species}/{sample}/` (the
directory is spelled `precomute` in the code).

| Step | Target | Runs per | Output |
|---|---|---|---|
| 1. Split the assembly | `faSplit_{species}_{sample}_diploid` | sample | `contigs/{contig}.fa` |
| 2. Mean coverage | `average_coverage_{species}_{sample}` | sample | `long_contigs.bed`, `avg_coverage.txt` |
| 3a. Unresolved haplotypes | `unresolved_haplotype_{species}_{sample}_{contig}_{start}_diploid` | 40 Mb window | `unresolved_haplotypes/chunks/{contig}_{start}_{end}.bed` |
| 3b. Join the windows | `merge_unresolved_{species}_{sample}_{contig}_diploid` | contig | `unresolved_haplotypes/{contig}.bed` |
| 4. Tissue table | `sorting_tissue_classification_{species}_{sample}_diploid` | sample in `tissue_samples` | `sorted_tissue_classification.tsv` |

**Split.** `faSplit byname` writes one FASTA per contig, for mreps in workflow 2.

**Mean coverage.** `samtools depth` runs over every contig longer than
`min_contig_length`, and the depths it reports are averaged (positions without
reads are not reported, so they do not lower the mean). The single number in
`avg_coverage.txt` sets the [coverage window](#coverage) of both callers.

**Unresolved haplotypes.** Where the assembly did not separate the two
haplotypes, or collapsed two copies of a duplicated sequence, reads from the
missing copy differ from the assembly at every site where the copies differ.
Those differences would be called as mutations. `resolve_haplotypes()` in
[utilities.py](scripts/utilities.py) finds such regions:

1. At every position of the window, it takes the reads whose base differs from
   the assembly and has base quality ≥ 40.
2. If at least two of them carry the same base, it marks a region from the
   leftmost start to the rightmost end of those reads' alignments.
3. A region that overlaps the previous one extends it.

A region can reach past its window, so the windows' BEDs are joined with
`bedtools merge`. The result is a headerless BED (`contig, start, end`) per
contig. It is part of workflow 3's mask; workflow 4 does not use it.

**Tissue table.** The read label file (read ID and tissue, whitespace-separated)
becomes a sorted two-column TSV with a `readID  tissue` header.

---

## Workflow 2: repeat annotation

[workflows/repetitive_dna.py](workflows/repetitive_dna.py) annotates tandem
repeats with periods k = 2 to 6 under three models. Homopolymers (k = 1) are not
annotated; workflow 3 handles them with its [homopolymer rule](#homopolymer-rule).
Outputs go to `output/repetitive_dna/{species}/{sample}/`.

The tracts serve both InDel workflows. Workflow 3 masks out the union of all
three models; workflow 4 scans each model's tracts on their own.

| Model | Tool | A tract is | Command | `trf_purity` | `model_purity` |
|---|---|---|---|---|---|
| M0 | mreps | an exact repeat: every copy identical, at least 3 copies | `mreps -minperiod k -maxperiod k -exp 3 -allowsmall` | 100 | 100 |
| M1 | TRF | an approximate repeat: copies may differ by substitutions and indels, and the alignment must score ≥ 50 (about 25 bp for a pure tract) | `trf window.fa 2 7 7 80 10 50 6 -ngs -h -l 7` | TRF's percent matches | Hamming match % |
| M2 | mreps | at least 3 copies, with up to one substitution between neighbouring copies | `mreps -minperiod k -maxperiod k -exp 3 -res 1 -allowsmall` | mreps' (1 − error rate) × 100 | Hamming match % |

| Step | Target | Runs per | Output |
|---|---|---|---|
| 1a. M0 tracts | `M0_mreps_k{k}_{contig}_{species}_{sample}_diploid` | contig × k | `per_k/{contig}_M0_k{k}.bed` |
| 1a. M2 tracts | `M2_mreps_k{k}_{contig}_{species}_{sample}_diploid` | contig × k | `per_k/{contig}_M2_k{k}.bed` |
| 1a. M1 tracts | `M1_trf_{contig}_{start}_{species}_{sample}_diploid` | 5 Mb window | `per_k/{contig}_M1_k{k}_{start}.bed`, one per k |
| 1b. Resolve overlaps | `{model}_merge_k_{contig}_{species}_{sample}_diploid` | contig × model | `{contig}_{model}.bed` |
| 2. Sample BED and summary | `cat_{model}_{species}_{sample}_diploid` | sample × model | `{model}.bed`, `results/repetitive_dna/{species}/{sample}/{model}_summary.tsv` |

**mreps (M0, M2).** mreps runs once per contig and period.
[parse_mreps_canonical.py](scripts/parse_mreps_canonical.py) turns its report into
the repeat BED:

- coordinates go from 1-based inclusive to 0-based half-open;
- `copies` is mreps' exponent, the tract length over the period, so it can be
  fractional;
- the motif is the consensus of the full copies (the most common base at each
  position, ties going to the first copy), since in an impure tract any copy can
  carry the mismatch;
- `trf_purity` holds mreps' own measure, (1 − error rate) × 100, and
  `model_purity` the Hamming measure described under *Purity* below.

**TRF (M1).** TRF's runtime explodes in long satellite arrays, and a whole
chromosome can run for days, so M1 runs per window:

- each contig is cut into 5 Mb core windows, and each window is fetched with
  1 Mb extra on both sides;
- a window keeps only the tracts that start in its core, so no tract is reported
  twice, and a tract crossing a core boundary is still found whole if it is
  shorter than the overlap (a longer one can be cut in two);
- TRF's maximum repeat length (`-l`, in millions of bp) is set to cover the whole
  fetched window, so no array inside it makes TRF restart its alignment;
- TRF runs once per window with maximum period 6, and its output is split by
  period.

[parse_trf.py](scripts/parse_trf.py) shifts the coordinates back onto the contig
and makes them 0-based half-open. `copies` is TRF's copy number, and `trf_purity`
its percent matches, which comes from an alignment and so forgives indels
between copies.

**Purity.** `model_purity` compares every base of the tract with the consensus
motif tiled from the tract's first base, and reports the percentage that match.
It has no notion of gaps: an indel inside the tract shifts the phase and lowers
it for every base after. For M0 both measures are 100.

**Canonical motif.** A motif is labelled by the smallest of its rotations and the
rotations of its reverse complement, so `AC`, `CA`, `GT` and `TG` are all `AC`:
the same repeat, read from either strand and starting anywhere. Both parsers use
`canonicalise()` from parse_trf.py.

**Overlap resolution.**
[resolve_repeat_overlaps.py](scripts/resolve_repeat_overlaps.py) leaves every
base in at most one tract per model:

- periods are handled smallest first, and a tract that overlaps one already kept
  is dropped whole;
- within one period, overlapping tracts are resolved longest first (ties go to
  the one listed first);
- tracts that only touch do not overlap.

Keeping the smallest period keeps a repeat's fundamental period: TRF often
reports a repeat again at multiples of its period, such as `(AT)n` at periods 4
and 6. A dropped tract is not trimmed, so its bases outside the kept tract lose
that model's annotation. Workflow 3 builds its mask from these resolved BEDs.

**Sample BED and summary.** The contigs' resolved BEDs are concatenated into
`{model}.bed`. The summary counts tracts per distinct combination of
`canonical_motif`, `k`, `copies`, `trf_purity` and `model_purity`, which keeps a
whole genome small enough to plot locally.

Every repeat BED has these 11 tab-separated columns and a header line:

| Column | Meaning |
|---|---|
| `contig`, `start`, `end` | tract coordinates, 0-based half-open |
| `canonical_motif` | motif label, see *Canonical motif* above |
| `k` | period |
| `copies` | copy number: tract length / period for mreps, TRF's estimate for M1 |
| `trf_purity` | TRF's percent matches (M1); mreps' (1 − error rate) × 100 (M0, M2) |
| `model_purity` | Hamming match % against the tiled consensus |
| `model` | `M0`, `M1` or `M2` |
| `species`, `sample` | |

---

## Workflow 3: InDels outside repeats

[workflows/indel_mutation_call.py](workflows/indel_mutation_call.py) calls
1–3 bp InDels at positions outside every repeat tract and every
unresolved-haplotype region. Intermediates go to
`output/nonrepetitive_indels/{species}/{sample}/`.

| Step | Target | Runs per | Output |
|---|---|---|---|
| 1. Build the mask | `non_rep_Mask_{species}_{sample}_diploid` | sample | `mask/{contig}.bed` |
| 2. Scan | `non_rep_Scan_{contig}_{start}_{species}_{sample}_diploid` | 10 Mb window | `calls/chunks/`, `opportunities/chunks/`, `fastq/chunks/` |
| 3. Join the windows | `non_rep_MergeWindows_{species}_{sample}_diploid` | sample | `calls/{contig}.tsv`, `opportunities/{contig}.tsv`, `fastq/{contig}.fastq` |
| 4. Remap to haploid | `non_rep_HaploidMapping_{species}_{sample}_diploid` | sample | `haploid/all.fastq`, `haploid/remapped.bam`, `haploid/locations.tsv` |
| 5. Haploid check | `non_rep_HaploidFilter_{species}_{sample}_diploid` | sample | `merged/indels.tsv` |
| 6. Merge opportunities, plot | `non_rep_MergeAndPlot_{species}_{sample}_diploid` | sample | `merged/opportunities.tsv`, `plots/nonrepetitive_indels/{species}/{sample}/plots.pdf` |
| 7. Final filtering | `non_rep_FinalFiltering_{species}_{sample}_diploid` | sample | `indels.tsv`, `opportunities.tsv`, `filter_report.tsv` in `results/nonrepetitive_indels/{species}/{sample}/` |
| 8. Germline/somatic | `non_rep_Classify_{species}_{sample}_diploid` | sample in `methylation_samples` | `indels_classified.tsv`, next to `indels.tsv` |

1. **Mask.** Per contig, the tracts of all three models (after overlap
   resolution) and the unresolved-haplotype regions, merged with
   `bedtools merge`.
2. **Scan.** [nonrepetitive_scan.py](scripts/nonrepetitive_scan.py) runs
   `nonrepetitive_scan()` from [utilities.py](scripts/utilities.py). At each
   position of the window that passes the coverage window and lies outside the
   mask, every read with a base there:
   - adds its opportunities: one per InDel size (1, 2 and 3 bp), under the
     canonical k-mer of that size starting at the read position, with its
     tissue, mapping quality and kmerQ (a read's first and last base are not
     counted);
   - becomes a call if it carries a 1–3 bp InDel anchored there that passes the
     homopolymer rule, and then goes to the FASTQ.
3. **Join.** Per contig, the windows' calls are concatenated, their opportunity
   counts summed per value of every other column, and their FASTQs
   concatenated keeping each read once.
4. **Remap.** All reads carrying a call are aligned to the haploid assembly with
   `pbmm2 align --preset HIFI --sort`. The alignments are reduced to read ID,
   flag, contig, position, mapping quality and CIGAR, sorted by read ID.
5. **Haploid check.** [filtering_haploid_mappping.py](scripts/filtering_haploid_mappping.py)
   checks every call against its read's haploid alignment (see
   [Haploid concordance](#haploid-concordance)) and writes one table for the
   sample, with `species` and `sample` columns added.
6. **Opportunities and plots.** The contigs' opportunities are concatenated, with
   `species` and `sample` added. [plot_indels.R](scripts/plot_indels.R) draws
   diagnostics of the calls before final filtering: kmerQ, the 20 most common
   InDel sequences, distance to the read end, the unique flag, the haploid
   check, tissue, kmerQ against haploid distance, and subread coverage.
7. **Final filtering.** Calls must pass the haploid check and have kmerQ ≥ 93.
   Opportunities must have kmerQ ≥ 93.
8. **Germline/somatic.** See [Tissue and cell-type labels](#tissue-and-cell-type-labels).

---

## Workflow 4: InDels inside repeats

[workflows/repetitive_indel_call.py](workflows/repetitive_indel_call.py) calls, for
each repeat model, the InDels anchored inside that model's tracts whose length is
a whole number of repeat units (a multiple of the tract's k), of any size.
Intermediates go to `output/repetitive_indels/{species}/{sample}/`.

| Step | Target | Runs per | Output |
|---|---|---|---|
| 1. Scan | `TRScan_{contig}_{start}_{species}_{sample}` | 10 Mb window | `{model}/calls/chunks/`, `{model}/opportunities/chunks/`, `fastq/chunks/` |
| 2. Join the windows | `TRMergeWindows_{contig}_{species}_{sample}` | contig | `{model}/calls/{contig}.tsv`, `{model}/opportunities/{contig}.tsv`, `fastq/{contig}.fastq` |
| 3. Remap to haploid | `TRHaploidMapping_{species}_{sample}` | sample | `haploid/all.fastq`, `haploid/remapped.bam`, `haploid/locations.tsv` |
| 4. Haploid check | `TRHaploidFilter_{species}_{sample}_{model}` | sample × model | `{model}/merged/indels.tsv` |
| 5. Merge opportunities, plot | `TRMergeAndPlot_{species}_{sample}_{model}` | sample × model | `{model}/merged/opportunities.tsv`, `plots/repetitive_indels/{species}/{sample}/{model}/plots.pdf` |
| 6. Final filtering | `TRFinalFiltering_{species}_{sample}_{model}` | sample × model | `indels.tsv`, `opportunities.tsv`, `filter_report.tsv` in `results/repetitive_indels/{species}/{sample}/{model}/` |
| 7. Germline/somatic | `TRClassify_{species}_{sample}_{model}` | sample × model, in `methylation_samples` | `indels_classified.tsv`, next to `indels.tsv` |

The steps mirror workflow 3, with two savings. One scan of a window serves all
three models ([tandem_repeat_scan.py](scripts/tandem_repeat_scan.py) runs
`tandem_repeat_scan()`), and the reads carrying a call in any model are remapped
to the haploid assembly once per sample. Each model's output is the same as
scanning that model alone, which [tests/test_integration.py](tests/test_integration.py)
checks.

At each position that passes the coverage window and lies inside a tract of a
model, every read with a base there:

- adds one opportunity under the tract's motif and k, with its tissue, mapping
  quality and kmerQ (a read's first and last base are not counted);
- becomes a call for that model if it carries an InDel anchored there whose
  length is a multiple of the tract's k. The call carries the tract's
  `canonical_motif`, `k`, `copies`, `trf_purity`, `model_purity` and `model`.

How the two callers differ:

| | Workflow 3, outside repeats | Workflow 4, inside repeats |
|---|---|---|
| Positions scanned | outside the mask | inside a tract of the model |
| Unresolved-haplotype regions | masked | not masked |
| InDel size | 1–3 bp | any multiple of the tract's k |
| Homopolymer rule | applied | not applied |
| `canonical_indel` | the smaller of the sequence and its reverse complement | the smallest rotation of the sequence or its reverse complement, like the tract motif |
| Opportunity unit | read position × InDel size, by canonical k-mer | read base inside a tract, by tract motif and k |
| Final filters | haploid check, kmerQ ≥ 93 | haploid check, kmerQ ≥ 93, ≥ 100 bp from the read end |
| Runs | once per sample | once per model per sample |

---

## Filtering in detail

Filters act at four points:

- **During the scan**, before anything is written. These decide which positions
  and reads count at all, and they apply to calls and opportunities alike.
- **At the haploid check**, which annotates every call and drops calls whose read
  has no usable haploid alignment.
- **At final filtering**, which applies thresholds to recorded columns and
  reports what each one removed.
- **In the analysis notebooks**, which apply a few more to the final tables.

### Where each filter acts

| Filter | Stage | Workflow 3 | Workflow 4 | Opportunities |
|---|---|---|---|---|
| Contig longer than 1 Mb | setup | ✓ | ✓ | ✓ |
| Coverage within 0.5×–2× the sample mean | scan | ✓ | ✓ | ✓ |
| Position outside the mask / inside a tract | scan | outside | inside | ✓ |
| InDel size | scan | 1–3 bp | multiple of k | – |
| Homopolymer rule | scan | ✓ | – | ✓ in workflow 3 |
| Read has a usable haploid alignment | haploid check | ✓ | ✓ | – |
| Same InDel in the haploid alignment | final | ✓ | ✓ | – |
| kmerQ ≥ 93 | final | ✓ | ✓ | ✓ |
| At least 100 bp from a read end | final | – | ✓ | – |
| No other InDel at the position (`unique`) | notebooks | ✓ | ✓ | – |

Nothing is filtered on read flags (secondary, supplementary and duplicate
alignments all enter the pileup), on mapping quality, or on base quality outside
the kmerQ window. Mapping quality to the diploid assembly is kept in the
opportunity tables (`mapQ`) but not in the calls, which carry the haploid
mapping quality (`map_mapq`) instead.

### The anchor base

pysam reports an InDel on the reference base just before the inserted or deleted
sequence: the anchor. Every positional test (coverage, mask, tract) is made at
the anchor, and a call's `position` column is the anchor's 1-based coordinate on
the contig.

No padding is added around masked tracts. If the aligner shifts an InDel in a
repeat to its leftmost equivalent position, the anchor is the base before the
tract. Workflow 4 then does not see the InDel, while workflow 3 sees it just
outside its mask and calls it if it is 1–3 bp and passes the homopolymer rule.
[repetitive_mutation_rate.ipynb](plotting/repetitive_mutation_rate.ipynb)
discusses what this does to the M0 calls.

### Coverage

A position is used only if the number of reads covering it (pysam's
`nsegments`, which includes reads with a deletion there) is between half and
twice the sample's mean coverage. This removes positions with too many reads,
such as collapsed duplications, and too few, such as poorly assembled or poorly
mapped ones. The test is made once per position, so calls and opportunities
share it.

### Mask and tracts

Workflow 3 skips every position inside its mask: the union of the M0, M1 and M2
tracts and the unresolved-haplotype regions. Workflow 4 uses only positions
inside a tract of the model being called, and does not mask unresolved
haplotypes. After overlap resolution a base lies in at most one tract per model,
so each position has a single tract, motif and k. (The code's tie-breaks, the
largest k for a call and the smallest for an opportunity, only matter for BEDs
that were not resolved.)

### InDel size

Workflow 3 keeps InDels of 1, 2 or 3 bp. Workflow 4 keeps an InDel if its length
is a multiple of the tract's k, with no upper limit: 2 or 4 bp in an `AC` tract,
but not 3 bp.

### Homopolymer rule

Homopolymers are not annotated as repeats, so workflow 3 removes them itself. It
takes the InDel's sequence with one flanking base on each side (the anchor and
the next base), from the assembly for a deletion and from the read for an
insertion, and drops the call if that run changes base at most once.

| Event | Run checked | Base changes | Result |
|---|---|---|---|
| deletion of `A` from `TAG` | `TAG` | 2 | called |
| deletion of one `A` from `TAAG` | `TAA` | 1 | dropped |
| insertion of `C` between `A` and `G` | `ACG` | 2 | called |
| insertion of `A` between `A` and `G` | `AAG` | 1 | dropped |
| deletion of `AC` from `AACC` | `AACC` | 1 | dropped |
| deletion of `AAA` from `TAAAT` | `TAAAT` | 2 | called |

The rule counts base changes, not run length, so a whole homopolymer deleted from
between two different bases (the last row) is called.

The opportunity counter applies the same rule to the k-mer at each read position
together with the base on each side of it. Workflow 4 has no homopolymer rule:
its periods start at 2, and a 1 bp InDel is never a multiple of k.

### kmerQ

`kmerQ` is the lowest base quality in a four-base window of the read around the
event:

```
deletion of CG      assembly  G A T C G T A C
                    read      G A T - - T A C
                    kmerQ       ^ ^     ^ ^       two read bases on each side of the gap

insertion of CG     assembly  G A T - - T A C
                    read      G A T C G T A C
                    kmerQ         ^ ^ ^ ^         anchor, first and last inserted base, next base

opportunity, 2 bp   read      G A T C G T A C     k-mer CG
                    kmerQ         ^ ^ ^ ^         base before, first and last k-mer base, base after
```

For a 1 bp event the window holds three distinct bases (the middle one counts
twice). For a 3 bp event the middle base of the sequence is not in it. Inside
repeats, the opportunity window uses the tract's k as the size. A deletion's
window holds only flanking bases, since the deleted bases are not in the read.

HiFi base qualities are capped at 93, so kmerQ ≥ 93 means every base in the
window is at the maximum. The threshold is applied at final filtering, to calls
(`--min-quality 93`) and opportunities (`--min-kmerq 93`) alike. The merged
tables in `output/` keep every value (opportunities are counted per kmerQ from
0 to 93), so another threshold can be applied to them without scanning again.

### Haploid concordance

A real InDel is part of the read, so it should appear wherever the read is
aligned. One that appears only in the alignment to the diploid assembly is more
likely an artefact of that alignment or assembly. Each call is therefore checked
against an independent alignment of its read to the haploid assembly, in
[filtering_haploid_mappping.py](scripts/filtering_haploid_mappping.py):

1. Every read carrying a call (before final filtering) is written to FASTQ in its
   original orientation and aligned with `pbmm2 align --preset HIFI`.
2. Calls are joined to their read's alignment records. The read's first record in
   the sorted table is used, which is its primary alignment when it has one. A
   call whose record is not a plain forward or reverse alignment (flag 0 or 16),
   because the read did not align or aligned only as secondary or supplementary,
   is dropped here. It does not appear in the filter report.
3. If the read lies on opposite strands in the two alignments, the call's read
   position is mirrored (`read_length − start − 1`).
4. The haploid CIGAR is walked from the start of the read, counting read bases
   through `M`, `=`, `X` and `S` operations. The first `I` (for an insertion) or
   `D` (for a deletion) at or after the call's read position decides:
   `haploid_check` is True if its length equals the call's, and False otherwise.
   With no such operation it is False.

Final filtering keeps calls with `haploid_check` True. The distance between the
call's read position and the matched haploid event (`distance_to_indel_in_hap`)
is recorded but not limited, so a call can pass on an InDel of the same size far
along the read. `max_hap_distance` in the notebooks can limit it.

### Distance to the read end

`ToreadEnd` is the distance from the anchor to the nearer end of the read's
alignment to the diploid assembly, in assembly bases (soft-clipped bases are not
counted). Workflow 4 requires at least 100 bp. Workflow 3 records it without
filtering; [nonrepetitive_indels.ipynb](plotting/nonrepetitive_indels.ipynb) can
apply a threshold (`min_read_end_distance`).

### The unique flag

During the scan, before any other filter, each position counts the reads that
have an InDel of any size anchored there or a deletion spanning it. A call is
`unique` if at most one read qualifies, and `non-unique` otherwise. The
workflows keep both. The notebooks drop non-unique calls by default
(`unique_only`), and `filtering_small_indels.py --unique-only` can do the same.

### Final filtering and the filter report

[filtering_small_indels.py](scripts/filtering_small_indels.py) filters the calls
and [filtering_opportunities.py](scripts/filtering_opportunities.py) the
opportunities. Each option adds one step, the steps run in a fixed order (the
order of the table below), and `--report` writes how many rows each step removed.

| Option (calls) | Keeps rows where | Used by |
|---|---|---|
| `--haploid-check True` | `haploid_check` is True | workflows 3 and 4 |
| `--softclip True` or `False` | `softclip_haploid` equals the value | – |
| `--min-distance N`, `--max-distance N` | `distance_to_indel_in_hap` ≥ N, ≤ N | – |
| `--indel-type insertion` or `deletion` | `indel_type` equals the value | – |
| `--tissue a,b` | `tissue` is in the list | – |
| `--unique-only` | `unique` is `unique` | – |
| `--min-quality N` | `kmerQ` ≥ N | workflows 3 and 4 (93) |
| `--min-read-end-distance N` | `ToreadEnd` ≥ N | workflow 4 (100) |

The opportunities take `--min-kmerq` (used, 93), `--max-kmerq`, `--min-mapq`,
`--max-mapq`, `--tissue`, `--context` and `--min-count`.

`filter_report.tsv` has one row per step (tab-separated):

```
filter                 rows_before  rows_after  rows_removed
haploid_check          …            …           …
min_kmerQ              …            …           …
min_read_end_distance  …            …           …             (workflow 4 only)
```

The first `rows_before` counts the calls that reached final filtering, so it
leaves out the calls the haploid check dropped for lack of a usable alignment.

Adding a filter takes one new option and one entry in `filter_steps`; the
comparisons (`eq`, `ne`, `gt`, `gte`, `lt`, `lte`, `in`, `notin`) already exist.
Values the callers write as `NA` stay `NA` in the output.

### What the opportunities go through

The opportunities pass the scan's filters (coverage, mask or tract, and the
homopolymer rule in workflow 3) and the kmerQ threshold. The haploid check, the
read-end distance and the unique flag remove calls only, with nothing matching
them in the denominator. Since they remove some real mutations along with the
artefacts, the rates are conservative.

### Filters in the notebooks

Set in the first cell of each rate notebook:

| Setting | Default | Notebook | Effect |
|---|---|---|---|
| `unique_only` | `TRUE` | both | drops calls flagged `non-unique` |
| `max_hap_distance` | `Inf` | both | drops calls whose haploid match is further away than this |
| `min_read_end_distance` | `0` | nonrepetitive_indels | read-end threshold for the calls outside repeats |
| `slippage_only` | `FALSE` | repetitive_mutation_rate | keeps only InDels that are whole copies of the tract's motif |

---

## Tissue and cell-type labels

**Tissue.** For samples in `tissue_samples`, every call and every opportunity
carries its read's tissue label. Both sides of the rate are split by tissue, so
per-tissue rates need no approximation. Unlabelled reads, and all reads of other
samples, are `NA`.

**Germline or somatic.** For samples in `methylation_samples`,
[join_read_classification.sh](scripts/join_read_classification.sh) adds five
columns to the final calls, in `indels_classified.tsv`. They come from the
methylation model's per-read log-likelihoods, one per cell type (higher means a
better fit):

| Column | Meaning |
|---|---|
| `max_germline` | highest likelihood among the germline cell types, the columns starting with `Sperm.` |
| `max_somatic` | highest likelihood among six somatic types: `Blood.Monocytes`, `Bladder.SmoothMuscle`, `Breast.BasalEpithelial`, `Lung.AlveolarMacrophages`, `Colon.RightEndocrine`, `Thyroid.Epithelial` |
| `classification_diff` | `max_germline − max_somatic`; positive means the read looks germline |
| `classification_bin` | the difference binned at −10, −8, −6, −4, −2, −1, 0, 1, 2, 4, 6, 8 and 10, labelled by the bin's lower edge; below −10 is labelled −12, and 10 or more is labelled 12 |
| `CpGs_onRead` | the CpGs behind the likelihoods. A read with few CpGs scores near 0 for lack of information, not because it is ambiguous. |

The germline prefix and the somatic list are the script's optional 4th and 5th
arguments. Every call is kept; a call whose read is missing from the likelihood
table gets `NA`. The opportunity tables carry no such label.

---

## Outputs

`output/` holds the large intermediates and stays on the cluster. `results/`
holds the small final tables the notebooks read, and `plots/` the diagnostic
PDFs.

### Directory layout

```
output/
├── precomute/{species}/{sample}/
│   ├── contigs/{contig}.fa
│   ├── long_contigs.bed
│   ├── avg_coverage.txt
│   ├── unresolved_haplotypes/{contig}.bed        per-window pieces in chunks/
│   └── sorted_tissue_classification.tsv          tissue_samples only
├── repetitive_dna/{species}/{sample}/
│   ├── per_k/                                    per period (and per window for M1)
│   ├── {contig}_{model}.bed                      overlaps resolved
│   └── {model}.bed                               all contigs
├── nonrepetitive_indels/{species}/{sample}/
│   ├── mask/{contig}.bed
│   ├── calls/  opportunities/  fastq/            per contig; per window in chunks/
│   ├── haploid/                                  all.fastq, remapped.bam, locations.tsv
│   └── merged/                                   indels.tsv, opportunities.tsv, before final filtering
└── repetitive_indels/{species}/{sample}/
    ├── fastq/  haploid/                          shared by the models
    └── {model}/
        ├── calls/  opportunities/
        └── merged/                               indels.tsv, opportunities.tsv

results/
├── repetitive_dna/{species}/{sample}/{model}_summary.tsv
├── nonrepetitive_indels/{species}/{sample}/
│   ├── indels.tsv
│   ├── indels_classified.tsv                     methylation_samples only
│   ├── opportunities.tsv
│   └── filter_report.tsv
└── repetitive_indels/{species}/{sample}/{model}/
    └── the same four files

plots/
├── nonrepetitive_indels/{species}/{sample}/plots.pdf
└── repetitive_indels/{species}/{sample}/{model}/plots.pdf
```

### Calls

`indels.tsv` has one row per read carrying a call. `start` and `end` are
positions in the read; the position on the contig is `position`.

| Column | Meaning |
|---|---|
| `contig` | contig of the diploid assembly |
| `start`, `end` | 1-based position in the read. For a deletion both are the anchor; for an insertion they run from the anchor to the last inserted base. |
| `ref`, `alt` | alleles including the anchor, as in a VCF. A deletion's `ref` is the anchor plus the deleted bases (from the assembly) and its `alt` the anchor; an insertion's `ref` is the anchor and its `alt` the anchor plus the inserted bases (from the read). |
| `indel_type` | `insertion` or `deletion` |
| `unique` | `unique` or `non-unique`, see [The unique flag](#the-unique-flag) |
| `indel_string` | the inserted or deleted sequence, on the assembly's forward strand |
| `canonical_indel` | `indel_string` grouped with its reverse complement (workflow 3), or with its rotations too (workflow 4) |
| `read_length` | length of the read |
| `np` | number of passes (PacBio `np` tag) |
| `ec` | effective coverage (PacBio `ec` tag) |
| `sm`, `sx` | subreads matching / mismatching at the read base after the anchor: the first inserted base, or the first base after a deletion (PacBio `sm`, `sx`) |
| `sm_pre`, `sx_pre` | the same at the anchor |
| `sa` | subread coverage at the anchor, decoded from PacBio's run-length-encoded `sa` tag |
| `ToreadEnd` | distance from the anchor to the nearer end of the read's alignment |
| `position` | 1-based coordinate of the anchor on the contig |
| `readID` | read name |
| `diploid_direction` | strand of the read's alignment to the diploid assembly |
| `canonical_motif`, `k`, `copies`, `trf_purity`, `model_purity`, `model` | workflow 4 only: the tract the call sits in |
| the keys | `tissue`, `kmerQ` and the other keys, the same columns as in the opportunity table: see [Keys](#keys) |
| `bitflag_hap`, `contig_hap`, `pos_1based_hap`, `map_mapq` | the read's alignment to the haploid assembly: flag, contig, position and mapping quality |
| `haploid_check` | whether the same InDel is in the haploid alignment, see [Haploid concordance](#haploid-concordance) |
| `distance_to_indel_in_hap` | read bases from the call's read position to the matched haploid event (0 means the same place). Without a match, the distance from a read end to the nearest event of the same type, or `NA` if there is none. |
| `softclip_haploid` | whether the haploid CIGAR has a soft clip in the part the check walked |
| `other_ops` | `HPN` if the haploid CIGAR has `H`, `P` or `N` operations, else `NA` |
| `species`, `sample` | |

### Opportunities

Workflow 3's `opportunities.tsv` has one row per contig, context, mapping
quality and value of the keys:

| Column | Meaning |
|---|---|
| `contig` | |
| `context` | canonical k-mer: the smaller of the k-mer and its reverse complement. Its length is the InDel size it counts for (1, 2 or 3), so 1 bp contexts are only `A` and `C`. |
| `mapQ` | the read's mapping quality to the diploid assembly |
| the keys | the same columns, with the same values, as at the end of each call: see [Keys](#keys) |
| `count` | read positions |
| `species`, `sample` | |

### Keys

Every feature a call might be filtered on is a key, defined once in
[scripts/strata.py](scripts/strata.py). Both scans write each key as a column
at the end of every call, and count the opportunities per value of it, so a
filter on a key has a matching denominator. Adding a key there adds it to both
tables.

| Key | Meaning |
|---|---|
| `tissue` | the read's tissue label, or `NA` |
| `kmerQ` | lowest base quality among the bases tested, see [kmerQ](#kmerq); only 93 remains after final filtering |

Every read position is counted once for each size, so the opportunities for one
size are the rows whose `context` has that length. `canonical_indel` in the
calls uses the same form, so the rate of a context is its calls divided by its
opportunities.

Workflow 4's table has `canonical_motif` and `k` (the tract's) in place of
`context`, and `count` is read bases inside tracts of that motif.

### Other tables

- `{model}_summary.tsv`: `canonical_motif`, `k`, `copies`, `trf_purity`,
  `model_purity`, `model`, `species`, `sample` and `count`, the number of tracts
  of each kind.
- `filter_report.tsv`: `filter`, `rows_before`, `rows_after`, `rows_removed`,
  see [the report](#final-filtering-and-the-filter-report).
- `indels_classified.tsv`: `indels.tsv` plus the five columns under
  [Germline or somatic](#tissue-and-cell-type-labels).

---

## Analysis notebooks

The R notebooks in [plotting/](plotting/) read `results/` and turn the tables into
rates with exact 95% Poisson intervals. They need R with readr, dplyr, tidyr and
ggplot2, and the IRkernel Jupyter kernel.

| Notebook | Reads | Shows |
|---|---|---|
| [tract_lengths.ipynb](plotting/tract_lengths.ipynb) | the `results/repetitive_dna/` summaries | how the M0 and M1 annotations (optionally M2) differ: tract counts, lengths, motifs and purity |
| [repetitive_mutation_rate.ipynb](plotting/repetitive_mutation_rate.ipynb) | `results/repetitive_indels/` and the summaries | rates inside repeats per period, motif, tract length and tissue, compared between models |
| [nonrepetitive_indels.ipynb](plotting/nonrepetitive_indels.ipynb) | `results/nonrepetitive_indels/` and `results/repetitive_indels/` | rates outside repeats per size, context, contig and tissue, and against one-unit InDels in repeat tracts |

Committed notebooks carry no outputs, while the copies on disk keep theirs.
[.gitattributes](.gitattributes) sends every notebook through
[scripts/strip_notebook_outputs.py](scripts/strip_notebook_outputs.py) on its
way into git. The filter is defined per clone, so set it up once after cloning:

```bash
git config filter.strip-notebook-outputs.clean "python3 scripts/strip_notebook_outputs.py"
git config filter.strip-notebook-outputs.required true
```

Without the first line, git commits notebooks with their outputs. The second
makes `git add` stop with an error if the filter fails, for instance on a
notebook that isn't valid JSON, rather than staging it unfiltered.

`results/` is small enough to copy from the cluster with rsync. From the
repository root of the local checkout, with `CLUSTER` the cluster login and
`REPO` the checkout there:

```bash
rsync -av CLUSTER:REPO/results/ results/
```

This shape has not been checked against the command used in practice.

Two approximations to keep in mind when reading the rates:

- **Rates by tract length or purity.** Calls carry their tract's `copies` and
  purity, but the opportunity table has only motif and k. The notebook splits
  each period's opportunities over tract-length bins in proportion to the repeat
  sequence in each bin, taken from the tract summaries, so rates by tract length
  are approximate.
- **Denominators outside and inside repeats differ.** Outside repeats, a read
  position counts for a motif when the k-mer starting there is that motif;
  inside, every read base of the tract counts. Both are per read base, so large
  differences between the two are meaningful and small ones are not.

---

## Tests

```bash
.venv/bin/python -m pytest
```

| File | Covers |
|---|---|
| [tests/test_pure_functions.py](tests/test_pure_functions.py) | canonicalisation, the homopolymer rule, reading the `sa` tag, Hamming purity, walking the haploid CIGAR |
| [tests/test_integration.py](tests/test_integration.py) | the callers on a small synthetic assembly and BAM with planted InDels: the right reads, positions, types and sizes; masking; windows adding up to the whole contig; one scan for three models matching three separate scans; the haploid check on several call tables |
| [tests/test_read_shapes.py](tests/test_read_shapes.py) | the scans on the read shapes real data has: reverse strand (with subread tags in sequencing order), base qualities below 93, reads ending early or crossing a window edge, soft clips, mapQ below 60, overlapping tracts |
| [tests/test_config.py](tests/test_config.py) | config.yaml: every sample belongs to a known species, override BAMs are well formed, the layout templates resolve, and every key is read by a workflow or a notebook |
| [tests/test_workflow_snapshots.py](tests/test_workflow_snapshots.py) | every target each workflow would run (name, inputs, outputs, options, command) for a fixture config, against the snapshots in `tests/snapshots/` |
| [tests/test_workflow_helper.py](tests/test_workflow_helper.py) | `workflows/pipeline.py`: paths from the layout templates, alignment overrides, the long contigs |
| [tests/test_strip_notebook_outputs.py](tests/test_strip_notebook_outputs.py) | the notebook clean filter, alone and through git |

The fixtures come from [tests/make_fixtures.py](tests/make_fixtures.py), which
`conftest.py` runs into a temporary directory at the start of every test
session; `python tests/make_fixtures.py` writes the same files to
`tests/fixtures/` to look at. The planted InDels and the expected outcome for
each are `PLANTED` in that file, and its `make_read()` builds the other read
shapes. The tests need no cluster and no external binaries.

The workflow snapshots show whether a change alters what the workflows run.
After a deliberate change, rewrite them and review their diff:

```bash
UPDATE_SNAPSHOTS=1 .venv/bin/python -m pytest tests/test_workflow_snapshots.py
git diff tests/snapshots/
```

### Replaying windows of real data

The mirror is for software testing only: it shows whether the code does on a
laptop what it does on the cluster, and no data conclusions are drawn from it.
[scripts/export_local_subset.py](scripts/export_local_subset.py) copies windows
of samples from the cluster, with everything workflows 3 and 4 need, and the
cluster's tables for each window. On the cluster, from the repository root:

```bash
python scripts/export_local_subset.py human/chk:h1tg000001l:0-10000000 \
    chimp/yoran human/ob006 chimp/carl --window-size 1000000
```

A sample without a window gets the `--window-size` window with the most
unique calls. Copy `mirror/` to the local checkout with rsync, then replay
every window:

```bash
for window in mirror/*/*/*/; do PYTHON=.venv/bin/python bash "$window/run_local.sh"; done
```

Each `run_local.sh` reruns both workflows' scans, haploid checks and final
filters on its window, and compares every table with the cluster's using
[scripts/compare_tables.py](scripts/compare_tables.py): by content, in any row
order, with numbers compared by value and the species and sample columns the
workflows append ignored. A window that is one whole scan job (10 Mb) has its
opportunities compared too. The script exits 1 if any table differs.

---

## Repository layout

```
config.yaml                        samples, input paths, tools
environment.yml                    conda environment for the cluster
requirements.txt                   Python packages for local work and tests
workflows/
  pipeline.py                      config loading and paths from the layout templates, shared by the workflows
  precompute.py                    workflow 1
  repetitive_dna.py                workflow 2
  indel_mutation_call.py           workflow 3
  repetitive_indel_call.py         workflow 4
scripts/
  utilities.py                     the pileup scans: calls, opportunities, unresolved haplotypes
  nonrepetitive_scan.py            workflow 3's scan job
  tandem_repeat_scan.py            workflow 4's scan job
  resolve_haplotypes.py            unresolved-haplotype regions of one window
  parse_mreps_canonical.py         mreps report -> repeat BED (M0, M2)
  parse_trf.py                     TRF output -> repeat BED (M1); canonicalise(), Hamming purity
  resolve_repeat_overlaps.py       one tract per base per model
  filtering_haploid_mappping.py    haploid concordance check
  filtering_small_indels.py        final filters on the calls, with the report
  filtering_opportunities.py       final filters on the opportunities
  join_read_classification.sh      germline/somatic columns
  plot_indels.R                    diagnostic PDF
  export_local_subset.py           copy windows of samples from the cluster, to replay locally
  compare_tables.py                compare a replayed table with the cluster's, by content
  strip_notebook_outputs.py        git clean filter: notebooks are committed without outputs
plotting/                          analysis notebooks (R)
tests/                             pytest suite and fixture generator
```

Not tracked by git: `results/`, `mirror/`, `plots/`, `figs/` and `papers/`
(local data, figures and papers), `TRF/` and `mreps/` (tool sources, see
*Requirements*), and `tests/fixtures/`, which `make_fixtures.py` writes.

The planned changes, and the reason for each, are in
[documents/decisions_2026-10-19.md](documents/decisions_2026-10-19.md). Where it
and this README disagree, the decisions document is the intent.

## License

GNU General Public License v2, see [LICENSE.md](LICENSE.md).
