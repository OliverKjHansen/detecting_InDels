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
   Each call records whether the same InDel shows up there.
5. **Filter.** Thresholds on base quality and on position in the read remove
   what is most likely sequencing error.

Tandem repeats mutate by replication slippage, at rates far above the rest of the
genome, so they are measured separately. The pipeline annotates repeats first,
masks them out for one caller, and scans only them with a second caller that has
its own size rule. Every base falls in one of three groups:

| Group | What's in it | Measured by |
|---|---|---|
| Repeat | tracts with periods 2–6 | workflow 4 |
| Excluded | homopolymer runs of 6 bp or more, and tracts with periods of 7 or more | not workflow 3, which masks them; workflow 4 drops the tracts that overlap a homopolymer run or a sub-window TRF didn't finish, but measures a tract that overlaps a long-period tract, those long-period bases included |
| Non-repetitive | everything else | workflow 3 |

Homopolymers are where HiFi makes the most errors, and repeats with longer
periods still slip faster than the background, so leaving them in the
non-repetitive group would raise its rate. Repeat bases workflow 4 can't
measure, the tracts that overlap a homopolymer run and the parts of tracts
that [overlap resolution](#workflow-2-repeat-annotation) drops, stay masked
from workflow 3, so they are in neither workflow either. So are the
[unresolved-haplotype regions](#workflow-1-precompute), whatever group their
bases belong to.

### The four workflows

| # | Workflow | What it does | Needs | Writes |
|---|---|---|---|---|
| 1 | [precompute.py](workflows/precompute.py) | Per-contig FASTAs, mean coverage, unresolved-haplotype regions, read → tissue table | the inputs | `output/precompute/`, `results/precompute/` |
| 2 | [repetitive_dna.py](workflows/repetitive_dna.py) | Tandem repeat tracts with periods 2–6, under three models (M0, M1, M2); the excluded group | 1 | `output/repetitive_dna/`, `results/repetitive_dna/` |
| 3 | [indel_mutation_call.py](workflows/indel_mutation_call.py) | 1–3 bp InDels outside every repeat tract, excluded region and unresolved region | 1, 2 | `results/nonrepetitive_indels/` |
| 4 | [repetitive_indel_call.py](workflows/repetitive_indel_call.py) | InDels inside the tracts of each model that overlap no homopolymer run, unfinished TRF sub-window or unresolved region, classed as slippage or other | 1, 2 | `results/repetitive_indels/` |

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
- **Crash-safe outputs.** Scan and repeat-annotation jobs, and the joins of
  windows and contigs, write to a `.tmp` file and rename it when done, so a job
  that dies never leaves a finished-looking file for GWF to skip.

---

## Getting started

### Requirements

| Kind | What |
|---|---|
| Workflow engine | [GWF](https://gwf.app/) ≥ 1.7, submitting to Slurm |
| Python | pysam, pandas, intervaltree, pyyaml, duckdb; pytest for the tests |
| Binaries | samtools, bedtools, pbmm2, faSplit, TRF 4.10, mreps |
| R | ggplot2, dplyr, tidyr for the diagnostic plots; also readr and the IRkernel Jupyter kernel for the notebooks, and duckdb to query the opportunity cube |

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
| `precompute`, `repeats`, `scan`, `strata`, `final_filters`, `classification` | The workflows' parameters, below. The scripts take them as arguments; the scans take the whole `strata` block as `--key-parameters`. |

No workflow reads `sperm_samples` or `sample_age`; the notebooks do.

Parameters:

| Parameter | Value | Key |
|---|---|---|
| Unresolved-haplotype scan window | 40 Mb | `precompute.unresolved_window_size` |
| Unresolved-haplotype evidence | ≥ 2 reads with the same base at quality ≥ 40 | `precompute.unresolved_min_reads`, `precompute.unresolved_min_base_quality` |
| How hap1's and hap2's contig names start | `h1tg`, `h2tg` | `precompute.haplotype_prefixes` |
| Repeat models | M0, M1, M2 | `repeats.models` |
| Repeat periods | 2, 3, 4, 5, 6 | `repeats.periods` |
| mreps minimum exponent (M0, M2) | 3 copies | `repeats.mreps_min_exponent` |
| mreps resolution (M2) | 1 | `repeats.m2_resolution` |
| TRF scoring (M1, and the long periods) | `2 7 7 80 10 50`, maximum period 6 | `repeats.trf_params` |
| TRF window and overlap (M1; a long-period job covers one such window) | 5 Mb, plus 1 Mb on each side | `repeats.trf_window_size`, `repeats.trf_overlap` |
| Long-period sub-windows, their overlap, and TRF's time limit on each | 500 kb, plus 10 kb on each side; 30 min | `repeats.long_period_sub_window_size`, `repeats.long_period_overlap`, `repeats.long_period_timeout` |
| Excluded homopolymer runs | 6 bp or more | `repeats.homopolymer_min_length` |
| Excluded periods | 7 (one above `repeats.periods`) to 2000 | `repeats.long_period_max` |
| InDel scan window, and workflow 3's windows per job (scanned at the same time, one core and the window's memory each) | 10 Mb; 4 | `scan.window_size`, `scan.windows_per_job` |
| Workflow 3: the contigs merged in one job, by their total length | 500 Mb (chk: 13 jobs of each merge) | `scan.merge_bp_per_job` |
| Coverage window | 0.5×–2× the sample mean | `scan.coverage_band` |
| InDel sizes outside repeats (workflow 3) | 1–3 bp | `scan.indel_sizes` |
| Minimum kmerQ, as the key `kmerQ_pass` | 93 | `strata.min_kmerq` |
| Workflow 3: the bins of `inner_flank_q` and `outer_flank_q`; the inserted bases' minimum quality | 60, 93; 93 | `strata.flank_quality_bins`, `strata.min_inserted_quality` |
| Workflow 3's cube (ticket 49): the keys it leaves out, and its own bins of the nearest other event, mapping quality and passes | `strand_passes_bin`; 10, 20, 50, 100, 1000 bp; 1, 20, 60; 7, 12, 17 | `strata.nonrepetitive_without`, `strata.nonrepetitive_bins` |
| Bins of read length and of distance to the read end, keys on calls and opportunities | 4, 8, 12, 16, 20 kb; 25, 100, 300, 500, 1000 bp | `strata.read_length_bins`, `strata.read_end_distance_bins` |
| Bins of mapping quality, passes, predicted accuracy, passes on the weaker strand, and distance to the nearest other event in the read | 1, 10, 20, 30, 40, 50, 60; 5, 7, 10, 12, 15, 18, 20, 25, 30; 0.99, 0.999, 0.9999; 2, 5, 10; 5, 10, 20, 25, 35, 50, 100, 300, 500, 1000 bp | `strata.mapq_bins`, `.passes_bins`, `.accuracy_bins`, `.strand_passes_bins`, `.nearby_event_bins` |
| Assembly bases on each side of a call's anchor | 25 | `scan.flank_width` |
| The definition of a final call: the values each key keeps, on calls and opportunities alike; then workflow 3's steps on the calls only | workflow 3: `inner_flank_q` ≥ 93, `outer_flank_q` ≥ 60, `base_after_masked` = `no`; then the inserted bases, the haploid check and the hap1/hap2 test. Workflow 4: `kmerQ_pass` = `pass`, `read_end_distance_bin` ≥ 100 bp | `final_filters.keep`, `.nonrepetitive_keep`, `.nonrepetitive_calls_only`, `.repeat_keep` |
| SMRT Cells taken to be capped, and dropped, if none of their reads' base qualities is above; reads looked at per sample | Q40; 100,000 | `final_filters.capped_quality`, `precompute.smrt_cell_reads` |
| Sequel II samples (the others are Revio) | da1, ob2v, stephan, carl, mbewe, samson | `sequel_ii_samples` |
| Bases a read must cover on each side of a tract to count for it (workflow 4) | 25 bp | `scan.min_tract_flank` |
| Bins of a tract's copies and purity, keys of workflow 4 | 5, 10, 20, 50; 90, 95, 99, 100% | `strata.tract_copies_bins`, `strata.tract_purity_bins` |
| Bins of subread support, a key on calls and opportunities | 0.8, 0.9, 1 | `strata.subread_support_bins` |
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
| Per-read cell-type likelihoods | `{classification}/{species}/{sample}/{sample}/05_ReadsClassification/{sample}_All_Reads_Likelihood_merged.tsv` | 1, 3, 4, `methylation_samples` only |

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
workflows share. Outputs go to `output/precompute/{species}/{sample}/` (spelled
`precomute` before the rerun; mirror windows exported then keep that name), and
the one summary to `results/precompute/{species}/{sample}/`.

| Step | Target | Runs per | Output |
|---|---|---|---|
| 1. Split the assembly | `faSplit_{species}_{sample}_diploid` | sample | `contigs/{contig}.fa` |
| 1b. Haplotypes alone | `haplotypes_{species}_{sample}_diploid` | sample | `haplotypes/hap1.fa`, `haplotypes/hap2.fa` |
| 2. Mean coverage | `average_coverage_{species}_{sample}` | sample | `long_contigs.bed`, `avg_coverage.txt` |
| 2b. SMRT Cells | `smrt_cells_{species}_{sample}_diploid` | sample | `smrt_cells.tsv`: each SMRT Cell's highest base quality, by [smrt_cell_qualities.py](scripts/smrt_cell_qualities.py), for final filtering's Q40 rule |
| 3a. Unresolved haplotypes | `unresolved_haplotype_{species}_{sample}_{contig}_{start}_diploid` | 40 Mb window | `unresolved_haplotypes/chunks/{contig}_{start}_{end}.bed` |
| 3b. Join the windows | `merge_unresolved_{species}_{sample}_{contig}_diploid` | contig | `unresolved_haplotypes/{contig}.bed` |
| 3c. Share masked | `unresolved_summary_{species}_{sample}_diploid` | sample | `results/precompute/{species}/{sample}/unresolved_summary.tsv` |
| 4. Tissue table | `sorting_tissue_classification_{species}_{sample}_diploid` | sample in `tissue_samples` | `sorted_tissue_classification.tsv` |
| 5. Read classification | `read_classification_{species}_{sample}_diploid` | sample in `methylation_samples` | `read_classification.tsv`: every read's germline/somatic score and bin, by [read_classification.sh](scripts/read_classification.sh) |
| 5b. Read bins, compact | `read_bins_{species}_{sample}_diploid` | sample in `methylation_samples` | `read_bins.npz`: the same bins, as the scans read them, by [read_bins.py](scripts/read_bins.py) |

**Split.** `faSplit byname` writes one FASTA per contig, for mreps in workflow 2.

**Haplotypes alone.** The diploid assembly holds both haplotypes' contigs,
told apart by how their names start (`precompute.haplotype_prefixes`: hifiasm's
`h1tg` and `h2tg`). Each haplotype's contigs are copied to a FASTA of their
own, for the [hap1/hap2 test](#the-hap1hap2-test).

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
3. A region that overlaps or touches the one being built joins it, which
   becomes their union: a later region never cuts it short.

A region can reach past its window, so the windows' BEDs are joined with
`bedtools merge`. The result is a headerless BED (`contig, start, end`) per
contig. Both InDel workflows leave the regions out: they are part of workflow
3's mask, and workflow 4 drops every tract that overlaps one. In such a region
a read can align to the wrong haplotype, and in a tandem repeat that turns a
germline length difference into a false call. `unresolved_summary.tsv` gives
the share of the genome they mask, per sample (see [Other tables](#other-tables)).

**Tissue table.** The read label file (read ID and tissue, whitespace-separated)
becomes a sorted two-column TSV with a `readID  tissue` header.

**Read bins.** Every scan job of workflows 3 and 4 looks up the germline/somatic
bins of its window's reads. Read whole in each job, `read_classification.tsv`
(about 600 MB for a 30x sample) made the jobs wait on the storage, as hundreds
of them read the same file at once (ticket 42). So step 5b writes the bins
once more, compact: per read a 64-bit key, a hash of its name, sorted, and the
bin as one byte, about 100 MB for 10 M reads, loaded in under a second
(ticket 44). A read missing from it gets `NA`, as from the table.

---

## Workflow 2: repeat annotation

[workflows/repetitive_dna.py](workflows/repetitive_dna.py) annotates tandem
repeats with periods k = 2 to 6 under three models, and the excluded group:
homopolymer runs of 6 bp or more and tracts with periods of 7 or more. Outputs
go to `output/repetitive_dna/{species}/{sample}/`.

The tracts serve both InDel workflows. Workflow 3 masks out the union of all
three models and the excluded group; workflow 4 scans each model's tracts on
their own, without those that overlap a homopolymer run or a sub-window TRF
didn't finish.

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
| 1b. Resolve overlaps | `{model}_merge_k_{contig}_{species}_{sample}_diploid` | contig × model | `{contig}_{model}.bed`, `bases/{contig}_{model}.bed` |
| 2. Sample BED and summaries | `cat_{model}_{species}_{sample}_diploid` | sample × model | `{model}.bed`, `results/repetitive_dna/{species}/{sample}/{model}_summary.tsv` and `{model}_dropped.tsv` |
| 3a. Long-period tracts | `long_period_trf_{contig}_{start}_{species}_{sample}_diploid` | 5 Mb window, TRF on its 500 kb sub-windows in turn | `excluded/chunks/{contig}_long_periods_{start}.bed` |
| 3b. Excluded group | `excluded_{contig}_{species}_{sample}_diploid` | contig | `excluded/{contig}_homopolymers.bed`, `excluded/{contig}_long_periods.bed`, `excluded/{contig}.bed` |
| 3c. Excluded summary | `excluded_summary_{species}_{sample}_diploid` | sample | `results/repetitive_dna/{species}/{sample}/excluded_summary.tsv` |

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
and 6. A dropped tract is not trimmed, so its bases outside the kept tract are
in no tract of the model, and workflow 4 does not see them. Workflow 3 masks
them all the same, because its mask is built from every tract before
resolution: the same step writes the bases in any of the model's tracts, as
merged intervals, to `bases/{contig}_{model}.bed`. So those bases are in
neither workflow, like the [excluded group](#workflow-2-repeat-annotation), and
`{model}_dropped.tsv` counts them.

**Sample BED and summaries.** The contigs' resolved BEDs are concatenated into
`{model}.bed`. The summary counts tracts per distinct combination of
`canonical_motif`, `k`, `copies`, `trf_purity` and `model_purity`, which keeps a
whole genome small enough to plot locally. `{model}_dropped.tsv` gives the bases
overlap resolution leaves out of the model (see [Other tables](#other-tables)).

**The excluded group.** Two kinds of sequence are in neither InDel workflow,
except that workflow 4 measures a period 2–6 tract that overlaps a long-period
tract, its long-period bases included (ticket 34):

- **Homopolymer runs** of at least `homopolymer_min_length` (6) bp, the shortest
  period-2 tract M0 annotates. [find_homopolymers.py](scripts/find_homopolymers.py)
  reads the contig's FASTA, ignoring case; runs of `N` are not homopolymers.
  Shorter runs are left to workflow 3's [homopolymer rule](#homopolymer-rule).
- **Tracts with periods from 7 to `long_period_max` (2000)**, from a TRF run of
  their own with M1's scoring (minimum score 50) and maximum period 2000. M1
  keeps its own run with maximum period 6, so M1 does not change.
  [parse_trf.py](scripts/parse_trf.py) keeps the periods 7–2000, but skips a
  tract whose consensus motif is whole copies of a motif shorter than 7: TRF can
  report a short repeat again at a multiple of its period, `(AC)n` as
  `ACACACAC`, and that is M1's tract. These tracts are only masked, so a window
  keeps every tract it reports rather than only those starting in its core: an
  array longer than the overlap is then covered by the windows' pieces.
- **Satellite arrays TRF can't finish.** Up to period 2000, TRF slows down with
  the length of a satellite array in its window: on the cluster some 5 Mb
  windows ran for more than 8 hours. So each job, one per M1 window, runs TRF
  on 500 kb sub-windows (10 kb more on each side) in turn, each within 30
  minutes (`repeats.long_period_timeout`; away from satellites TRF takes
  seconds). A sub-window it doesn't finish is excluded whole, which errs on the
  side of masking. Its window's BED in `excluded/chunks/` gives it with
  `timed_out` in the motif column, and `excluded_summary.tsv` sums them.

Per contig, `excluded/{contig}_homopolymers.bed` lists the runs (contig, start,
end, base), `excluded/{contig}_long_periods.bed` the long-period tracts as
merged intervals, and `excluded/{contig}.bed` their union, which workflow 3
masks. Workflow 4 drops the tracts that overlap a homopolymer run or a
sub-window TRF didn't finish, but measures a tract that overlaps a long-period
tract (the shorter period wins, ticket 34). The merging is done by
[merge_regions.py](scripts/merge_regions.py), which joins intervals that overlap
or touch, as `bedtools merge` does. `excluded_summary.tsv` gives the bases of
each kind and of their union per sample (see [Other tables](#other-tables)).

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
1–3 bp InDels at positions outside every repeat tract, every excluded region
and every unresolved-haplotype region. Intermediates go to
`output/nonrepetitive_indels/{species}/{sample}/`.

| Step | Target | Runs per | Output |
|---|---|---|---|
| 1. Build the mask | `non_rep_Mask_{species}_{sample}_diploid` | sample | `mask/{contig}.bed` |
| 2. Scan | `non_rep_Scan_{contig}_{start}_{species}_{sample}_diploid`, named after its first window | `scan.windows_per_job` 10 Mb windows, scanned at the same time | `calls/chunks/`, `opportunities/chunks/`, `fastq/chunks/` |
| 3. Join the windows | `non_rep_MergeWindows_{contig}_{species}_{sample}_diploid`, named after its first contig | contigs adding up to `scan.merge_bp_per_job`, one after the other | `calls/{contig}.tsv`, `fastq/{contig}.fastq.gz` |
| 3b. Sum the opportunities | `non_rep_MergeOpportunities_{contig}_{species}_{sample}_diploid`, named after its first contig | the same contigs; a sum that fails its check is retried up to three times in the job (ticket 41) | `opportunities/{contig}.parquet` |
| 4. Remap to haploid | `non_rep_HaploidMapping_{species}_{sample}_diploid` | sample | `haploid/all.fastq.gz`, `haploid/remapped.bam`, `haploid/locations.tsv` |
| 4b. Remap to each haplotype | `non_rep_HaplotypeMapping_{haplotype}_{species}_{sample}_diploid` | sample × haplotype | `haploid/{haplotype}.bam`, `haploid/{haplotype}_locations.tsv` |
| 5. Haploid check | `non_rep_HaploidFilter_{species}_{sample}_diploid` | sample | `merged/indels.tsv` |
| 6. Merge opportunities, plot | `non_rep_MergeAndPlot_{species}_{sample}_diploid` | sample | `merged/opportunities.parquet`, `plots/nonrepetitive_indels/{species}/{sample}/plots.pdf` |
| 7. Final filtering | `non_rep_FinalFiltering_{species}_{sample}_diploid` | sample | `indels.tsv`, `opportunities.parquet`, `filter_report.tsv` in `results/nonrepetitive_indels/{species}/{sample}/` |
| 8. Germline/somatic | `non_rep_Classify_{species}_{sample}_diploid` | sample in `methylation_samples` | `indels_classified.tsv`, next to `indels.tsv` |

1. **Mask.** Per contig, every tract of the three models before overlap
   resolution (`bases/{contig}_{model}.bed`), the unresolved-haplotype regions
   and the excluded group, merged with [merge_regions.py](scripts/merge_regions.py).
2. **Scan.** [nonrepetitive_scan.py](scripts/nonrepetitive_scan.py) runs
   `nonrepetitive_scan()` from [utilities.py](scripts/utilities.py). At each
   position of the window, every read with a base there (each event's anchor
   tested against the coverage window and its bases against the mask, see
   [Coverage](#coverage) and [The anchor base](#the-anchor-base)):
   - adds its opportunities: a deletion of each size (1, 2 and 3 bp), under
     the canonical k-mer of that size starting at the read position, and an
     insertion of each size in the gap before the read position, under the
     canonical pair of bases around it, each with the [keys](#keys) of the
     read and of the site (an event without two read bases on each side is
     not counted);
   - becomes a call if it carries a 1–3 bp InDel anchored there that passes the
     homopolymer rule and touches no masked base (see [The anchor
     base](#the-anchor-base)), and then goes to the FASTQ.
3. **Join.** Per contig, the windows' calls are concatenated, their opportunity
   counts summed per value of every other column into the contig's part of the
   [opportunity cube](#opportunities), and their FASTQs concatenated keeping
   each read once. The counts are summed
   ([opportunity_tables.py](scripts/opportunity_tables.py)) in a job of their
   own (3b): the calls and reads, and the remapping after them, don't depend on
   how the opportunities are kept. A contig can have hundreds of millions of
   keys (374 M in one of ob007's), more than a sum can hold in memory, so
   `sort` puts the windows' lines in key order, with the job's 4 cores and half
   its memory, spilling to its `$TMPDIR` on the node's local disk, and the
   lines of each key, then next to each other, are added up one key at a time
   (ticket 50). DuckDB only writes the summed lines to Parquet. The contig's
   cube is kept only if it holds the windows' total, counted as they are added
   up (and `sort` must give back as many lines as it was given), with no row
   without a count and no key with a NUL byte (ticket 41):
   DuckDB 2.0 has misread tables at random, into such rows or by losing lines,
   which only the total shows. A sum that fails the check leaves no cube, so the
   next `gwf -f workflows/<name>.py run` submits it again; one that fails
   again with the same numbers points at the windows themselves. Workflow 4's
   sum (step 2b) is the same. Each contig is its own job, so the contigs merge
   side by side and a failure costs one contig.
4. **Remap.** All reads carrying a call are aligned to the haploid assembly with
   `pbmm2 align --preset HIFI --sort`. The alignments are reduced to read ID,
   flag, contig, position, mapping quality and CIGAR, sorted by read ID. The
   same reads are aligned the same way to hap1 alone and to hap2 alone (4b),
   for the [hap1/hap2 test](#the-hap1hap2-test).
5. **Haploid check.** [filtering_haploid_mapping.py](scripts/filtering_haploid_mapping.py)
   checks every call against its read's haploid alignment and its sequence in
   `haploid/all.fastq.gz` (see [Haploid concordance](#haploid-concordance)) and
   writes one table for the sample, with `species` and `sample` columns added.
6. **Opportunities and plots.** The contigs' opportunity cubes are joined into
   the sample's, with `species` and `sample` added.
   [plot_indels.R](scripts/plot_indels.R) draws diagnostics of the calls
   before final filtering (it doesn't read the opportunities): kmerQ, the 20
   most common InDel sequences, distance to the read end, the unique flag, the haploid
   check, tissue, kmerQ against haploid distance, and subread coverage.
7. **Final filtering.** [final_filtering.py](scripts/final_filtering.py)
   applies the definition of a final call to the calls and the opportunities
   alike (see [Final filtering and the filter
   report](#final-filtering-and-the-filter-report)). The haploid check stays a
   column of the calls.
8. **Germline/somatic.** See [Tissue and cell-type labels](#tissue-and-cell-type-labels).

---

## Workflow 4: InDels inside repeats

[workflows/repetitive_indel_call.py](workflows/repetitive_indel_call.py) calls, for
each repeat model, every InDel in that model's tracts, of any size: anchored
inside a tract or just before one (see [The anchor base](#the-anchor-base)).
Each is classed as slippage or other (see [InDel size](#indel-size)).

Intermediates go to `output/repetitive_indels/{species}/{sample}/`.

| Step | Target | Runs per | Output |
|---|---|---|---|
| 0. Tracts measured | `TRTracts_{species}_{sample}` | sample | `{model}/tracts/{contig}.bed`, `results/repetitive_indels/{species}/{sample}/tracts_summary.tsv` |
| 1. Scan | `TRScan_{contig}_{start}_{species}_{sample}` | 10 Mb window | `{model}/calls/chunks/`, `{model}/opportunities/chunks/`, `fastq/chunks/` |
| 2. Join the windows | `TRMergeWindows_{contig}_{species}_{sample}` | contig | `{model}/calls/{contig}.tsv`, `fastq/{contig}.fastq.gz` |
| 2b. Sum the opportunities | `TRMergeOpportunities_{contig}_{species}_{sample}` | contig | `{model}/opportunities/{contig}.parquet` |
| 3. Remap to haploid | `TRHaploidMapping_{species}_{sample}` | sample | `haploid/all.fastq.gz`, `haploid/remapped.bam`, `haploid/locations.tsv` |
| 3b. Remap to each haplotype | `TRHaplotypeMapping_{haplotype}_{species}_{sample}` | sample × haplotype | `haploid/{haplotype}.bam`, `haploid/{haplotype}_locations.tsv` |
| 4. Haploid check | `TRHaploidFilter_{species}_{sample}_{model}` | sample × model | `{model}/merged/indels.tsv` |
| 5. Merge opportunities, plot | `TRMergeAndPlot_{species}_{sample}_{model}` | sample × model | `{model}/merged/opportunities.parquet`, `plots/repetitive_indels/{species}/{sample}/{model}/plots.pdf` |
| 6. Final filtering | `TRFinalFiltering_{species}_{sample}_{model}` | sample × model | `indels.tsv`, `opportunities.parquet`, `filter_report.tsv` in `results/repetitive_indels/{species}/{sample}/{model}/` |
| 7. Germline/somatic | `TRClassify_{species}_{sample}_{model}` | sample × model, in `methylation_samples` | `indels_classified.tsv`, next to `indels.tsv` |

**Tracts measured.** [drop_excluded_tracts.py](scripts/drop_excluded_tracts.py)
drops every tract of a model that shares a base with a homopolymer run
(`excluded/{contig}_homopolymers.bed` from workflow 2), a sub-window the
long-period TRF didn't finish (marked `timed_out` in the window BEDs in
`excluded/chunks/`) or an unresolved-haplotype region (from workflow 1). A tract
that overlaps a long-period tract is measured: the shorter period wins, as in
overlap resolution (ticket 34; the long-period bases stay in workflow 3's mask).
A tract is dropped whole, so it is measured whole or not at all, for calls and
opportunities alike, and its bases stay in workflow 3's mask, so they are in
neither workflow. A tract that only touches such a region is kept.
`tracts_summary.tsv` gives each model's tracts and their bases before and after
(see [Other tables](#other-tables)), and
[tract_drop_causes.py](scripts/tract_drop_causes.py) counts, per model, the
tracts each cause drops and what the rules ticket 34 weighed would measure.

Steps 1–7 mirror workflow 3, with two savings. One scan of a window serves all
three models ([tandem_repeat_scan.py](scripts/tandem_repeat_scan.py) runs
`tandem_repeat_scan()`), and the reads carrying a call in any model are remapped
to the haploid assembly once per sample. Each model's output is the same as
scanning that model alone, which [tests/test_integration.py](tests/test_integration.py)
checks.

The unit is one tract in one read that spans it (Q22): the read's alignment
covers the tract and at least `scan.min_tract_flank` (25) bases on each side.
A read ending inside the tract or its flank can't show a change in the
tract's length, so it counts for neither calls nor opportunities. For each
measured tract of a model that begins in the window, and whose first base
passes the coverage window, every read that spans it:

- adds one opportunity under the tract's motif and k, with the tract keys
  (its copies and purity, the read's flank on each side) and the read and
  site keys;
- gives a call for each of its InDels counted in the tract ([The anchor
  base](#the-anchor-base)), with its class and size and the same keys as its
  opportunity. The call also carries the tract's
  `canonical_motif`, `k`, `copies`, `trf_purity`, `model_purity` and `model`.

A tract's call and opportunity on a read are tested on the same bases (R7 in
[the workflow 4 review](documents/workflow4_review.md)): kmerQ and subread
support over the read's bases from two before the tract to two after it,
whatever the InDel's size and place. Its distance to the read end is the
nearer flank, and its nearest event the read's nearest InDel or mismatch other
than the tract's own InDels. Workflow 2 leaves one tract per base per model;
the scan stops with an error if a model's tracts overlap (R8).

How the two callers differ:

| | Workflow 3, outside repeats | Workflow 4, inside repeats |
|---|---|---|
| Events counted | touching no masked base, from the anchor to k bases on | touching a tract of the model, from the anchor to 6 bases on, on a read spanning the tract |
| Excluded group | masked | tracts overlapping a homopolymer run or an unfinished TRF sub-window dropped; not those overlapping a long-period tract |
| Unresolved-haplotype regions | masked | tracts overlapping one dropped |
| InDel size | 1–3 bp | any, classed slippage or other |
| Homopolymer rule | applied | not applied |
| `canonical_indel` | the smaller of the sequence and its reverse complement | the smallest rotation of the sequence or its reverse complement, like the tract motif |
| Opportunity unit | read position × deletion size, by canonical k-mer; gap between two read bases × insertion size, by the canonical pair around it | tract × read spanning it with 25 bp on each side, by tract motif, k, copies, purity and flanks |
| Final filters, on calls and opportunities alike | capped SMRT Cells, kmerQ ≥ 93 | capped SMRT Cells, kmerQ ≥ 93, tract ≥ 100 bp from the read's ends |
| Runs | once per sample | once per model per sample |

---

## Filtering in detail

Filters act at four points:

- **During the scan**, before anything is written. These decide which positions
  and reads count at all, and they apply to calls and opportunities alike.
- **At the haploid check**, which annotates every call and drops none; where it
  cannot check a call, it leaves `NA`.
- **At final filtering**, which applies thresholds to recorded columns and
  reports what each one removed.
- **In the analysis notebooks**, which apply a few more to the final tables.

### Where each filter acts

| Filter | Stage | Workflow 3 | Workflow 4 | Opportunities |
|---|---|---|---|---|
| Contig longer than 1 Mb | setup | ✓ | ✓ | ✓ |
| Excluded group (homopolymer runs ≥ 6 bp, periods ≥ 7) | setup | masked | tracts overlapping a homopolymer run or an unfinished TRF sub-window dropped | ✓ |
| Unresolved-haplotype regions | setup | masked | tracts overlapping one dropped | ✓ |
| Coverage within 0.5×–2× the sample mean | scan | ✓ | ✓ | ✓ |
| Position outside the mask / inside a tract | scan | outside | inside | ✓ |
| InDel size | scan | 1–3 bp | any (a class) | – |
| Read spans the tract and 25 bp on each side | scan | – | ✓ | ✓ in workflow 4 |
| Homopolymer rule | scan | ✓ | – | ✓ in workflow 3, deletions |
| Same InDel in the haploid alignment | not applied, a column (`haploid_check`) | – | – | – |
| kmerQ ≥ 93 | final | ✓ | ✓ | ✓ |
| SMRT Cell not capped at Q40 | final | ✓ | ✓ | ✓ |
| Tract at least 100 bp from the read's ends | final | – | ✓ | ✓ in workflow 4 |
| No other InDel at the position (`unique`) | notebooks | ✓ | ✓ | – |

Nothing is filtered on read flags (secondary, supplementary and duplicate
alignments all enter the pileup), on mapping quality, or on base quality outside
the kmerQ window. Mapping quality to the diploid assembly is kept in the
opportunity tables (`mapQ`) but not in the calls, which carry the haploid
mapping quality (`map_mapq`) instead.

### The anchor base

pysam reports an InDel on the reference base just before the inserted or deleted
sequence: the anchor. The coverage test is made at the anchor, and a call's
`position` column is the anchor's 1-based coordinate on the contig. An
opportunity's anchor is the base before its k-mer or gap.

pbmm2 places an InDel that can slide at its leftmost equivalent position, so
slippage in a pure tract is anchored on the base before the tract, or before a
partial unit in front of it (`…A|CACACA…`, deleting `AC`), and slippage of a
tract that follows another on the first tract's last base. So the mask and the
tracts are tested on every base an event touches, not the anchor alone (L1 in
[the workflow 3 review](documents/nonrepetitive_indel_review.md), R1 and R2 in
[the workflow 4 review](documents/workflow4_review.md)):

- **Workflow 3** counts an event of k bases anchored at base a, call or
  opportunity, deletion or insertion, only if bases a to a+k are all unmasked.
- **Workflow 4** counts a call in the tract of the first base after the anchor,
  up to k and at most the longest period (6) bases on, that lies in a tract of
  the model, or else in the anchor's tract.

So an event of 1–3 bp near a tract is counted by exactly one of the two
workflows, as long as its read spans the tract (workflow 4 counts only those).
A longer event that only reaches a tract more than 6 bases after its anchor is
counted by neither: workflow 3 calls only 1–3 bp.

### Coverage

A position is used only if the number of reads covering it (pysam's
`nsegments`, which includes reads with a deletion there) is between half and
twice the sample's mean coverage. This removes positions with too many reads,
such as collapsed duplications, and too few, such as poorly assembled or poorly
mapped ones. Workflow 3 makes the test at the anchor, for calls and
opportunities alike: a call's anchor is the base before its gap, an
opportunity's the base before its k-mer or gap (ticket 31; before, an
opportunity was tested one base later). Workflow 4 makes it once per tract, at
the tract's first base, for the tract's calls and opportunities alike.

### Mask and tracts

Workflow 3 skips every event that touches its mask (see [The anchor
base](#the-anchor-base)): the union of the M0, M1 and M2 tracts (every tract,
before overlap resolution), the excluded group and the unresolved-haplotype
regions. Workflow 4 uses only the measured tracts of the model being called: a
tract that overlaps no homopolymer run, no sub-window TRF didn't finish and no
unresolved-haplotype region (it may overlap a long-period tract). After overlap
resolution a base lies
in at most one tract per model, so each position has a single tract, motif and
k. (The code's tie-breaks, the
largest k for a call and the smallest for an opportunity, only matter for BEDs
that were not resolved.)

### InDel size

Workflow 3 keeps InDels of 1, 2 or 3 bp. Workflow 4 keeps every InDel counted
in a tract, with no upper limit, and gives it a class (Q23), `indel_class`:

- **slippage**: whole copies of the tract's motif, in any rotation, on either
  strand: `AC`, `CA`, `GT`, `TG` or `ACAC` in an `AC` tract;
- **other**: everything else, such as `C`, `ACA` or `GG` in an `AC` tract.

Each class gets its own rate over the same opportunities. Before, InDels whose
size isn't a multiple of k were counted nowhere, and multiples were kept
whatever their sequence (R5 in [the workflow 4 review](documents/workflow4_review.md)).

### Homopolymer rule

Homopolymer runs of 6 bp or more are in workflow 3's mask, with the rest of the
[excluded group](#workflow-2-repeat-annotation). Shorter runs are not masked, so
workflow 3 removes InDels in them itself. It
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
together with the base on each side of it. It counts an insertion opportunity
at every gap: whatever the bases around a gap, some inserted sequence passes
the rule there (a base differing from both). Workflow 4 has no homopolymer rule:
its periods start at 2, and a 1 bp InDel is never slippage.

### kmerQ

`kmerQ` is the lowest base quality of the two read bases on each side of the
event, never of the event's own bases:

```
deletion of CG      assembly  G A T C G T A C
                    read      G A T - - T A C
                    kmerQ       ^ ^     ^ ^

insertion of CG     assembly  G A T - - T A C
                    read      G A T C G T A C
                    kmerQ       ^ ^     ^ ^

opportunity,        read      G A T C G T A C     deleting the k-mer CG
2 bp deletion       kmerQ       ^ ^     ^ ^

opportunity,        read      G A T C G T A C     inserting into the gap between T and C
insertion           kmerQ       ^ ^ ^ ^
```

A deletion can't test the bases it removes, so neither does its opportunity
(D4 in [the review](documents/nonrepetitive_indel_review.md)); an insertion's
opportunity has no inserted bases to test, so neither does the call. An event
without two read bases on a side, at the very end of a read, is neither called
nor counted (D3). Subread support is tested on the same bases. Inside repeats,
workflow 4 tests a tract's calls and opportunity on a read over the whole
tract and two bases on each side (see [Workflow 4](#workflow-4-indels-inside-repeats)).

HiFi base qualities are capped at 93, so kmerQ ≥ 93 means every base in the
window is at the maximum. The scans write the threshold (`strata.min_kmerq`)
as the key `kmerQ_pass` on calls and opportunities alike. Workflow 4's final
filtering keeps `kmerQ_pass` = `pass` on both. The calls keep the raw `kmerQ`
as well.

Workflow 3 tests kmerQ by position (ticket 49,
artefact_descriptions/kmerq_results.md). The keys `inner_flank_q` (L1 and R1,
the bases next to the event) and `outer_flank_q` (L2 and R2) bin each pair's
lowest quality at 60 and 93, and final filtering keeps L1 and R1 at Q93 and L2
and R2 at Q60 or more, on calls and opportunities alike. An insertion call's
inserted bases must also all be at `strata.min_inserted_quality` (93): that
test is on the calls only, as an opportunity has no inserted bases, so the
insertions' rate is low by the share of real insertions it removes (s in
kmerq_results.md). The calls carry each flank base's quality and the inserted
bases' lowest.

### Haploid concordance

A real InDel is part of the read, so it should appear wherever the read is
aligned. One that appears only in the alignment to the diploid assembly is more
likely an artefact of that alignment or assembly. Each call is therefore checked
against an independent alignment of its read to the haploid assembly, in
[filtering_haploid_mapping.py](scripts/filtering_haploid_mapping.py):

1. Every read carrying a call (before final filtering) is written to FASTQ in its
   original orientation and aligned with `pbmm2 align --preset HIFI`. The check
   reads the same FASTQ for the reads' sequence.
2. Each call is joined to one record of its read: its plain forward or reverse
   alignment (flag 0 or 16), the primary, when it has one. A call whose read
   has no such record (it did not align, or aligned only as secondary or
   supplementary), or is not in the FASTQ at the call's read length, is kept
   with `NA` in the check's columns.
3. The two alignments are compared in the read as sequenced. Each InDel is an
   interval of the read: its inserted bases, or the empty interval at a
   deletion's gap. The haploid CIGAR is walked counting read bases through
   `M`, `=`, `X`, `I`, `S` and `H`. If the read lies on opposite strands in the
   two alignments, the InDel's interval is mirrored, not its anchor.
4. An InDel in a repeat can sit anywhere in the stretch of read it can slide
   over, and aligners put it at the leftmost place on their own assembly's
   strand: on opposite strands, at opposite ends. The call's stretch is worked
   out from the read and, for a deletion, the deleted bases. `haploid_check` is
   True if the haploid alignment has an InDel of the same type and size
   anywhere in that stretch, False otherwise, and `NA` where the haploid
   alignment clips the call off or spans another number of bases than the
   read.

`distance_to_indel_in_hap` is the number of read bases from the call to the
nearest matching haploid InDel: 0 at the same place, more where the two
alignments put an InDel in a repeat at different places. Without a match it is
`NA`.

The check is a column of the calls. Workflow 3's final filtering drops the
calls it fails, keeping `NA`, as a step on the calls only (ticket 49), and the
filter report counts the calls it removes and those without a check. Which of
it and the [hap1/hap2 test](#the-hap1hap2-test) a final call must pass is
ticket 43; each step can be switched off in `final_filters.nonrepetitive_calls_only`.
Tables made before ticket 18 dropped the calls the check failed, with four bugs
in it that failed real InDels (H1–H4 in
[the review](documents/nonrepetitive_indel_review.md)).

### The hap1/hap2 test

A call can also be a germline difference between the two haplotypes, on a read
aligned to the wrong one in the diploid assembly. A real mutation is in the
read, so it shows against both haplotypes; a germline difference disappears
against the haplotype the read came from. So every read carrying a call is
also aligned to hap1 alone and to hap2 alone (see [Haplotypes
alone](#workflow-1-precompute)), and the same check as above, the call's
InDel anywhere in the stretch it can slide over, gives `hap1_check` and
`hap2_check`: True, False, or `NA` where it cannot be made (no flag 0 or 16
record, the call clipped off). True in both is what a mutation looks like;
True in one only, a germline difference. Like the haploid check it is a
column, and workflow 3's final filtering drops the calls False in either,
keeping `NA`, on the calls only (ticket 49).

### Distance to the read end

`ToreadEnd` is the distance from the anchor to the nearer end of the read's
alignment to the diploid assembly, in assembly bases (soft-clipped bases are not
counted). The key `read_end_distance_bin` bins it for calls and opportunities
alike; in workflow 4 it is the distance from the tract, the nearer flank.
Workflow 4's final filter keeps the bins from 100 bp on. Workflow 3 records it
without filtering; [nonrepetitive_indels.ipynb](plotting/nonrepetitive_indels.ipynb)
can apply a threshold (`min_read_end_distance`).

### The unique flag

During the scan, before any other filter, each position counts the reads that
have an InDel of any size anchored there or a deletion spanning it. A call is
`unique` if at most one read qualifies, and `non-unique` otherwise. The
workflows keep both. The notebooks drop non-unique calls by default
(`unique_only`).

### Final filtering and the filter report

There is one definition of a final call (Q7), applied in the pipeline by
[final_filtering.py](scripts/final_filtering.py) to a sample's calls and its
opportunities alike, so a filter never lowers the rate by what it removes from
the calls alone. DuckDB filters the tables, the opportunity cube as Parquet and
the calls as TSV, spilling to disk past the job's memory. Its steps, in order,
each on both tables:

1. **Capped SMRT Cells** (Q12, Q17). A SMRT Cell none of whose reads has a base
   quality above `final_filters.capped_quality` (40) had its qualities capped
   there; its calls and opportunities are dropped. Precompute's
   `smrt_cells.tsv` gives each cell's highest base quality, over the sample's
   first `precompute.smrt_cell_reads` alignments, so every table of a sample
   drops the same cells; a cell not among those reads is kept. The report's
   note lists the cells dropped.
2. **The kept values.** Each key under `final_filters.keep` (both workflows),
   then `final_filters.nonrepetitive_keep` (workflow 3) or
   `final_filters.repeat_keep` (workflow 4), keeps the rows whose value is one
   of those listed. Workflow 3: kmerQ by position (`inner_flank_q` = `>=93`,
   `outer_flank_q` from `60-93` on) and `base_after_masked` = `no`. Workflow 4:
   `kmerQ_pass` = `pass` and `read_end_distance_bin` from `100-300` on. A
   filter on a key that is `NA` throughout the sample, a tag its platform
   lacks (Sequel II reads have no subread tags, and mbewe's and samson's no
   `np`), is skipped for it and noted (Q15). Nothing is filtered on the read
   end in workflow 3: `read_end_distance_bin` is a key to filter on later.
3. **On the calls only** (workflow 3, `final_filters.nonrepetitive_calls_only`,
   ticket 49): an insertion's inserted bases (`inserted_bases_pass`), the
   haploid check (`haploid_check`) and the hap1/hap2 test (`hap1_check`,
   `hap2_check`). These are properties of the call's read, which the
   opportunities don't have (Q7), so each lowers the rate by the real calls it
   removes too, and the report says so. `NA`, a call that can't be checked, is
   kept. Which haplotype check a final call must pass is ticket 43.

Every row of both outputs gets a `platform` column: `Sequel II` for the samples
under `sequel_ii_samples`, else `Revio`. Changing a filter or a threshold is a
change to `final_filters` and a rerun of this step: a labelled sensitivity
analysis.

`filter_report.tsv` has one row per step (tab-separated): the tables it applies
to (`both`, or `calls`), the calls and the opportunities (summed counts) before
and after it and removed by it, the calls left without a haploid check, and a
note:

```
filter                     applies_to  calls_before  calls_after  calls_removed  opportunities_before  opportunities_after  opportunities_removed  haploid_check_missing  note
input                      both        …             …            0              …                     …                    0                      …
capped_smrt_cells          both        …             …            …              …                     …                    …                      …                      the cells dropped
inner_flank_q=>=93         both        …             …            …              …                     …                    …                      …
haploid_check=True,NA      calls       …             …            …              …                     …                    0                      …                      calls only: …
```

`haploid_check_missing` counts the calls left after the step that have no
haploid check (`NA`): their read had no usable haploid alignment, or the
alignment clipped the call off. Values the callers write as `NA` stay `NA` in
the output.

### What the opportunities go through

The opportunities pass the scan's filters (coverage, mask or tract, and the
homopolymer rule in workflow 3) and the kmerQ threshold. The read-end distance
and the unique flag remove calls only, with nothing matching them in the
denominator. Since they remove some real mutations along with the
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
├── precompute/{species}/{sample}/
│   ├── contigs/{contig}.fa
│   ├── long_contigs.bed
│   ├── avg_coverage.txt
│   ├── unresolved_haplotypes/{contig}.bed        per-window pieces in chunks/
│   ├── haplotypes/                               hap1.fa, hap2.fa
│   ├── smrt_cells.tsv                            each SMRT Cell's highest base quality
│   ├── sorted_tissue_classification.tsv          tissue_samples only
│   ├── read_classification.tsv                   methylation_samples only
│   └── read_bins.npz                             the same, compact, for the scans; methylation_samples only
├── repetitive_dna/{species}/{sample}/
│   ├── per_k/                                    per period (and per window for M1)
│   ├── {contig}_{model}.bed                      overlaps resolved
│   ├── bases/{contig}_{model}.bed                the bases in any tract, before resolution
│   ├── {model}.bed                               all contigs
│   └── excluded/                                 {contig}_homopolymers.bed, {contig}_long_periods.bed
│                                                 and their union {contig}.bed; TRF windows in chunks/
├── nonrepetitive_indels/{species}/{sample}/
│   ├── mask/{contig}.bed
│   ├── calls/  opportunities/  fastq/            per contig; per window in chunks/
│   ├── haploid/                                  all.fastq.gz, remapped.bam, locations.tsv,
│   │                                             {haplotype}.bam, {haplotype}_locations.tsv
│   └── merged/                                   indels.tsv, opportunities.parquet, before final filtering
└── repetitive_indels/{species}/{sample}/
    ├── fastq/  haploid/                          shared by the models
    └── {model}/
        ├── tracts/{contig}.bed                   the tracts measured
        ├── calls/  opportunities/
        └── merged/                               indels.tsv, opportunities.parquet

results/
├── precompute/{species}/{sample}/unresolved_summary.tsv
├── repetitive_dna/{species}/{sample}/
│   ├── {model}_summary.tsv
│   ├── {model}_dropped.tsv
│   └── excluded_summary.tsv
├── nonrepetitive_indels/{species}/{sample}/
│   ├── indels.tsv
│   ├── indels_classified.tsv                     methylation_samples only
│   ├── opportunities.parquet
│   └── filter_report.tsv
└── repetitive_indels/{species}/{sample}/
    ├── tracts_summary.tsv
    └── {model}/                                  the same four files

plots/
├── nonrepetitive_indels/{species}/{sample}/plots.pdf
└── repetitive_indels/{species}/{sample}/{model}/plots.pdf
```

The scans' opportunity tables and the FASTQs are gzipped: they are big (one
1 Mb window of opportunities is hundreds of MB as text, about 20 times less
gzipped). From the merge on, the opportunities are a cube in Parquet (see
[Opportunities](#opportunities)), about a third the size of gzipped text.

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
| `context` | workflow 3 only: the opportunity context the call is counted against. For a deletion `canonical_indel`; for an insertion the anchor and the read base after the inserted ones, grouped with their reverse complement. |
| `read_length` | length of the read |
| `np` | number of passes (PacBio `np` tag) |
| `ec` | effective coverage (PacBio `ec` tag) |
| `rq`, `fn`, `rn` | predicted accuracy, and full passes on the forward and reverse strands of the molecule (PacBio tags); `NA` without them |
| `mapQ` | the read's mapping quality to the diploid assembly |
| `sm`, `sx` | subreads matching / mismatching at the read base after the anchor: the first inserted base, or the first base after a deletion (PacBio `sm`, `sx`). PacBio stores these tags in the order the read was sequenced, which runs backwards along the alignment of a read aligned to the reverse strand; they are read in the alignment's order (S1). |
| `sm_pre`, `sx_pre` | the same at the anchor |
| `sa` | subread coverage at the anchor, decoded from PacBio's run-length-encoded `sa` tag |
| `ToreadEnd` | distance from the anchor to the nearer end of the read's alignment |
| `kmerQ` | lowest base quality in the window around the event, see [kmerQ](#kmerq) |
| `position` | 1-based coordinate of the anchor on the contig |
| `readID` | read name |
| `diploid_direction` | strand of the read's alignment to the diploid assembly |
| `flank_sequence` | the assembly from `scan.flank_width` bases before the anchor to as many after it, for comparing sites between individuals |
| `canonical_motif`, `k`, `copies`, `trf_purity`, `model_purity`, `model` | workflow 4 only: the tract the call is counted in |
| `indel_class`, `size` | workflow 4 only: `slippage` or `other` (see [InDel size](#indel-size)), and the InDel's length |
| `tract_copies_bin`, `tract_purity_bin`, `left_flank_bin`, `right_flank_bin` | workflow 4 only, keys: see [Keys](#keys) |
| the keys | `tissue`, `kmerQ` and the other keys, the same columns as in the opportunity table: see [Keys](#keys) |
| `qual_L2`, `qual_L1`, `qual_R1`, `qual_R2` | workflow 3 only (ticket 49): the base quality of the two read bases before the event and the two after it, the bases kmerQ tests |
| `sx_L1`, `sx_R1` | workflow 3 only: subreads mismatching at L1 and R1, the read bases next to the event (the anchor, and the first base after the deletion or after the inserted bases); `NA` without subread tags |
| `inserted_min_quality`, `inserted_bases_pass` | workflow 3 only: an insertion's lowest inserted base quality, and `pass` if it reaches `strata.min_inserted_quality`, else `fail`; `NA` for a deletion |
| `mask_distance` | workflow 3 only: bases from the event (the anchor to the base after) to the nearest masked base: 0 if the base after a deletion is masked, 1 if a masked base is right next to it; `NA` if none is within 1,000 bp |
| `bitflag_hap`, `contig_hap`, `pos_1based_hap`, `map_mapq` | the read's alignment to the haploid assembly: flag, contig, position and mapping quality |
| `haploid_check` | whether the same InDel is in the haploid alignment, anywhere it can slide; `NA` where that cannot be checked. See [Haploid concordance](#haploid-concordance) |
| `distance_to_indel_in_hap` | read bases from the call to the matching haploid InDel (0 means the same place; more where the two alignments put an InDel in a repeat at different places). `NA` without a match. |
| `softclip_haploid` | whether the haploid CIGAR has a soft clip; `NA` without a usable haploid alignment |
| `other_ops` | `HPN` if the haploid CIGAR has `H`, `P` or `N` operations, else `NA` |
| `hap1_check`, `hap2_check` | whether the same InDel is in the read's alignment to hap1 alone, and to hap2 alone; `NA` where that cannot be checked. See [The hap1/hap2 test](#the-hap1hap2-test) |
| `species`, `sample` | |

### Opportunities

The opportunities are a cube (ticket 40): the read positions that could carry
an InDel, counted for every combination of contig, InDel type, size, context
and the keys, so the count for any set of them is a filter and a sum.
Workflow 3's `opportunities.parquet` has one row per combination seen
(about 1.5 billion for a Revio sample), `size` and `count` as numbers and the
rest as text:

| Column | Meaning |
|---|---|
| `contig` | |
| `indel_type` | `deletion` or `insertion` |
| `size` | the InDel size it counts for: 1, 2 or 3 |
| `context` | for a deletion, the canonical k-mer: the smaller of the k-mer and its reverse complement, so 1 bp contexts are only `A` and `C`. For an insertion, the two read bases around the gap, grouped the same way. |
| the keys | the same columns, with the same values, as at the end of each call: see [Keys](#keys) |
| `count` | read positions |
| `species`, `sample` | |

Workflow 4's have `canonical_motif`, `k` (a number) and the tract keys in place
of `indel_type`, `size` and `context`. DuckDB reads only the columns a query
needs, so a lookup over a whole sample takes seconds. In R, with dplyr:

```r
library(duckdb); library(dplyr)
con <- dbConnect(duckdb())
opportunities <- tbl(con, sql("SELECT * FROM read_parquet('results/nonrepetitive_indels/*/*/opportunities.parquet')"))
opportunities |>
  filter(mapQ_bin == ">=60", read_end_distance_bin != "<25") |>
  group_by(species, sample, size) |>
  summarise(opportunities = sum(count, na.rm = TRUE)) |>
  collect()
```

and in Python:

```python
import duckdb
duckdb.sql("""
    SELECT species, sample, size, sum(count) AS opportunities
    FROM 'results/nonrepetitive_indels/*/*/opportunities.parquet'
    WHERE mapQ_bin = '>=60' AND read_end_distance_bin <> '<25'
    GROUP BY ALL""").df()
```

The calls carry the same keys, so the same filter on them gives the rate.

### Keys

Every feature a call might be filtered on is a key, defined once in
[scripts/strata.py](scripts/strata.py). Both scans write each key as a column
at the end of every call, and count the opportunities per value of it, so a
filter on a key has a matching denominator. Adding a key there adds it to both
tables.

| Key | Meaning |
|---|---|
| `tissue` | the read's tissue label, or `NA` |
| `methylation_bin` | the read's germline/somatic bin, one of the 14 of the `classification` block, from precompute's `read_classification.tsv` (read by the scans as `read_bins.npz`); `NA` for reads without a score and for samples not under `methylation_samples` |
| `read_length_bin` | the read's length, binned by `strata.read_length_bins` |
| `aligned_strand` | the strand of the assembly the read is aligned to: `+`, or `-` for flag 16 |
| `smrt_cell` | the SMRT Cell (movie): the first part of the read's name |
| `kmerQ_pass` | `pass` if every base tested has at least `strata.min_kmerq`, else `fail`; final filtering keeps `pass` on both tables. Calls also carry the raw `kmerQ`. |
| `read_end_distance_bin` | distance from the anchor to the nearer end of the read's alignment, binned by `strata.read_end_distance_bins`. An opportunity's anchor is the base before its k-mer or gap. |
| `mapQ_bin` | the read's mapping quality to the diploid assembly, binned by `strata.mapq_bins` (N1); in workflow 3 by `strata.nonrepetitive_bins.mapq_bins` |
| `np_bin`, `rq_bin` | the read's passes (`np`) and predicted accuracy (`rq`), binned by `strata.passes_bins` (in workflow 3 `strata.nonrepetitive_bins.passes_bins`) and `strata.accuracy_bins`; `NA` without the tag |
| `strand_passes_bin` | full passes on the weaker strand of the molecule, the smaller of `fn` and `rn`, binned by `strata.strand_passes_bins`; `NA` without the tags. Workflow 4 only since ticket 49 (`strata.nonrepetitive_without`); workflow 3's calls keep `fn` and `rn` |
| `nearby_event_bin` | distance from the anchor to the read's nearest other InDel or mismatch with the assembly, binned by `strata.nearby_event_bins` (in workflow 3 `strata.nonrepetitive_bins.nearby_event_bins`); `none` if the read has none. Events at the anchor itself, such as the call's own InDel, don't count. |
| `subread_support` | lowest share of covering subreads that match (`sm`/`sa`), among the bases tested, binned by `strata.subread_support_bins`; `NA` for reads without subread tags (Sequel II) |
| `inner_flank_q`, `outer_flank_q` | workflow 3 only (ticket 49): the lowest base quality of L1 and R1, the read bases next to the event, and of L2 and R2, the two outside them, binned by `strata.flank_quality_bins` (`<60`, `60-93`, `>=93`). Final filtering keeps `inner_flank_q` = `>=93` and `outer_flank_q` from `60-93` on: kmerQ by position, see [kmerQ](#kmerq) |
| `junction_disagreement` | workflow 3 only: `yes` if a subread mismatches the read at L1 or R1, else `no`; `NA` without subread tags |
| `repeats_next` | workflow 3 only: for a deletion, `yes` if the deleted bases repeat right before or right after them (as AC in ACAC), else `no`; from the assembly for a call, which lacks them, and from the read for an opportunity. `NA` for an insertion |
| `base_after_masked` | workflow 3 only: `yes` if the base after a deletion is masked. The mask test covers the anchor and the deleted bases, not the base after; the class the expected rates are for excludes these deletions. Always `no` for an insertion, whose mask test covers it. Final filtering keeps `no` |
| `tract_copies_bin`, `tract_purity_bin` | workflow 4 only: the tract's copies and purity (`model_purity`), binned by `strata.tract_copies_bins` and `strata.tract_purity_bins` (`>=100` is a pure tract) |
| `left_flank_bin`, `right_flank_bin` | workflow 4 only: how many bases of the read's alignment lie before the tract and after it, binned by `strata.read_end_distance_bins`; at least `scan.min_tract_flank` |

In workflow 4 the site keys are those of the tract in the read: kmerQ and
subread support over the tract and two bases on each side, the distance to the
read end from the tract, and the nearest event other than the tract's own
InDels.

Every read position is counted once for each size of deletion and of
insertion. The calls' `context` uses the same form, so the rate of a context
is its calls divided by its opportunities, joined on `indel_type`, size (the
length of the calls' `indel_string`), `context` and the keys. Insertions have
their own opportunities, at the gaps between read bases, rather than the
deletions' (D1 in [the review](documents/nonrepetitive_indel_review.md)).

Workflow 4's table has `canonical_motif` and `k` (the tract's) in place of
`indel_type`, `size` and `context`, then the tract keys, and `count` is tracts
in reads that span them (R6 in [the workflow 4
review](documents/workflow4_review.md)), so rates by tract length and purity
need no approximation.

### Other tables

- `{model}_summary.tsv`: `canonical_motif`, `k`, `copies`, `trf_purity`,
  `model_purity`, `model`, `species`, `sample` and `count`, the number of tracts
  of each kind.
- `{model}_dropped.tsv`: one row, the `model`, with `bp_before`, the bases in a
  tract of the model before overlap resolution; `bp_kept`, those in the tracts
  it keeps; `bp_dropped`, the difference, the bases of dropped tracts that no
  kept tract covers, which are in neither InDel workflow; and `share_dropped`,
  `bp_dropped` over `bp_before`; then `species` and `sample`.
- `excluded_summary.tsv`: one row per `class` of the excluded group,
  `homopolymer`, `long_period`, `long_period_timed_out` (the sub-windows TRF
  didn't finish, excluded whole, part of `long_period`) and `excluded` (the
  union), with `regions`
  (homopolymer runs, or merged intervals), `bp`, and `share`, the bases over
  those of the contigs analysed; then `species` and `sample`.
- `tracts_summary.tsv` (workflow 4): one row per `model`, with `tracts` and `bp`,
  the model's tracts from workflow 2 and their bases, and `tracts_measured` and
  `bp_measured`, those left once the tracts overlapping a homopolymer run, an
  unfinished long-period TRF sub-window or an unresolved-haplotype region are
  dropped; then `species` and `sample`.
- `unresolved_summary.tsv` (workflow 1): one row, with `regions`, the merged
  unresolved-haplotype regions; `bp`, their bases; `contigs_bp`, the bases of
  the contigs analysed; and `share`, `bp` over `contigs_bp`, the share of the
  genome both InDel workflows leave out for them; then `species` and `sample`.
- `filter_report.tsv`: `filter`, `applies_to`, `calls_before`, `calls_after`,
  `calls_removed`, `opportunities_before`, `opportunities_after`,
  `opportunities_removed`, `haploid_check_missing`, `note`, see
  [the report](#final-filtering-and-the-filter-report).
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
git config filter.strip-notebook-outputs.smudge cat
git config filter.strip-notebook-outputs.required true
```

Without the first line, git commits notebooks with their outputs. The second
checks notebooks out as they are; a required filter needs it, or git can't
check out a notebook (in a new worktree, say). The third makes `git add` stop
with an error if the filter fails, for instance on a notebook that isn't valid
JSON, rather than staging it unfiltered.

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
| [tests/test_pure_functions.py](tests/test_pure_functions.py) | canonicalisation, the homopolymer rule, reading the `sa` tag, Hamming purity |
| [tests/test_haploid_check.py](tests/test_haploid_check.py) | the haploid check: findings H1–H4 of the review on its CIGAR cases, InDels in a repeat on flipped reads, the calls it cannot check kept with `NA`, the primary record chosen, and the filter report's count of unchecked calls |
| [tests/test_haplotype_check.py](tests/test_haplotype_check.py) | the hap1/hap2 test: a planted mutation shows against both haplotypes and a germline difference against one, on either strand; precompute's split of the diploid assembly |
| [tests/test_strata.py](tests/test_strata.py) | the stratification keys: one set, in [strata.py](scripts/strata.py), columns of calls and opportunities in both scans |
| [tests/test_subreads.py](tests/test_subreads.py) | subread support (S1): the tags read from the other end on reverse-strand reads, and the binned support of the bases tested |
| [tests/test_read_geometry.py](tests/test_read_geometry.py) | the read-geometry keys: read length, distance to the read end, aligned strand, SMRT Cell, and kmerQ pass or fail |
| [tests/test_read_quality.py](tests/test_read_quality.py) | the read-quality keys: diploid mapQ (N1), passes, predicted accuracy, full passes on the weaker strand, the nearest other event; the assembly around each call |
| [tests/test_methylation_bins.py](tests/test_methylation_bins.py) | the germline/somatic bin as a key of calls and opportunities, and the per-read table binning reads as the classification scripts do |
| [tests/test_read_bins.py](tests/test_read_bins.py) | the compact read bins: every read's bin as the classification table gives it, `NA` for a read missing from it, 9 bytes a read, more bins than a byte holds refused, and workflow 1's target writing them |
| [tests/test_final_filtering.py](tests/test_final_filtering.py) | the final filtering stage: a filter on a feature spread alike in calls and opportunities leaves the rate unchanged, capped SMRT Cells (from their reads' base qualities) dropped and listed, a filter on a missing tag skipped and noted, the platform on every row, the report; Parquet, plain and gzipped tables |
| [tests/test_opportunity_tables.py](tests/test_opportunity_tables.py) | the opportunity cube: a contig's windows summed per key into Parquet, sorted by key as text, counts past 4 billion exact, `size`, `k` and `count` as numbers, an empty cube for a contig without opportunities, a sample's contigs joined with its species and sample, nothing left behind when a sum or a join fails, a count that isn't a whole number, a key with a NUL byte or windows with different columns failing the sum |
| [tests/test_tract_unit.py](tests/test_tract_unit.py) | workflow 4's unit: one tract per spanning read, reads ending in the tract or its flank counting for neither, the tract keys on both tables, calls and opportunities of a tract tested on the same bases, overlapping tracts refused, windows adding up by tract; the slippage and other classes |
| [tests/test_tract_edges.py](tests/test_tract_edges.py) | InDels at a tract's edge (L1, R1, R2): slippage anchored before a tract or a partial unit is the tract's, slippage of the next tract is that tract's, every event near a tract counted by one workflow, opportunities tested on the same bases, windows still adding up |
| [tests/test_denominators.py](tests/test_denominators.py) | denominators that match the calls (D1, D3, D4): the bases calls and opportunities are tested on, insertion opportunities per gap, events at read ends |
| [tests/test_classification_bins.py](tests/test_classification_bins.py) | the two classification scripts binning reads the same way |
| [tests/test_compare_tables.py](tests/test_compare_tables.py) | the mirror replay's table comparison: by content, whatever the row order, number format or appended columns |
| [tests/test_integration.py](tests/test_integration.py) | the callers on a small synthetic assembly and BAM with planted InDels: the right reads, positions, types and sizes; masking; windows adding up to the whole contig; one scan for three models matching three separate scans; the haploid check on several call tables |
| [tests/test_read_shapes.py](tests/test_read_shapes.py) | the scans on the read shapes real data has: reverse strand (with subread tags in sequencing order), base qualities below 93, reads ending early or crossing a window edge, soft clips, mapQ below 60, overlapping tracts |
| [tests/test_config.py](tests/test_config.py) | config.yaml: every sample belongs to a known species, override BAMs are well formed, the layout templates resolve, and every key is read by a workflow or a notebook |
| [tests/test_workflow_snapshots.py](tests/test_workflow_snapshots.py) | every target each workflow would run (name, inputs, outputs, options, command) for a fixture config, against the snapshots in `tests/snapshots/`; workflow 4's step 0 reading the long-period window BEDs workflow 2 writes; workflow 3 merging each contig in its own job; both workflows summing the opportunities in jobs of their own |
| [tests/test_workflow_helper.py](tests/test_workflow_helper.py) | `workflows/pipeline.py`: paths from the layout templates, alignment overrides, the long contigs; the joins of calls and reads, which stop the job at a window they can't read and leave no file behind |
| [tests/test_unresolved_haplotypes.py](tests/test_unresolved_haplotypes.py) | the unresolved-haplotype regions: the evidence rule, and regions that overlap or touch merging to their union |
| [tests/test_masks.py](tests/test_masks.py) | what the InDel workflows leave out, from the workflows' own targets run with bash on a one-contig sample (TRF and samtools are stand-ins, from [tests/targets.py](tests/targets.py)): the excluded group, workflow 3's mask, the tracts workflow 4 measures (dropped for a homopolymer run, measured over a long-period tract), and their summaries |
| [tests/test_tract_drop_causes.py](tests/test_tract_drop_causes.py) | why workflow 4 drops tracts: each tract counted once under its first cause, and the tracts each rule ticket 34 weighed would measure, option 2 matching the workflow's own step 0 |
| [tests/test_mirror_export.py](tests/test_mirror_export.py) | the mirror exporter and the replay, end to end on a fake cluster checkout, one of its samples with read bins, which a window exported before them makes from the text table |
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
opportunities compared too. The script exits 1 if any table differs. Workflow
4's scan reads the window's measured tracts, or, in a window exported before
workflow 4 dropped tracts (step 0), all of the model's tracts, as the cluster
scanned them then.

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
  parse_trf.py                     TRF output -> repeat BED (M1, the long periods); canonicalise(), Hamming purity
  resolve_repeat_overlaps.py       one tract per base per model
  find_homopolymers.py             homopolymer runs of a FASTA (the excluded group)
  merge_regions.py                 the merged union of BEDs: the mask, the excluded group
  drop_excluded_tracts.py          the tracts workflow 4 measures
  tract_drop_causes.py             why workflow 4 drops tracts, under each rule ticket 34 weighed
  filtering_haploid_mapping.py    haploid concordance check
  final_filtering.py               the definition of a final call, on calls and opportunities, with the report
  opportunity_tables.py            the opportunity cube: a contig's windows summed, a sample's contigs joined
  smrt_cell_qualities.py           each SMRT Cell's highest base quality (the Q40 rule)
  read_bins.py                     every read's germline/somatic bin, compact, for the scans
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
