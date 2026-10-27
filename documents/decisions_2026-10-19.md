# Decisions, 19 October 2026

These are the decisions from a structured review of the whole project on 19 October 2026. Each decision carries its reason and the number of the question that settled it (Q1–Q43), so it can be traced and reopened deliberately rather than by drift.

Bug labels (S1, L1, D1–D4, H1–H4, U1, F1, N1, M1, E1) refer to [nonrepetitive_indel_review.md](nonrepetitive_indel_review.md). Floor figures come from [rate_floor.md](rate_floor.md).

This document replaces `overview.md`, `InDel_call_plan.md` and `tandem_repeat_plan.md`, which contradict the code and are to be deleted (Q39). Where this document and the README disagree, this document is the intent and the README is out of date until it's rewritten after the rerun.

## Taken as given

The foundations are not under review (Q4):
- PacBio HiFi reads;
- each individual's reads aligned to that individual's own diploid assembly;
- calls from the pysam pileup;
- gwf as the workflow engine.

Reopen one only if it's tied to a known problem.

## The goal

**The paper is about germline vs somatic InDel rates** (Q5, Q9). Repeat context, species and age are breakdowns of that comparison, not separate headline results.

**The current goal is to lower the false-positive rate.** No germline or somatic claim is made until it's low enough (Q11).

### Why the false-positive rate comes first

A measured rate is the true rate plus the rate of false calls that survive every filter, the "floor". As of 7 October:

| Quantity | Value per bp | Source |
|---|---|---|
| Summed rate outside repeats, median Revio sample, before filtering | 1.6×10⁻⁷ | rate_floor.md |
| After the cleanest-subset filters | 8.7×10⁻⁸ | rate_floor.md |
| Estimated with S1 fixed | 5.3×10⁻⁸ | rate_floor.md |
| Expected germline rate per sperm, InDels outside tandem repeats | about 1–2×10⁻⁹ | literature, below |

With a floor of about 5×10⁻⁸, about 97% of calls in a germline sample are false. Germline and somatic reads in the same testis samples currently give the same rate (E1), which is what a floor-dominated measurement looks like.

Subtracting the floor doesn't work either, because the floor isn't constant: Sequel II runs about 2.6× higher than Revio, and SMRT Cells within a sample differ by up to about 3×.

### The gate (Q11)

| Claim | Floor needed |
|---|---|
| Germline or somatic rates, or germline vs somatic | ≤ about 1×10⁻⁹ per bp; about 1×10⁻¹⁰ keeps false calls under 10% |
| Sperm-age effect | Below about 4×10⁻¹⁰ per bp, the expected change over ages 25–38, or estimated per sample well enough to subtract |

The uncited target of 1×10⁻⁹ in `rate_floor.md` is replaced by these per-class values. 1×10⁻⁹ is the size of the germline signal, not a floor low enough to measure it.

### Evidence that the floor is low (Q14b, Q30)

The direct test, requiring both strands of the molecule to carry the InDel, is impossible with these data. The only subread information is the pooled `sa`/`sm`/`sx` tags: no subreads and no by-strand consensus were kept (Q14). No calibration run will be made.

The evidence is two indirect tests:
1. **The hap1/hap2 mismapping test.** Each call-carrying read is aligned to hap1 and to hap2 separately, and a call counts only if it appears in both. A real mutation is in the read, so it shows against both. A germline difference disappears against the haplotype the read came from.
2. **The rate plateaus as filters tighten**, with sensitivity tracked by the control calls (the same InDel in ≥ 3 reads).

The paper must say plainly that no strand-level test was possible. Agreement between SMRT Cells was dropped as a proof because outlier cells are out of scope (Q26, Q30).

The floor must not be "proven" by sperm rates landing on the literature value. That would be circular.

### Literature values

These come from a literature search on 19 October 2026. Values marked *derived* are arithmetic on published figures. Check every value against the paper before citing it in the manuscript.

- **Germline InDel rate, human, per bp per generation:**
  - 9.3×10⁻¹⁰ (95% CI 8.8–9.8×10⁻¹⁰), InDels < 35 bp. Besenbacher et al. 2016, *PLoS Genet*, doi:10.1371/journal.pgen.1006315.
  - About 0.9×10⁻⁹, *derived* from 6.8% of de novo mutations being InDels. Jónsson et al. 2017, *Nature*, doi:10.1038/nature24018.
  - 7.4 InDels outside tandem repeats per transmission, so about 1.4×10⁻⁹ per bp per gamete (*derived*). Porubsky et al. 2025, *Nature*, doi:10.1038/s41586-025-08922-2.
  - No published rate is restricted to 1–3 bp outside repeats.
- **Sperm, direct measurement:** 0.10 InDels per year per haploid genome (95% CI 0.06–0.15), about 3×10⁻¹¹ per bp per year (*derived*). Neville et al. 2025, *Nature*, doi:10.1038/s41586-025-09448-3.
- **Somatic, per diploid genome per year:**
  - Blood stem cells: 0.71 InDels, about 6×10⁻⁹ per bp at age 50 (*derived*). Mitchell et al. 2022, *Nature*, doi:10.1038/s41586-022-04786-y.
  - Neurons: 3.1 InDels. Abascal et al. 2021, *Nature*, doi:10.1038/s41586-021-03477-4.
  - These cover all contexts. The share outside repeats is not published.
- **Error floors other methods reach:**
  - NanoSeq (Illumina duplex): InDel error < 3×10⁻⁹, an upper bound (Abascal et al. 2021).
  - HiDEF-seq (PacBio): leaves InDels out, stating it does not reach single-molecule fidelity for them.
  - No PacBio HiFi method has published an InDel floor.

## Samples and platforms

- **Revio is the reference platform (Q12, Q15).** Filters are developed and validated on Revio samples, using all their tags.
- **Sequel II samples stay in the analysis:** carl, stephan, mbewe, samson, da1 and ob2v.
  - No validation is based on them.
  - They get only the filters their tags support. They carry no `sa`/`sm`/`sx`, and mbewe and samson also lack `np`.
- **Q40 rule (Q12, Q17):** SMRT Cells whose base qualities top out at Q40 are dropped. The rule is applied per SMRT Cell, so a sample loses only its capped cells.
  - la4 loses every call at the base-quality 93 filter (4.7M to 0, `results/nonrepetitive_indels/human/la4/filter_report.tsv`), so it goes entirely.
  - tsaed keeps 0.1% of its calls at that filter against 0.8% for chk, which suggests some of its cells are capped.
- **Outlier SMRT Cells are out of scope for now (Q21, Q26).** All cells are kept and treated alike; none are dropped or investigated for running high. The test to run if this comes back into scope is parked in `.scratch/false-positive-reduction/issues/01-diagnose-outlier-smrt-cells.md`.
- **Gorilla override BAMs (Q42).** mbewe and samson use `*_plus_unmapped_no_clips_diploid_v2.bam`. They are treated like the other Sequel II samples. If "no clips" means soft clips were removed, their read-end distances aren't comparable with other samples'. How these BAMs were made is not recorded (see *Open*).
- **The mirror is for software testing only.** `mirror/` holds a small subset of real data. It's used to test code and to see the technical shape of the data. No data conclusions are drawn from it.

## Germline and somatic reads

- **Labels decide a read's class first (Q10).** Sperm samples and sperm-labelled reads are germline. Reads from non-gonadal tissues are somatic.
  - The methylation classifier is used only for testis reads, where tissue labels can't separate germ cells from somatic cells.
  - Testis reads without methylation scores belong to neither class. That covers the tissue-labelled samples, since no sample has both kinds of information.
  - Reason: labels have no classifier error and define exact read sets, so their denominators are exact.
- **Exact per-bin denominators (Q18).** Opportunities carry each read's methylation bin, all 14 bins, alongside its tissue label. This replaces the approximation in `workflows/classification_denominator.py`, which multiplies the total opportunities by each bin's share of reads or CpGs, counted over all reads rather than the reads that pass the filters.
- **Methylation thresholds are parked (Q25).** The notebook calls a read germline at a score ≥ −1; the collaborator's notebook calls it somatic at ≤ −4. With the bins stored, the cutoffs can be chosen after the floor gate without a rerun.
- **The sperm-age effect is exploratory (Q20).** It's reported with its statistical power stated.
  - Revio human sperm with ages: chk (25), la (27), ph (35), tr (38) and possibly tsaed (35).
  - da1 (49) is Sequel II, la4 (24) is dropped, and aa_sperm has no age.
  - The expected change over 25–38 is about 4×10⁻¹⁰ per bp, too small for 4–6 donors.

## What counts as a repeat

- **A tandem repeat has a period of 2–6 (Q8).** Each base falls into exactly one of three groups:

  | Group | What's in it | Measured by |
  |---|---|---|
  | Repeat | Tracts with periods 2–6 | Workflow 4 |
  | Excluded | Homopolymers and tracts with periods ≥ 7 | Neither; masked from both workflows (since ticket 34, workflow 4 measures a period 2–6 tract that overlaps a long-period tract, see Q16) |
  | Non-repetitive | Everything else | Workflow 3 |

  Reason: homopolymers are where HiFi makes the most errors, and repeats with periods ≥ 7 still slip faster than background. Leaving them in the non-repetitive group inflates the background rate (M1). They are an explicit scope limit of the paper.
- **The excluded group is defined as (Q16):**
  - homopolymer runs ≥ 6 bp, matching the shortest period-2 tract M0 masks;
  - tracts with periods ≥ 7 from TRF with minimum score 50, as for M1.
  - Shorter homopolymer runs are still handled by the existing homopolymer rule in workflow 3.
  - **A period 2–6 tract that overlaps a long-period tract is measured (ticket 34, 27 October 2026).** Workflow 4 drops a tract whole only for a homopolymer run, a sub-window the long-period TRF couldn't finish (most likely a satellite array) or an unresolved region: the shorter period wins, as in overlap resolution. Workflow 3 still masks the long-period bases, so none is counted by both workflows. Reason: dropping every tract that touched the excluded group cost chk 22% of its M0 tracts, 44% of M1's and 18% of M2's (da1 more, with more unresolved regions). Long-period tracts were the largest cause the rule could avoid (9–17% of tracts). Ignoring overlaps of up to 5 bases at a tract's edge would have kept 1–3% more of the M0 and M2 tracts and 6–8% more of M1's, but almost none of M1's bases, and it would bring homopolymer bases at a tract's ends into workflow 4, so such overlaps still drop a tract. Measured with `scripts/tract_drop_causes.py` on chk and da1. M1 still keeps only about 30% of its bases: homopolymer runs inside its tracts drop half of them.
- **M1 (TRF) is the main repeat model (Q28).** M0 and M2 are computed for sensitivity checks only. Reasons: TRF is the standard tool reviewers will expect, and M1 records purity. The non-repetitive mask stays the union of all three models.
- **Unresolved haplotype regions are masked in both workflows (Q19).** Before, only workflow 3 masked them. In those regions a read can align to the wrong haplotype, and in a tandem repeat that turns a real germline length difference into a false call.
- **The definition of unresolved regions stays (Q32).** Any mismatch shared by ≥ 2 reads at base quality ≥ 40 masks the reads' full span. U1 is fixed, and the masked share of the genome is recorded per sample. The definition is aggressive but on the safe side for false positives.

## Calls and opportunities inside repeats

- **One opportunity per tract per spanning read (Q22).** A read counts only if it spans the whole tract plus a flank on each side; partial reads can't reveal a change in tract length.
  - The tract's length and purity are stored on each record, so rates by tract length and purity are exact. This is the fix deferred in September.
  - The flank length on each side is stored as a binned column, so the flank requirement can be tuned after the rerun.
- **Every InDel anchored in a tract is written out and classified (Q23):**
  - **Slippage:** the size is a multiple of the period, and the sequence is the motif, in any rotation, repeated.
  - **Other:** everything else, including sizes that aren't a multiple of the period. Today those are counted nowhere.

  Each class gets its own rate, and nothing is dropped silently.

## Calls outside repeats

**Sizes are 1–3 bp (Q29).** HiFi errors and most small InDels fall in this range, and widening it would add sizes with no denominators and no review. Comparisons with published rates must be restricted to 1–3 bp or adjusted for the size distribution.

## Final filtering

- **There is one definition of a final call, applied in the pipeline (Q7).** The notebooks read only the pipeline's filtered output. Changing a threshold is a clearly labelled sensitivity analysis.
- **Thresholds are tuned after the rerun (Q13).** Every candidate filter feature is a column on calls and opportunities alike, so a call filter always has a matching denominator. This closes the trap described in `rate_floor.md`: lowering the rate by filtering calls but not opportunities.
  - This includes the read-end distance, which currently appears in the documents as 0, 100, 300 and 1000 bp.

## The rerun

There is one batch of fixes and one cluster rerun (Q6, Q31). Everything that changes calls or opportunities goes into it.

**Bug fixes:**
- The review's findings that change calls or opportunities: at least S1, L1, D1, D4, H1–H4 and U1.
- **L1 is fixed on both sides.** Workflow 3 checks the anchor and the k bases after it against the mask, for calls and opportunities alike. Workflow 4 assigns those calls to the tract.
- Whatever the workflow 4 review finds (workflow 4 has never been reviewed).

**Masks:**
- homopolymers ≥ 6 bp and period ≥ 7 tracts (Q16);
- unresolved haplotype regions in workflow 4 (Q19).

**Repeats:**
- one tract per spanning read, with tract length, purity and binned flank lengths (Q22);
- slippage and other classes (Q23);
- M0, M1 and M2 all computed, with M1 as the main model (Q28).

**New columns on calls and opportunities** (Q13, Q17, Q18):
- read length bin and read-end distance bin;
- strand-corrected subread support at the k-mer (`sa`/`sm`/`sx`);
- passes per strand (`fn`/`rn`), `np` and `rq`;
- strand;
- SMRT Cell (movie);
- diploid mapQ (N1);
- other InDels or mismatches nearby in the read;
- base quality 93 as pass/fail, instead of 94 levels, to keep the opportunity table small;
- methylation bin (14 bins) and tissue label.

**New columns on calls only:**
- ±25 bp of surrounding assembly sequence (Q24). Each individual has their own assembly, so coordinates can't be compared between individuals; a sequence key can. This lets a cross-individual recurrent-site filter, like those NanoSeq and Owl use, be tested later without a rerun.
- The hap1/hap2 test result and the old haploid check result, side by side (Q27). Whether the new test replaces the old one is decided from the data.

**Final filtering stage:**
- applies the Q40 rule per SMRT Cell;
- applies to Sequel II samples only the filters their tags support;
- writes the tables the notebooks read.

## Engineering

- **Git (Q33).**
  - Commit the code, tests, `documents/`, environment files and notebooks, with notebook outputs stripped.
  - Ignore `results/`, `mirror/`, `plots/`, `figs/`, `papers/`, `TRF/`, `mreps/` and `__pycache__/`.
  - Record the TRF and mreps versions and build steps in the README.
- **Tests (Q34).**
  - Delete the 14 string-match tests in `tests/test_config.py`. Keep a short config check: every sample has a species, override paths are well formed, and every config key is used.
  - Every bug fix starts with a failing test.
  - Extend the fixtures to reverse-strand reads, sub-93 qualities, reads ending mid-tract, soft clips and overlapping tracts. Today's fixtures would pass with S1 still broken.
- **Mirror (Q35).** Extend `scripts/export_local_subset.py`:
  - several samples: a Revio sperm sample, a tissue-labelled sample, a methylation sample and a Sequel II sample;
  - smaller windows;
  - a workflow 4 replay, which doesn't exist today.
- **Shared code and config (Q36).**
  - One helper module for config loading, path building (using the `layout` block, which nothing reads today) and `.fai` parsing. It replaces the four copied preambles.
  - The species goes into every gwf target name, so two species with the same sample name can't collide (Q41).
  - Hardcoded numbers move to `config.yaml`: window sizes, periods, TRF and mreps parameters, Q93, the homopolymer length, flank sizes and the 0.5–2× coverage band.
  - The `precomute` and `filtering_haploid_mappping` misspellings are fixed at the rerun, when every output path is regenerated anyway.
- **Dead code (Q37).**
  - Delete the scripts nothing calls: `small_indels.py`, `opportunities.py`, `tandem_repeat_indels.py`, `tandem_repeat_opportunities.py`, `small_indels_summary.py`, `large_indels.py`, `impurity_models.py`, `build_mutation_fastq.py`, `recanonicalise_motifs.py`, `plotting_quality.py` and `trf_diagnose.sh`.
  - Delete the dead functions in `utilities.py`.
  - Replace `fetch_results.py` with the `rsync` command actually used, documented in the README.
- **Notebooks (Q38).**
  - The notebooks read the pipeline's final output only.
  - Their own unique, haploid-distance and read-end filters go, as does the "usable share" correction, which the read-end distance bins replace.
  - The copied R helpers, the canonicalisation variants among them, move to one shared R file.
- **Documents (Q39).** Delete `overview.md`, `InDel_call_plan.md` and `tandem_repeat_plan.md`. Rewrite the README after the rerun.
- **Excluded samples (Q41).** Each excluded sample gets a one-line reason in `config.yaml` (see *Open*).

## Order of work (Q40)

1. Commit the current state as a baseline.
2. Delete the dead code, stale documents and string-match tests.
3. Extend the mirror; Oliver runs the export on the cluster.
4. Add the shared helper and consolidate the config, **with no change in behaviour**. Confirm the mirror reproduces today's cluster outputs exactly.
5. Review workflow 4 locally.
6. For each bug, a failing test, then the fix.
7. Fix the misspelled paths, then make the one cluster rerun.
8. Update the notebooks and README, and delete `workflows/classification_denominator.py`.

The refactor comes before the fixes because, until the fixes land, the mirror can confirm the refactor changed nothing by comparing outputs byte for byte. Once the fixes are in, outputs change on purpose and that check is gone.

## Parked

- Methylation thresholds for germline and somatic (Q25).
- Which tissues serve as the somatic comparison.
- Outlier SMRT Cells (Q21, Q26), in `.scratch/false-positive-reduction/issues/01-diagnose-outlier-smrt-cells.md`.
- A cross-individual recurrent-site filter (Q24). The flank sequence is stored, but the filter isn't built.
- Germline–somatic and sperm-age claims, until the floor gate is passed (Q11).

## Open: needed from Oliver

**Reasons for excluded samples (Q41).** These go into `config.yaml`:

| Sample | Current state | Reason |
|---|---|---|
| chimp/martin | Commented out of `samples` (uncommitted change); still in `sample_age` | Removal intended; reason to add |
| long_tailed_macaque/J17017 | Commented out of `samples` and `tissue_samples` | to add |
| orangutan/kevin | Commented out of `samples`; still in `tissue_samples` | to add |
| baboon/bs5v | Not in `samples` | "contaminated" (config comment) |
| chimp/benjie | Only in `sample_age` | to add |
| human/ob003h | In `samples`, commented out of `methylation_samples` | to add |
| gibbon/gibbon_testis, gibbon/gibbon_sperm | Only in `sample_age` | to add |

**The gorilla v2 BAMs (Q42).** How `*_plus_unmapped_no_clips_diploid_v2.bam` was made for mbewe and samson, and specifically whether soft clips were removed.
