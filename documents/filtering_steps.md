# Every filter on the way to a rate outside repeats

*7 October 2026. Code at commit 542fa20, which the links to lines of code
point at; the notebook was not yet committed.
A companion to [nonrepetitive_indel_review.md](nonrepetitive_indel_review.md):
the finding IDs (L1, S1, H1, …) refer to it, and
[candidate_filters.md](candidate_filters.md) lists what could be added.*

Every step that removes data between the BAM and a plotted rate, in the order
it happens. In the tables, **Calls** and **Opps** say whether a step applies to
the calls (the numerator), to the opportunities (the denominator), or to both.
A step that touches only one side biases the rate unless something else
balances it.

## 1. Inputs

| Step | Where | Calls | Opps | Notes |
|---|---|---|---|---|
| Only samples listed under `samples` in config.yaml | every workflow | ✓ | ✓ | |
| Only contigs longer than 1 Mb (`min_contig_length`, strictly greater) | [indel_mutation_call.py:88](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/indel_mutation_call.py#L88) | ✓ | ✓ | |
| Whatever the BAM already lacks | upstream | ✓ | ✓ | Not visible from here; the gorilla override BAMs are named `plus_unmapped_no_clips`, so something was done to them |

## 2. Precomputed inputs

### Mean coverage, which sets the coverage window

`samtools depth` over the long contigs, averaged over the positions it reports
([precompute.py:82](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/precompute.py#L82)). With its defaults it:
- skips unmapped, secondary, QC-fail and duplicate reads;
- leaves out positions with no reads;
- does not count reads with a deletion at the position.

The scan's coverage test counts all of those reads (`nsegments`), so the two
depths are measured slightly differently.

### Unresolved-haplotype regions, part of the mask

[`resolve_haplotypes()`](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L1163), per 40 Mb window, then
merged with `bedtools merge` ([precompute.py:128](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/precompute.py#L128)).

| Step | Where | Notes |
|---|---|---|
| A read counts at a position if its base differs from the assembly with base quality ≥ 40 | [utilities.py:1206](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L1206) | Substitutions only: reads with an InDel at the position are skipped (`query_position` None), so an unresolved InDel never triggers a region |
| A region is called where at least 2 such reads carry the same base | [utilities.py:1220](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L1220) | |
| The region spans the leftmost start to the rightmost end of those reads | [utilities.py:1222](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L1222) | |
| A region overlapping the current one extends it | [utilities.py:1228](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L1228) | **U1**: it sets the end to the new region's end, not the larger of the two, and never moves the start left, so regions can shrink |
| Not applied: coverage window, mapQ, alignment flags | [utilities.py:1190](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L1190) | The coverage test is commented out |

### The repeat tracts, the rest of the mask

[repetitive_dna.py](../workflows/repetitive_dna.py), per contig.

| Step | Where | Notes |
|---|---|---|
| Periods 2–6 only | [repetitive_dna.py:40](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/repetitive_dna.py#L40) | Homopolymers (k = 1) and periods of 7 or more are never masked (**M1**) |
| M0: mreps, exact, at least 3 copies, `-allowsmall` | [repetitive_dna.py:115](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/repetitive_dna.py#L115) | 25.9 million tracts in human/da1, 93% of them under 5 copies |
| M2: mreps, one substitution allowed between neighbouring copies (`-res 1`) | [repetitive_dna.py:132](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/repetitive_dna.py#L132) | 52.7 million tracts in human/da1, 86% under 5 copies |
| M1: TRF `2 7 7 80 10 50`, maximum period 6; output cut to k = 2–6 | [repetitive_dna.py:52](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/repetitive_dna.py#L52), [:168](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/repetitive_dna.py#L168) | Score ≥ 50, about 25 bp for a pure tract |
| M1 windows: a tract belongs to the window its start is in | [parse_trf.py:151](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/parse_trf.py#L151) | A tract longer than the 1 Mb overlap can be cut |
| Parsers drop rows without a full copy, or with missing fields | [parse_mreps_canonical.py:79](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/parse_mreps_canonical.py#L79), [parse_trf.py:130](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/parse_trf.py#L130) | |
| Overlap resolution, per model: smallest k first, a tract overlapping a kept one dropped whole; within one k, the longest kept | [resolve_repeat_overlaps.py:56](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/resolve_repeat_overlaps.py#L56) | **M1**: the mask is built from these resolved BEDs, so a dropped tract's bases stay unmasked unless another model covers them |
| Mask = M0 ∪ M1 ∪ M2 ∪ unresolved regions, merged | [indel_mutation_call.py:132–134](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/indel_mutation_call.py#L132-L134) | |

## 3. The scan

[`nonrepetitive_scan()`](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L548), one job per 10 Mb
window. The tests run in this order for every pileup position, then for every
read at it.

| # | Step | Where | Calls | Opps | Notes |
|---|---|---|---|---|---|
| 1 | A read enters the pileup: `stepper="nofilter"`, base quality ≥ 0, mapQ ≥ 0 | [utilities.py:623–625](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L623-L625) | ✓ | ✓ | No alignment-flag filter: secondary, supplementary, duplicate and QC-fail alignments enter (**F1**). htslib leaves out unmapped reads; pysam's `max_depth` of 8000 is never reached |
| 2 | Position in the window `[start, end)` | `truncate=True` | ✓ | ✓ | Each position belongs to one window |
| 3 | Coverage window: mean/2 ≤ reads at the position ≤ 2 × mean | [utilities.py:629](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L629) | ✓ | ✓ | Tested at the anchor for a call, at the k-mer's first base for an opportunity |
| 4 | Position outside the mask | [utilities.py:634](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L634) | ✓ | ✓ | Tested at the anchor only, so an InDel can remove masked bases (**L1**). pbmm2 places an InDel that can slide at its leftmost position (verified on real reads), so every InDel of whole units in a pure tract is anchored just before it |
| 5 | The read has a base at the position | [utilities.py:651](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L651) | ✓ | ✓ | Reads with a deletion across the position are skipped |
| 6 | Opportunity only: 1 ≤ read position ≤ read length − 2 | [utilities.py:660](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L660) | – | ✓ | The same test for all three sizes, so 2–3 bp k-mers at the last bases are cut short (**D3**) |
| 7 | Opportunity only: the k-mer plus one read base on each side changes base at least twice | [utilities.py:688–691](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L688-L691) | – | ✓ | Correct for deletions, too strict for 1 bp insertions (**D1**) |
| 8 | Call only: InDel of 1–3 bp | [utilities.py:701](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L701) | ✓ | – | |
| 9 | Call only, deletion: anchor + deleted bases + next base, from the assembly, change base at least twice | [utilities.py:731–736](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L731-L736) | ✓ | – | |
| 10 | Call only, insertion: anchor + inserted bases + next read base change base at least twice | [utilities.py:757–759](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L757-L759) | ✓ | – | |

The scan records these without filtering on them:
- **`unique`** ([utilities.py:644–648](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L644-L648)): `non-unique` when more than one read at the anchor has an InDel of any size or a deletion across it. All reads count, whatever their quality.
- **`kmerQ`**: calls at [lines 743 and 766](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L743), opportunities at [line 671](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L671). The windows differ (**D4**).
- **`ToreadEnd`** ([line 708](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L708)): the distance in assembly bases to the nearer end of the alignment.
- **The read's tags:** tissue, `np`, `ec`, and `sm`/`sx`/`sa`. The subread tags are read at the wrong base on reverse-strand reads (**S1**).
- **Diploid mapQ**, for opportunities only (**N1**).

## 4. Joining the windows

Nothing is filtered: calls are concatenated, opportunities summed, and the FASTQ
keeps each read once
([indel_mutation_call.py:179–203](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/indel_mutation_call.py#L179-L203)).

## 5. Remapping to the haploid assembly

Every read carrying a call, at any kmerQ, is aligned with
`pbmm2 --preset HIFI`. `samtools view` keeps every record it writes, sorted by
read name ([indel_mutation_call.py:223–224](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/indel_mutation_call.py#L223-L224)).
Reads that do not align are absent from the table.

## 6. The haploid check

[filtering_haploid_mappping.py](../scripts/filtering_haploid_mappping.py), calls
only.

| Step | Where | Notes |
|---|---|---|
| Join each call to its read's haploid records | [line 169](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/filtering_haploid_mappping.py#L169) | A read without a record gets no flag and is dropped at the flag test |
| Keep one row per read and call position | [line 170](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/filtering_haploid_mappping.py#L170) | The first record is kept, which is the primary alignment, as it sorts first |
| Keep flag 0 or 16 only | [lines 191–197](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/filtering_haploid_mappping.py#L191-L197) | Calls dropped here are silent: they are not in `filter_report.tsv` |
| Mirror the call's position when the strands differ | [lines 201–202](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/filtering_haploid_mappping.py#L201-L202) | Maps the anchor, not the event (**H2**) |
| Walk the haploid CIGAR to the first event of the call's type at or after its position; `haploid_check` is whether its size matches | [lines 68–148](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/filtering_haploid_mappping.py#L68) | Insertions don't advance the walk (**H1**); repeat InDels land at opposite ends on flipped reads (**H3**); a later, unrelated event can pass (**H4**) |

## 7. Final filtering

[indel_mutation_call.py:271–290](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/indel_mutation_call.py#L271-L290);
every step is reported in `filter_report.tsv`.

| Step | Where | Calls | Opps |
|---|---|---|---|
| `haploid_check` is True | [filtering_small_indels.py:127](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/filtering_small_indels.py#L127) | ✓ | – |
| kmerQ ≥ 93 | [filtering_small_indels.py:149](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/filtering_small_indels.py#L149), [filtering_opportunities.py:129](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/filtering_opportunities.py#L129) | ✓ | ✓ |

The scripts have options the workflow doesn't use:
- for calls: soft clip, haploid distance, InDel type, tissue, unique only and
  read-end distance;
- for opportunities: mapQ, tissue, context and count.

## 8. The notebook

[nonrepetitive_indels.ipynb](../plotting/nonrepetitive_indels.ipynb), by
section.

| Where | Step | Calls | Opps |
|---|---|---|---|
| Setup and loading | Samples whose `indels.tsv` or `opportunities.tsv` has no rows drop out of counts or rates; `species_keep` is NULL | ✓ | ✓ |
| Loading (`indels`), so every section | `distance_to_indel_in_hap` ≤ `max_hap_distance` (0) or NA | ✓ | – |
| Loading (`indels`) | `ToreadEnd` ≥ `min_read_end_distance` (0, so no effect) | ✓ | – |
| Unique and non-unique, size, context, contigs, groups | The two classes are kept apart; the context plot shows unique calls only | ✓ | – |
| *Where the calls are* | Contigs with an excess are flagged, not removed; the plot shows contigs with at least one call | – | – |
| *Distance to the read end* | Distance bins with fewer than 10⁵ positions are not shown; the slide stops at 100 kb | – | – |
| *Rate per tissue, read ends removed* | Unique, `ToreadEnd` > 300, tissue not NA | ✓ | tissue not NA; × `usable_share` |
| *Rate per bin* (germline/somatic) | Read in the likelihood table, unique, `ToreadEnd` > 300 | ✓ | split by read or CpG share; × `usable_share` (**N2**) |
| *Sperm against age* | Whole sperm samples, or reads labelled sperm; unique; `ToreadEnd` > 300; samples without unique calls dropped (la4, tsaed) | ✓ | × `usable_share` |
| *The human sperm samples pooled* | Human sperm samples with calls; unique; `ToreadEnd` > 300; samples without unique calls dropped | ✓ | × `usable_share` |
| *Compared with the repeat models* | Sizes 2–3; unique; outside, motifs with k bases that are not homopolymers; inside, one-unit InDels of the tract's motif; motifs with fewer than 10⁵ opportunities not plotted | ✓ | ✓ |

`usable_share` is the share of read positions more than 300 bp from a read end,
estimated from the reads that carry calls.

## 9. What reaches the sperm plot

**A call is counted if it:**
1. is on a contig longer than 1 Mb;
2. has its anchor in a position whose coverage is within 0.5–2× the mean;
3. has its anchor outside the mask;
4. is a 1–3 bp InDel on a read with a base at the anchor;
5. passes the homopolymer rule;
6. is on a read whose primary haploid alignment exists (flag 0 or 16);
7. passes `haploid_check`;
8. has kmerQ ≥ 93;
9. is matched in the haploid alignment at distance 0;
10. is `unique`;
11. lies more than 300 bp from a read end.

**A read position is counted if it:**
1. is on a contig longer than 1 Mb;
2. is at a position whose coverage is within 0.5–2× the mean;
3. is outside the mask;
4. is not the first or last base of its read;
5. starts a k-mer that passes the homopolymer rule in the read;
6. has kmerQ = 93.

It is then scaled by `usable_share`.

**Nothing on the opportunity side matches** call steps 6, 7, 9 and 10, which
only remove calls, so the rate errs low on that account. The read-end cut
(step 11) is matched only approximately, through `usable_share`. Steps 2–5 are
matched, apart from the different bases tested (L1), the different homopolymer
rules (D1), the read-end positions (D3) and the different kmerQ windows (D4).
