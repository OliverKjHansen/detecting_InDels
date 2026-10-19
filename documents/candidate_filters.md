# What else could be filtered on

*7 October 2026. A companion to [nonrepetitive_indel_review.md](nonrepetitive_indel_review.md)
(finding IDs) and [filtering_steps.md](filtering_steps.md) (what is filtered
now). Every number comes from [check_candidate_filters.py](review_checks/check_candidate_filters.py)
or [check_rate_floor.py](review_checks/check_rate_floor.py), run on the tables
in `results/` copied on 30 September 2026. The real-read figures come from
[check_mirror_window.py](review_checks/check_mirror_window.py), on one human/chk
window. Links to lines of code point at commit 542fa20, the code these
numbers describe.*

## Summary

The calls carry several properties that mark artefacts. The strongest, subread
support, works on only half the calls today, because the tags are read at the
wrong base on reverse-strand reads (S1).

- **With S1 fixed, subread support cleans up the insertions.** Estimated from
  forward-strand reads, where the tags are read correctly, the summed rate
  falls from 1.6×10⁻⁷ to 6.9×10⁻⁸. Adding reads ≥ 12 kb and calls ≥ 1 kb from
  a read end brings it to 5.3×10⁻⁸. The insertion classes then sum to about
  5×10⁻⁹.
- **What is left is mostly deletions.** Under the strictest filter, 1, 2 and
  3 bp deletions stay at 8.2×10⁻⁹, 2.5×10⁻⁸ and 7.3×10⁻⁹. The subread tags are
  read at a base beside the gap, so they barely test a deletion. 2–3 bp
  deletions carry the repeat leak (L1).
- **So the next gains need new information from the scan:** the L1 mask fix, a
  deletion-specific test, and a test for mismapped reads.

## How a candidate is judged

The opportunity table cannot be split by read properties, so each candidate is
measured against a yardstick. That yardstick is control calls that must be real
sequence in the read: the same InDel at the same site in at least three reads.
A value of 1 means "the sample's average rate"; 3 means three times as many
unique calls per real-sequence-equivalent read position. Medians over samples
are reported. [rate_floor.md](rate_floor.md) describes the method and its limits.

On real reads the yardstick holds approximately. In one chk window, 28.0% of
Q93 read positions pass the subread filter (read correctly), against about 32%
of the controls (review, *Tested on real reads*). So the corrected rates below
are slightly optimistic, by about an eighth.

**The trap.** A filter applied to the calls but not to the opportunities always
lowers the rate, artefacts or not. Each candidate below says how its
denominator would be matched:
- **exact**: the opportunity table already has the field;
- **yardstick**: divide by the share of controls the filter keeps;
- **scan**: add the field to the opportunity key.

## Candidates, by expected gain

| # | Filter | Evidence | Unique calls affected | Denominator |
|---|---|---|---|---|
| 1 | **Subread support, after fixing S1**: no disagreeing subread at the anchor or the base after it, and every covering subread agreeing | Forward strand: under 80% support at **24×**, full support at 0.49×. Insertions fall 5–8× (table below). On real reads, 2 of 19 unique kmerQ ≥ 93 calls have full support, against 43–45% of the controls | about half, once both strands are read correctly | yardstick, or a support flag at the opportunities (scan) |
| 2 | **Mask the bases an InDel touches** (L1): the anchor and the k bases after it | 2 bp deletions are the largest class left (46% after the clean filters, 2.5×10⁻⁸ under the strictest); one-unit InDels inside M0 k = 2 tracts at 0.46× the outside rate. On real reads, 93% of kmerQ ≥ 93 2 bp calls touch a masked base | most 2–3 bp calls | the same mask test on both sides (scan) |
| 3 | **Read length ≥ 12 kb** | < 4 kb at 3.8×, 4–8 kb at 2.1×, 8–12 kb at 1.2×; ≥ 12 kb at 0.71–0.82× | 41% (on reads < 12 kb) | read-length bin (scan), or yardstick |
| 4 | **≥ 1 kb from a read end** | Up to 28× within 25 bp of an end (notebook, *Distance to the read end*) | 31% | end-distance bin (scan), `usable_share`, or yardstick |
| 5 | **Drop contigs with an excess of unique calls** (one-sided Poisson, p < 0.001 / contigs) | A median of 10 contigs per sample hold 26% of the unique calls on 15% of the read positions; dropping them lowers the rate to 0.86× | 26% | **exact**: opportunities are counted per contig |
| 6 | **Haploid mapQ 60** | mapQ < 60 at 2.4× | 4% | yardstick |
| 7 | **No soft clip in the haploid alignment** | soft-clipped at 2.2× | 4% | yardstick |
| 8 | **Drop SMRT Cells at 1.5× their sample or more** | cells run from 0.45× to 2.9× of their sample; one cell at 1.9–2.9× in marmoset/lambiek, rhesus R21008, bonobo/matos and rhesus R22013 | depends on the sample | yardstick, or a cell field (scan) |
| 9 | **Drop clustered calls** (another unique call within 1 kb) | 2.2% have one, against 0.13% expected at random; half of those within 100 bp | 2% | mask those windows on both sides (exact if done in the scan) |
| 10 | **Analyse Sequel II separately** | the 1 bp rate is 2.6× the Revio median; no subread tags; mbewe and samson have no `np` either | 6 samples | – |

Checked and not useful:
- **Passes (`np`, `ec`).** These rise at 30 or more (2.8–2.9×), but many passes
  means a short insert, so they mirror read length. Filter on length instead.
- **Diploid strand, haplotype (h1/h2) and H/P/N operations.** All flat (0.93–1.1×).
- **Supplementary alignments.** No read carries calls at two read lengths, so
  no calls come from hard-clipped alignments. On real reads, 0.1% of calls sit
  on soft-clipped supplementary alignments.

### What a working subread filter does

Forward-strand calls only, where the tags are read at the right base, rates
corrected per class with forward-strand controls. These are medians over the 25
Revio samples.

| class | no filter | + strand, subreads agree | and ≥ 12 kb, ≥ 1 kb from an end |
|---|---|---|---|
| 1 bp deletion | 2.1×10⁻⁸ | 1.2×10⁻⁸ | 8.2×10⁻⁹ |
| 1 bp insertion | 3.2×10⁻⁸ | 3.9×10⁻⁹ | 2.0×10⁻⁹ |
| 2 bp deletion | 5.7×10⁻⁸ | 3.2×10⁻⁸ | 2.5×10⁻⁸ |
| 2 bp insertion | 1.8×10⁻⁸ | 3.7×10⁻⁹ | 2.3×10⁻⁹ |
| 3 bp deletion | 1.5×10⁻⁸ | 1.1×10⁻⁸ | 7.3×10⁻⁹ |
| 3 bp insertion | 6.6×10⁻⁹ | 1.1×10⁻⁹ | 7.5×10⁻¹⁰ |
| **sum** (median of per-sample sums) | **1.6×10⁻⁷** | **6.9×10⁻⁸** | **5.3×10⁻⁸** |

The filters keep 15% and 10% of the controls, so these estimates rest on few
calls. The lowest samples under the strictest filter sum to 2.3–2.9×10⁻⁸
(bonobo/matos, chimp/yoran, rhesus R15079).

## Fixing S1

The subread tags `sm`, `sx` and `sa` are per-base arrays. PacBio stores
per-base tags in the read's own orientation, and pbmm2 does not reverse them for
a reverse-strand alignment. The scan, however, indexes them with the position in
the aligned, reverse-complemented read
([utilities.py:716–722](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L716-L722), and
`subread_array()` for `sa`). On a reverse-strand read it therefore looks at a
base at the mirrored position, far from the call. The bug changes no call and no
rate, only these columns; the review's S1 section explains what it does and
does not affect.

The evidence, genome-wide:
- **The real-sequence controls look the same on both strands:** 45% and 43%
  have every subread agreeing.
- **The unique calls do not.** On forward-strand reads 21% of them have full
  support and 37% have two or more disagreeing subreads; on reverse-strand
  reads, 34% and 16%.
- **Calls with under 80% support** run at 24× the average on forward-strand
  reads, against 1.3× on reverse-strand ones.

To fix it:
- On a reverse-strand alignment, index each per-base tag at
  `query_length − 1 − i`. For `sa`, decode the run-length encoding first, then
  index the same way.
- For an insertion, test every inserted base, not only the first.
- Record the same support test for the k-mer at each opportunity, so the filter
  gets an exact denominator.

Real reads settle the orientation. A low base quality goes with disagreeing
subreads, and pysam returns base qualities in the aligned order. In a chk window
mirrored from the cluster, subread support correlates with base quality at a
median of 0.32 per read on forward-strand reads. On reverse-strand reads it is
0.000 with the tags indexed the scan's way, and 0.32 once they are read in the
read's own order. At the window's call sites, reverse-strand calls look fully
supported in 70.0% of cases as recorded, against 37.7% once read correctly,
close to forward-strand calls (39.3%).
[check_mirror_window.py](review_checks/check_mirror_window.py) reruns these
tests on any mirrored window.

## Ideas that need new information from the scan

In the order they are likely to matter:

1. **Mismapped reads.** A read from the other haplotype, or from a paralog the
   assembly collapsed, carries real germline differences, so it passes every
   filter above. Two steps:
   - Write the read's diploid mapQ and its mismatch count near the call into
     each call; the opportunities already carry mapQ (N1).
   - Align each call-carrying read to hap1 and hap2 separately, and keep the
     call only if it appears against both.
2. **A test made for deletions.** The tags at the base after the gap don't say
   whether the subreads agree on the deletion. That needs the subreads
   themselves, realigned to the read, or a separate consensus per strand
   (by-strand CCS) with both strands required to carry the change. Check
   whether the subreads were kept.
3. **Sequence context, on both sides.** The rules would cover:
   - distance to the nearest homopolymer run (homopolymers are not masked, M1);
   - local sequence complexity. The collaborators' SNV notebook already filters
     on entropy and low complexity (`entropy`, `low_complexity` in
     [test_mutation_rate.ipynb](../test/test_mutation_rate.ipynb)).

   Apply the same rules to the k-mer at each opportunity.
4. **Distance to the nearest other InDel or mismatch in the same read.**
   Alignment artefacts come in clusters. A helper that measured this from the
   CIGAR, `closest_indel()`, was never called by the scan and has since been
   deleted. For the denominator, record the same distance bin at each
   opportunity.
5. **A wider all-Q93 window**, say five bases on each side, for calls and
   opportunities alike. The current window is three or four bases, and not the
   same ones on both sides (D4).
6. **Closing the mask holes.**
   - Build the mask from every tract before overlap resolution (M1).
   - Fix the region bounds of the unresolved-haplotype scan (U1).
   - Consider masking long homopolymers, periods of 7 or more, satellites and
     segmental duplications.
7. **Primary alignments only** (F1), for calls and opportunities alike. Little
   gain is expected, since no supplementary calls were seen, but it removes a
   source of doubt.
