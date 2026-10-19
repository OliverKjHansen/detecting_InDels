# How far filtering lowers the rate outside repeats

*2 October 2026. A follow-up to [nonrepetitive_indel_review.md](nonrepetitive_indel_review.md),
whose finding IDs (L1, E1, D1, N1, …) are used here. Data: the tables in
`results/` copied from the cluster on 30 September 2026. Every number comes
from [check_rate_floor.py](review_checks/check_rate_floor.py), except the
forward-strand estimates, which come from
[check_candidate_filters.py](review_checks/check_candidate_filters.py).*

*Corrected on 7 October 2026. The subread tags are read at the wrong base on
reverse-strand reads (S1 in the review), so the subread figures below describe
a filter that works on half the calls. [candidate_filters.md](candidate_filters.md)
estimates what it does when it works; the summary and *What it would take*
include that estimate.*

## Summary

The expected base rate of InDels outside repeats is about 1×10⁻⁹ per bp, all
sizes and both types together. Summed over the six classes (1, 2 and 3 bp,
insertions and deletions), the median Revio sample is at 1.6×10⁻⁷, and the
Sequel II samples are higher still.

- **Filtering on what the calls record gains about 2× as the code stands.**
  Keeping long reads, calls far from read ends and calls that every subread
  agrees with brings the summed rate to a median of 8.7×10⁻⁸, and to 3.8×10⁻⁸
  in the best sample. That is still 40–90× the target.
- **With the subread tags read correctly (S1), the same filters gain about 3×.**
  On forward-strand reads, where the tags are read at the right base, the
  summed rate falls to 5.3×10⁻⁸. Almost all of the extra gain is in the
  insertions, which come down to 0.75–2.3×10⁻⁹ per class.
- **What survives is mostly deletions that look like real sequence in the
  read.** These calls sit on long reads, far from the ends, and every subread
  agrees with them. For a deletion, though, the tags test a base beside the gap,
  not the deletion.
- **It is not mutation either.** Sperm samples still show 5–10×10⁻⁹ for 1 bp
  deletions alone. Germline and somatic reads in the same testis samples have
  the same rate (E1).
- **Four sources remain, and each needs a new test:** the mask-boundary leak
  (L1; 2 bp deletions are 46% of what is left), deletion artefacts that the
  subread tags cannot see, mismapped reads that carry real germline
  differences, and single-strand events merged into one consensus. See *What
  it would take*.

## Method

The opportunity table counts read positions by context, tissue, mapQ and
kmerQ. It cannot be split by read length, passes or subread support, so there
is no direct rate for, say, reads shorter than 4 kb. Instead, unique calls are
compared with **control calls** that must be real sequence in the read: the
same InDel at the same site in at least three reads. These are germline
differences or assembly errors. Real sequence is present whatever the read's
quality, so the controls are spread over reads the way Q93 positions are, and
they stand in for the opportunities.

- **Relative rate of a group of reads.** Within a sample, the share of its
  unique calls in the group, divided by the share of its controls in the group.
  1 is the sample's average; the tables give the median over samples.
- **Rate after a filter.** The calls the filter keeps, over the opportunities,
  divided by the share of controls the filter keeps. That share estimates the
  share of real mutations the filter keeps.

Both sets are matched to the haploid alignment at distance 0 and lie at least
100 bp from a read end. The controls fit the opportunities only roughly: they
sit at variant sites, which may differ from typical positions in ways this
analysis cannot see.

On real reads the fit is close. In one chk window, 28.0% of Q93 read positions
pass the subread filter, read correctly, against about 32% of the controls
(review, *Tested on real reads*). Rates corrected with the controls are
therefore slightly low, by about an eighth.

## What the read-level features show

Rate of unique calls relative to the sample's average, median over samples:

| read length | < 4 kb | 4–8 kb | 8–12 kb | 12–16 kb | 16–20 kb | ≥ 20 kb |
|---|---|---|---|---|---|---|
| relative rate | 3.8 | 2.1 | 1.2 | 0.82 | 0.71 | 0.73 |

| subread support at the call (`sm`/`sa`) | < 0.8 | 0.8–0.9 | 0.9–0.95 | 0.95–<1 | 1 |
|---|---|---|---|---|---|
| relative rate | 6.4 | 1.2 | 0.80 | 1.3 | 0.68 |

| subreads disagreeing (`sx`) | 0 | 1 | 2 | ≥ 3 |
|---|---|---|---|---|
| relative rate | 0.78 | 1.1 | 1.9 | 6.8 |

- **Short reads carry an artefact.** A real mutation does not depend on read
  length, yet reads under 4 kb carry 3.8× the sample's average rate. The rate
  levels off from 12 kb.
- **Passes (`np`) follow read length.** The rate is lowest at 10–14 passes
  (0.66) and rises to 2.9 at 30 or more. Many passes means a short insert, so
  this mostly repeats the read-length effect.
- **Subread disagreement marks artefacts, more than these figures show.** Calls
  with under 80% subread support, or with three or more disagreeing subreads,
  run at 6–7× the average. On reverse-strand reads, though, the tags describe
  the wrong base (S1). On forward-strand reads alone, under 80% support runs at
  24× and covers a third of the unique calls.
- **The tags test insertions better than deletions.** `sm` and `sx` are read at
  the base after the anchor. For an insertion that is the first inserted base;
  for a deletion it is an ordinary base beside the gap. Full support separates
  unique insertions from the controls (17% against 44% for 1 bp) more sharply
  than deletions (30% against 46%).
- **SMRT Cells differ by up to about 3×.** Within a sample, single cells run
  from 0.45 to 2.9 times the sample's average. In bonobo/matos,
  marmoset/lambiek and rhesus R21008 and R22013, one cell runs at 1.9–2.9×.
- **Sequel II runs 2.6× higher than Revio.** The median 1 bp rate is
  1.5×10⁻⁷ in the 6 Sequel II samples with calls (carl, stephan, mbewe,
  samson, da1, ob2v), against 5.6×10⁻⁸ in the 25 Revio samples. The Sequel II
  BAMs carry none of the subread tags, and those of mbewe and samson lack `np`
  too.

## The cleanest subset

Revio samples only, because the filters use the subread tags. A call is kept
if its read is at least 12 kb long, it lies at least 1 kb from a read end, no
subread disagrees at the tested base, and every subread covering that base
matches it. These filters keep a median of 29% of the controls. The share is
lowest in aa_sperm and ob007, whose reads are short, so their filtered rates
are uncertain.

The subread test here uses the tags as the code reads them, which is right on
forward-strand reads only (S1). Repeated on forward-strand reads, the same
filters bring the summed rate to 5.3×10⁻⁸ rather than 8.7×10⁻⁸, with the
insertions at 0.75–2.3×10⁻⁹ per class ([candidate_filters.md](candidate_filters.md)).

| class | before | after |
|---|---|---|
| 1 bp deletion | 2.1×10⁻⁸ | 9.8×10⁻⁹ |
| 1 bp insertion | 3.2×10⁻⁸ | 5.5×10⁻⁹ |
| 2 bp deletion | 5.7×10⁻⁸ | 3.8×10⁻⁸ |
| 2 bp insertion | 1.9×10⁻⁸ | 8.1×10⁻⁹ |
| 3 bp deletion | 1.5×10⁻⁸ | 8.4×10⁻⁹ |
| 3 bp insertion | 6.6×10⁻⁹ | 3.2×10⁻⁹ |

These are medians over the 25 Revio samples. Summed per sample, the median
falls from 1.6×10⁻⁷ to 8.7×10⁻⁸ (1.8×). The lowest sample,
long_tailed_macaque/J12013, reaches 3.8×10⁻⁸.

- **1 bp insertions fall most (6×).** That fits tags that test the inserted
  base itself. With their own denominator (D1) they would fall about 2× more.
- **2 bp deletions fall least, and are 46% of what is left.** They are slippage
  events from masked tracts (L1): real sequence that no quality filter removes.
- **The four sperm samples with enough calls (chk, la, ph, tr) are left with
  4.8–9.8×10⁻⁹ for 1 bp deletions alone.** That is above the whole germline
  InDel rate expected outside repeats.

## What it would take

Ordered by how much each is likely to remove.

1. **Fix S1 and filter on subread support.** This is the cheapest step: index
   the per-base tags at `query_length − 1 − i` on reverse-strand alignments.
   On the forward-strand estimate it cuts the insertions 5–8×.
2. **Fix L1.** 2 bp deletions are about half of what is left after the
   filters, and they are not read artefacts. Testing the anchor and the k bases
   after it against the mask, for calls and opportunities alike, takes them
   out of the outside rate (review, L1).
3. **Test for mismapped reads.** A read from the other haplotype, or from a
   paralog the assembly collapsed, carries real germline differences, so it
   passes every filter above. The haploid check cannot catch it reliably:
   wherever the haploid assembly carries the other allele, the read shows the
   same InDel there too. Two steps, the cheap one first:
   - Write the read's diploid mapQ and the number of mismatches near the call
     into the call table (the opportunities already carry mapQ; N1). Then see
     whether the rate falls at mapQ 60 with few mismatches.
   - Align each call-carrying read to hap1 and hap2 separately, and keep a call
     only if the InDel appears in both alignments. A real mutation is in the
     read and shows against both. A germline difference disappears against the
     haplotype the read came from.
4. **Require both strands.** A HiFi read is one consensus of both strands of
   the molecule, so a change on one strand, such as damage or a mismatch, can
   be called with high quality. Published single-molecule HiFi methods,
   HiDEF-seq among them, build a consensus per strand and require both to
   carry the change. That needs the subreads, so check whether they were kept;
   that is more likely for the Sequel II runs. It is also the most direct test
   for the deletions, which the subread tags cannot check.
5. **Give the read filters real denominators.** Add read-length and
   read-end-distance bins to the opportunity table in the scan, so that these
   filters no longer rest on the control yardstick. Subread support can get one
   too: record the same test for the k-mer at each opportunity.
6. **Drop or separate the Sequel II samples,** and SMRT Cells running at 1.5×
   their sample's average or more.

## A trap

Filtering calls on read length, passes or subread support lowers the rate
whether or not it removes artefacts, as long as the opportunities stay
unfiltered. 1×10⁻⁹ can be reached that way, and it would be wrong. Every call
filter needs either the same filter on the opportunities or a correction for
the real calls it removes, as in *The cleanest subset*.

## Reproducing the numbers

From the repository root, with `results/` copied from the cluster:

```bash
.venv/bin/python documents/review_checks/check_rate_floor.py
.venv/bin/python documents/review_checks/check_candidate_filters.py
```
