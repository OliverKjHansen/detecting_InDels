# Review: InDel rates outside repeats

*2 October 2026. Covers workflow 3
([indel_mutation_call.py](../workflows/indel_mutation_call.py)), the haploid
check ([filtering_haploid_mappping.py](../scripts/filtering_haploid_mappping.py))
and [nonrepetitive_indels.ipynb](../plotting/nonrepetitive_indels.ipynb). Code
at commit 542fa20, which the links to lines of code point at (the notebook was
not yet committed). Data: the tables in
`results/` copied from the cluster on 30 September 2026, 35 samples. Every
number here can be reproduced with the scripts in [review_checks/](review_checks/).*

*Fact-checked on 7 October 2026 against the code and the data. That pass
corrected a link and a few statements, and added findings S1, U1, M1, F1 and
D4. Every filtering step is listed in [filtering_steps.md](filtering_steps.md),
and what else could be filtered on is in
[candidate_filters.md](candidate_filters.md).*

*Tested on real reads on 7 October 2026, in one window of human/chk copied from
the cluster. The findings that could be tested held up; the details are under
[Tested on real reads](#tested-on-real-reads).*

## Summary

The 1–3 bp InDel rates outside repeats came out 10–50× higher than expected.
Across samples the 1 bp rate runs from 2.1e-8 to 3.5e-7 per read position
(median 6.1e-8), and in 27 of 31 samples the 2 bp rate is higher than the 1 bp
rate.

The opportunity table, the denominator, is not the cause. It holds about the
number of read positions it should, and the largest problem in its rules (D1)
makes 1 bp insertion rates about 2× too high. The excess is in the calls:

- **2–3 bp:** InDels inside repeat tracts are reported at the base before the
  tract, which is outside the mask. Workflow 3 counts them, while its
  opportunities leave the tract out (L1).
- **All sizes:** a floor of single-read sequencing errors passes every filter,
  kmerQ ≥ 93 included (E1).

The best filter against that floor, subread support, is blind on half the
calls: the tags are read at the wrong base on reverse-strand reads (S1). That
bug leaves the set of calls and the rates untouched; it corrupts five columns
that nothing filters on.

The haploid check has bugs of its own (H1–H4). They remove real calls. On real
reads, about three quarters of all 2–3 bp calls are lost on every contig, and
most of the rest on contigs whose haploid counterpart runs the other way.
Fixing them will raise the 2–3 bp rates several-fold, not lower them.

How far filtering on read-level features lowers the rate, and what it would
take to reach the expected 1×10⁻⁹, is in [rate_floor.md](rate_floor.md).

| ID | Finding | Where | Effect on the rate | On real reads |
|---|---|---|---|---|
| L1 | Repeat InDels counted at the base before a masked tract | [utilities.py:634](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L634) | too high, 2–3 bp | confirmed: 93% of the window's kmerQ ≥ 93 2 bp calls |
| E1 | Sequencing errors still reach kmerQ 93 | the method | too high, all sizes | not testable in one window |
| S1 | The subread tags are read at the wrong base on reverse-strand reads | [utilities.py:716](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L716-L722) | none directly; the best artefact filter fails on half the calls | confirmed |
| D1 | 1 bp insertions get the 1 bp deletion denominator | [utilities.py:669](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L669-L693) | 1 bp insertions about 2× too high | confirmed: 2.0× |
| H1 | Insertions do not advance the haploid CIGAR walk | [filtering_haploid_mappping.py:106](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/filtering_haploid_mappping.py#L106) | too low: about three quarters of 2–3 bp calls lost | confirmed |
| H2 | A strand flip maps the anchor, not the event | [filtering_haploid_mappping.py:201](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/filtering_haploid_mappping.py#L201-L202) | too low | confirmed across samples |
| H3 | On flipped reads, InDels in repeats land at opposite ends | [filtering_haploid_mappping.py:87](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/filtering_haploid_mappping.py#L87) | too low, mostly 2–3 bp | confirmed across samples |
| H4 | `haploid_check` passes on a later, unrelated event; the distance mixes two measures | [filtering_haploid_mappping.py:87](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/filtering_haploid_mappping.py#L87) | both ways | confirmed: 0.2% of passing calls once H1–H3 are fixed |
| D2 | Samples whose base qualities never reach 93 drop out silently | [indel_mutation_call.py:57](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/indel_mutation_call.py#L57) | two samples lost | not testable with chk |
| D3 | k-mers cut short at read ends are counted under shorter contexts | [utilities.py:660](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L660) | negligible | confirmed: 1×10⁻⁴ of opportunities |
| D4 | Calls and opportunities test kmerQ over different bases | [utilities.py:671](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L671) | small | confirmed: 1 bp deletions low by 12% |
| U1 | Unresolved-haplotype regions can shrink as they are built | [utilities.py:1228](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L1228) | too high, slightly | confirmed: 1.4% of unresolved bp lost |
| M1 | Mask holes: tracts dropped by overlap resolution; repeat classes never annotated | [indel_mutation_call.py:132](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/indel_mutation_call.py#L132-L134) | too high, slightly | 0.7% of exact-repeat bases unmasked |
| F1 | No alignment-flag filter | [utilities.py:623](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L623) | probably small | 0.1% of calls on supplementary alignments |
| N1–N5 | Smaller problems in the calls, the notebook and the README | see *Other issues* | small | not testable in one window |

## How the numbers were made

Unless stated otherwise, calls are `unique`, matched to the haploid alignment
at distance 0, and at least 100 bp from a read end. The first two are the
notebook's defaults. The read-end cut is the one workflow 4 applies, so rates
inside and outside repeats are filtered alike. A rate is calls divided by the
Q93 opportunities of the same size.

A read is **flipped** when its alignment to the haploid assembly is on the other
strand from its alignment to the diploid assembly. 17% of the final calls are on
flipped reads.

## The denominator is about the right size

For the 12 samples whose read count is known (from `classification_reads.tsv`),
the Q93 opportunities of size 3 are 14–55% of the bases sequenced, or 18–49%
when the bases are estimated from CpGs rather than read length. Size 3 is the
fairest comparison, because the homopolymer rule removes only about 5% of
positions at that size. Those shares are what the scan's own rules leave: the
mask (in the three samples checked, M2 tracts alone hold 12–15% as many read
bases as the outside), the coverage window, the homopolymer rule, and four
bases that must all sit at Q93. A 10–50× undercount is impossible: even the
lowest share, multiplied by 10, would exceed the bases sequenced.

Apart from L1, D1, D3 and D4, calls and opportunities pass the same tests in
the same pass over the BAM. [filtering_steps.md](filtering_steps.md) lists
every test and which side it applies to. The definitive check needs the
unfiltered opportunity table on the cluster (check 2 below).

## Findings that make the rate too high

### L1. Repeat InDels counted at the base before a masked tract

The mask is tested at the anchor only
([utilities.py:634](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L634)). An InDel of whole repeat
units in a pure tract can be placed at any copy, and pbmm2 places it at the
leftmost copy. On real reads none of about 48,000 InDels in 3,000 alignments
could still slide left (*Tested on real reads*). So the anchor is the base before
the tract. That base is not masked, but the bases the InDel removes or extends
are. The call is counted, and no opportunity matches it, because the k-mer it
would be counted under starts inside the tract. The README notes the missing
padding (*The anchor base*); the data show that the effect is large.

The mask has a great many such edges. human/da1 has 25.9 million M0 tracts and
52.7 million M2 tracts, most of them a few copies of a dinucleotide, so the
sequence outside the mask comes in short pieces.

- **Exact tracts mutate less than unique sequence.** One-unit InDels inside M0
  tracts with k = 2 occur at 0.46× the rate of 2 bp InDels outside repeats
  (median over 30 samples, and below 1 in 27 of them). Exact dinucleotide
  repeats should mutate far faster than unique sequence. In M1 tracts, where an
  impurity stops the leftward shift, the rate is 66× the outside rate. The
  events missing from M0 are in workflow 3's calls.
- **The outside spectrum looks like slippage.** Pooled over samples, 2 bp
  deletions outside repeats run at 0.9–1.1e-7 in the contexts that form
  dinucleotide repeats (GA, TA, AT, AG, AC, CA), against 3.0e-8 for GC and
  1.3e-8 for CC.
- **Most 2–3 bp calls sit where the gap can slide.** See H3.
- **On real reads the leak shows directly.** In the chk window, by chance 0.92%
  of unmasked positions sit right before a mask interval. The first base an
  InDel touches is masked for 17.4% of unique 2 bp calls (19× chance), 93% of
  non-unique 2 bp calls, and 93% of the window's 128 kmerQ ≥ 93 2 bp calls.

How much it explains: exact k = 2 tracts hold about 4% as many read bases as
the outside. If they mutated at 3%, 10% or 30% of the M1 per-base rate, the
events missing from them would be 6%, 27% or 86% of the 2 bp calls outside
(median over samples). That is M0 with k = 2 alone. Pure stretches at the left
end of M1 and M2 tracts, k = 3 tracts (for 3 bp calls), and short repeats below
every model's threshold, such as (CA)2, add to it.

L1 is the largest lever left. After the read-level filters in
[rate_floor.md](rate_floor.md), 2 bp deletions are still 46% of what remains,
because slippage events are real sequence that no quality filter removes.

**Suggested fix.** Test the same bases on both sides. A call passes the mask if
its anchor and the k bases after it are unmasked. An opportunity for the k-mer
at read position p (anchor p − 1) passes if p − 1 to p + k − 1 are unmasked. In
workflow 4, assign an InDel anchored at the base before a tract to that tract,
so the events move to the repeat rates rather than vanish.

### E1. Sequencing errors still reach kmerQ 93

kmerQ ≥ 93 is the strongest filter. It keeps 0.07–2.7% of the haploid-checked
calls (median 0.4%). What it keeps still sums to 1.6×10⁻⁷ per read position
over the six classes, from 6.6×10⁻⁹ for 3 bp insertions to 5.7×10⁻⁸ for 2 bp
deletions (medians over the Revio samples). Four signs point to sequencing
artefacts rather than mutations:

- **Germline and somatic reads have the same rate.** The germline/somatic rate
  ratio is 0.83–1.23 in all 12 classified samples (notebook, *Rate per bin*).
- **Tissues barely differ.** In chimp/yoran, nine of ten tissues lie between
  1.8e-8 and 2.9e-8 for 1 bp InDels (notebook, *Rate per tissue, read ends
  removed*).
- **Samples differ more than biology allows.** The 1 bp rate varies 17× across
  samples and up to 6× within a species, with gorilla mbewe and samson, chimp
  carl and stephan, and human ob007, ob2v and aa_sperm at the top. Five of
  those seven are Sequel II runs.
- **Unique calls crowd the read ends.** Up to 28× as many as expected fall
  within 25 bp of a read end (notebook, *Distance to the read end*). The core
  sections of the notebook (`size_rates`, `context_rates`) keep these calls,
  because `min_read_end_distance <- 0`.

This is a limit of the method rather than a bug: one read's base qualities
cannot push the error rate below this floor. [rate_floor.md](rate_floor.md)
measures how far the other read-level features go. Subread support (`sm`/`sx`,
recorded for the 25 Revio samples but not for the Sequel II ones), read length
and distance from the read ends together lower the summed rate only about 2× as
the code stands. With S1 fixed, the estimate is 3×, to 5.3×10⁻⁸, and the
insertions come down to 0.75–2.3×10⁻⁹ per class. The deletions stay high. Both
documents list what would go further: a test made for deletions, a test for
mismapped reads, and a consensus per strand with both strands required to carry
the change.

### S1. The subread tags are read at the wrong base on reverse-strand reads

PacBio writes three per-base arrays into each HiFi read:
- `sm`, the subreads that match the consensus at each base;
- `sx`, the subreads that mismatch it;
- `sa`, the subreads that cover it, run-length encoded.

They are stored in the order the read was sequenced. When a read aligns to the
reverse strand, pbmm2 reverse-complements the sequence and reverses the base
qualities, but leaves these arrays as they are. The scan finds a call's position
in the aligned read and uses that index into the arrays
([utilities.py:716–722](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L716-L722), and
`subread_array()` for `sa`). On a forward-strand read that is the right base. On
a reverse-strand read it is base `read_length − 1 − i`: the mirror image of the
call's position, usually thousands of bases away. The repeat scan does the same
([utilities.py:1069–1073](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L1069-L1073)).

An example from the chk window. A 1 bp insertion of A at h1tg000001l:1,502,222
sits on a reverse-strand read of 5,983 bases, anchored at aligned base 1,408.
- **What the scan recorded:** the arrays at index 1,409, 13 of 13 subreads
  matching and none mismatching, so the call looks perfectly supported.
- **The inserted base:** sits at index 4,573. There, 8 of 15 subreads match and
  5 disagree.

**How it was established.**
- **Genome-wide, from `results/`.** The real-sequence controls (calls whose
  InDel is in at least three reads) have every subread agreeing in 45% of
  forward-strand and 43% of reverse-strand cases. The unique calls do not
  behave alike: 21% have full support on forward-strand reads but 34% on
  reverse-strand ones. Calls with under 80% support run at 24× the sample's
  average on forward-strand reads, against 1.3× on reverse-strand ones.
- **On real reads.** In the chk window's reads, subread support follows base
  quality: quality is low where subreads disagree. On forward-strand reads the
  median correlation is 0.32. On reverse-strand reads it is 0.000 with the
  arrays indexed the scan's way, and 0.32 once they are read in the read's own
  order. Recomputing the scan's values reproduces the cluster's recorded ones
  for all 10,487 calls, so the indexing above is what the pipeline does.

**What S1 changes in your call set.**
- **Five columns on about half the rows.** `sm`, `sx`, `sa`, `sm_pre` and
  `sx_pre` are wrong on every call from a reverse-strand read: 49% of the
  window's calls. In the window the recorded value is wrong in 76% (`sm`), 44%
  (`sx`), 51% (`sa`), 69% (`sm_pre`) and 41% (`sx_pre`) of those rows. Where it
  agrees, that is chance: `sx` is 0 at most bases.
- **Where.** Every call table of workflows 3 and 4: `calls/`, `merged/`,
  `results/…/indels.tsv` and `indels_classified.tsv`.
- **Plot 8 of `plots.pdf`**, the histogram of `sa`.
- **Subread-based numbers in [rate_floor.md](rate_floor.md)**, which are marked
  there.

**What S1 does not change.** No filter and no notebook reads these columns. So
these are all right:
- which calls are made, and their positions, sequences and sizes;
- kmerQ (pysam returns base qualities in the aligned order, which is right);
- `unique`, `ToreadEnd`, `np`, `ec`, tissue and the haploid columns;
- the opportunities and every rate.

**Why it matters anyway.** Read correctly, the arrays show that most unique
calls are poorly supported by their subreads. Of the window's 19 unique kmerQ ≥
93 calls, 2 have every subread agreeing at the base after the anchor, against
43–45% of the real-sequence controls genome-wide. As recorded, 6 of the 10 on
reverse-strand reads looked fully supported; in truth 1 is. Missing support
does not prove a call false, but it is the best artefact signal the reads
carry, and S1 hid it on half the calls. Fixing it removes no call by itself,
but it makes a subread filter usable
([candidate_filters.md](candidate_filters.md)).

**To get the right values** for an existing call: on a reverse-strand read
(`diploid_direction` "-"), read the arrays at `read_length − 1 − i` instead of
`i`. Here `i` is `start` for `sm` and `sx` (the base after the anchor), and
`start − 1` for `sa`, `sm_pre` and `sx_pre` (the anchor). Decode `sa` before
indexing it.

**Suggested fix.** In both scans, index the arrays that way on reverse-strand
alignments. Test every inserted base, not only the first, and record the same
test at the opportunities so the filter gets an exact denominator.

### D1. 1 bp insertions get the 1 bp deletion denominator

A 1 bp opportunity is a read base that differs from both of its neighbours
([utilities.py:669–693](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L669-L693)), which is where a
1 bp deletion passes the homopolymer rule. A 1 bp insertion passes the rule
when the inserted base differs from the anchor and from the next base
([utilities.py:757](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L757)), which is possible at nearly
every gap between two bases. On real reads only 50.6% of bases differ from
both neighbours (the 1 bp opportunities are 60–70% of the size-3 ones, which
also need one more base at Q93), so 1 bp insertion rates are about 2× too high.
That explains much of why 1 bp insertions outnumber deletions. Per context it is worse: a C insertion
is divided by the read positions holding a C, which have nothing to do with
where a C can be inserted.

**Suggested fix.** Count insertion opportunities separately, one per gap, with
the anchor and the next base as the kmerQ window. Drop or redefine the
per-context insertion rates.

## Findings that make the rate too low: the haploid check

The haploid check confirms a call by finding the same InDel in the read's
alignment to the haploid assembly.
[find_indel_position_in_read()](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/filtering_haploid_mappping.py#L68)
walks the haploid CIGAR from the start of the read and stops at the first event
of the call's type at or after the call's read position. All four problems
below make a real InDel fail or match the wrong event.

The function was run on CIGARs that carry exactly the called InDel, so the
right answer is always "found, distance 0" (a 40 bp read, InDel after read
base 17):

| case | haploid CIGAR | `haploid_check` | distance |
|---|---|---|---|
| same strand, 1 bp deletion | `17M1D23M` | True | 0 |
| same strand, 1 bp deletion, 2 bp insertion before it | `5M2I10M1D23M` | **False** | 17 |
| as above, plus a 1 bp deletion further on | `5M2I10M1D10M1D13M` | True | **8** |
| flipped strand, 1 bp deletion | `23M1D17M` | True | **1** |
| as above, 1 inserted base before it | `10M1I12M1D17M` | True | 0 |
| as above, 2 inserted bases before it | `10M2I11M1D17M` | **False** | 17 |
| flipped strand, 1 bp insertion | `22M1I17M` | True | 0 |
| flipped strand, 2 bp insertion | `21M2I17M` | **False** | 19 |

### H1. Insertions do not advance the walk

`I` operations consume read bases, but the walk does not count them: the line
`current_pos += length` is commented out
([line 106](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/filtering_haploid_mappping.py#L106)). Every insertion
before the call puts the walk behind by its length. The walk then reaches the
true event too early, skips it, and decides on whatever event of that type
comes next. The README describes the walk this way (*Haploid concordance*,
step 4). The tests in [test_pure_functions.py](../tests/test_pure_functions.py)
only use CIGARs with no insertion before the call.

In the data, same-strand 1 bp deletions are matched at distance 0 in 79% of
cases in the first tenth of a read, falling to 45% in its second half as
insertions accumulate upstream. The share recovers in the last tenth, where
little read is left for a wrong event to match. So most of the 39% of
same-strand 1 bp deletions in `results/` whose match lies more than 10 bp away
are probably real events that the walk skipped and that then passed on an
unrelated 1 bp deletion further on. For 2–3 bp InDels a later event of the same
size is rare, so a skipped call usually fails outright and leaves no trace in
`results/`.

On real reads (*Tested on real reads*) the pipeline's own function reproduces
the recorded `haploid_check` for all 10,487 calls of the chk window. The haploid
alignments carry a median of 33 InDels per read, so the walk falls behind
almost at once. Counting inserted bases is the only change needed to raise the
pass rate of same-strand calls:

| same-strand calls | pass now | pass with H1 fixed |
|---|---|---|
| 1 bp deletions | 93.7% | 99.9% |
| 2 bp deletions | 24.6% | 100% |
| 3 bp deletions | 22.8% | 100% |
| 1 bp insertions | 89.9% | 100% |
| 2 bp insertions | 14.1% | 100% |
| 3 bp insertions | 17.0% | 100% |

With H1 fixed the check passes 99.9% of calls. As coded, it removes about three
quarters of 2–3 bp calls, real read events, and very little else.

### H2. A strand flip maps the anchor, not the event

On a flipped read the call's position becomes `read_length − start − 1`
([lines 201–202](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/filtering_haploid_mappping.py#L201-L202)). That maps
the anchor base, but the event sits on one side of the anchor, and on the other
strand that side is reversed.

- **Deletions land 1 bp off.** Flipped 1 bp deletions are matched at distance 1
  in 48% of cases and at distance 0 in only 16%. Distance 0 happens when
  exactly one inserted base precedes them, so the two bugs cancel (H1).
- **An insertion of n bases lands n − 1 bases past the event**, so it is
  skipped for n ≥ 2. 16,036 flipped 2 bp insertions pass, against 198,958
  same-strand ones, and their median match lies 1,281 bp away.

With `max_hap_distance <- 0`, the notebook keeps 26% of flipped deletions
against 90% of same-strand ones.

### H3. On flipped reads, InDels in repeats land at opposite ends

Both alignments place a gap at its leftmost position along their own reference
(verified on real reads, see L1). When the haploid contig runs the other way,
its leftmost is the diploid's rightmost. So in a repeat the two alignments put the same InDel at
opposite ends, and the walk misses it. This hits every InDel that can slide,
not only those in masked tracts.

Flipping is a property of the read, so without this effect the flipped share of
calls would not depend on InDel size. Instead, across samples, flipped reads
carry 20–52% of the 1 bp calls and only 2–19% of the 2–3 bp calls. Under the
homopolymer rule a 1 bp call cannot slide. For deletions H2 only shifts the
distance by one, so the shortfall of flipped 2–3 bp deletions (11% of them,
against 38–40% for 1 bp InDels) comes from this effect. Most 2–3 bp calls
therefore sit in sequence where the gap can slide, which fits L1.

H2 and H3 together, measured per contig. Whether a read is flipped depends on
its contig: a contig whose haploid counterpart runs the other way flips every
read on it. Such contigs hold 33–48% of read positions. On them, the rate of
unique calls at distance 0 is 0.66–0.85× that of the other contigs for 1 bp,
0.14–0.25× for 2 bp, and 0.07–0.15× for 3 bp, in each of 7 samples
([check_haploid_check.py](review_checks/check_haploid_check.py), section 5).
For a direct test in reads, export a window of such a contig, e.g. human/chk
`h2tg000001l:0-10000000`.

**Suggested fix for H1–H3.** Walk the haploid CIGAR counting `I` as read bases.
Convert each event, and the call, to an interval of the read in its original
orientation, an interval rather than the anchor. Then accept an event of the
same type and size anywhere in the stretch of read the InDel can slide over.

### H4. `haploid_check` passes on a later event, and the distance mixes two measures

`haploid_check` is True if the first event of the call's type at or after its
position has the same size, however far away that event is. Final filtering
applies it without a distance limit, so `results/*/indels.tsv` and
`filter_report.tsv` count such calls as concordant.

`distance_to_indel_in_hap` is the smaller of two numbers: the distance from the
call to that event, and the gap between that event and the previous one of its
type ([lines 92 and 114](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/filtering_haploid_mappping.py#L92)). The
second is not measured from the call, so only a distance of 0 means what it
says. The notebook's `max_hap_distance <- 0` removes the wrong matches, but
because of H1–H3 it removes real calls too.

## Other issues

- **D2. Samples whose base qualities never reach 93 drop out.** human/la4 and
  human/ob003h have no read position at Q93: zero opportunities, and kmerQ ≥ 93
  removes all 4.7 and 6.5 million of their calls. Their quality scale differs
  (binned, or capped lower), so the fixed threshold
  ([indel_mutation_call.py:57](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/indel_mutation_call.py#L57)) means
  something else for them. The notebook calls la4's missing calls "an upstream
  problem"; this is the cause.
- **D3. k-mers cut short at read ends.** The opportunity test
  `1 <= read_pos <= query_length - 2`
  ([utilities.py:660](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L660)) is the same for all three
  sizes. At the last positions of a read, 2 and 3 bp k-mers run past its end
  and are counted under a 1 or 2 bp context, with a short kmerQ window. That is
  about three positions per read, so negligible (at most 1.2×10⁻⁴ of the
  opportunities in the chk window), but the test should be
  `read_pos + k <= query_length - 1` for each size.
- **D4. Calls and opportunities test kmerQ over different bases.** A 1 bp
  opportunity takes the lowest quality of three read bases: the one before, the
  base itself and the one after
  ([utilities.py:671](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L671)). A deletion call takes four
  read bases, two on each side of the gap
  ([utilities.py:743](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L743)). So 1 bp deletions need one
  more base at Q93 than their opportunities. On real reads, 34.9% of three-base
  windows and 30.7% of four-base windows are all Q93, so 1 bp deletion rates
  are low by about 12%. For 2–3 bp both windows have four bases, but not the
  same four.
- **U1. Unresolved-haplotype regions can shrink as they are built.** When a new
  region overlaps the current one, the code sets the end to the new region's
  end, not the larger of the two, and never moves the start left
  ([utilities.py:1228–1229](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L1228-L1229)). A later site
  whose reads end sooner cuts the region short, and the cut part is never
  masked. The README describes overlapping regions as extending each other.
  Fix: take the smaller start and the larger end, and treat a region that only
  touches the current one as overlapping. In the chk window, rebuilding the
  regions as coded reproduces the cluster's 769,672 bp exactly. With the bounds
  fixed they cover 10,566 bp more (1.4%), which would mask 23 more calls, 2 of
  them at kmerQ ≥ 93.
- **M1. Holes in the mask.** Overlap resolution drops a tract whole when it
  overlaps a kept one ([resolve_repeat_overlaps.py:62](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/resolve_repeat_overlaps.py#L62)).
  The mask is built from the resolved BEDs
  ([indel_mutation_call.py:132–134](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/indel_mutation_call.py#L132-L134)),
  so the dropped tract's bases outside the kept one stay unmasked unless another
  model covers them. The mask also never covers homopolymers (k = 1), periods of
  7 or more, satellites or segmental duplications, because TRF output is cut to
  k = 2–6 ([repetitive_dna.py:168](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/workflows/repetitive_dna.py#L168)).
  Building the mask from every tract before resolution closes the first kind of
  hole. In the chk window, exact repeats (k = 2–6, at least 3 copies) are
  0.7% unmasked, touching 9 calls and none at kmerQ ≥ 93, so the first kind is
  small. The never-annotated classes were not tested.
- **F1. No alignment-flag filter.** With `stepper="nofilter"`
  ([utilities.py:623](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/scripts/utilities.py#L623)), secondary, supplementary,
  duplicate and QC-fail alignments enter the calls and the opportunities alike.
  No read in `results/` carries calls at two different read lengths, so no
  calls come from hard-clipped supplementary alignments. In the chk window,
  13 of 10,487 calls (0.1%) sit on supplementary alignments, which are
  soft-clipped there, and none on secondary ones. Keeping primary alignments
  only, on both sides, removes the doubt.
- **N1. Calls carry no diploid mapQ.** The opportunities are split by diploid
  mapQ, and 0.5–40% of them (median 18%) are at mapQ 0 or 1. The calls carry
  only the haploid mapQ, so no mapQ filter can be applied to both. Writing
  `alignment.mapping_quality` into each call fixes this. In the chk window, 122
  of the 136 kmerQ ≥ 93 calls are on reads with diploid mapQ below 60, but with
  only 3 unique 1 bp calls there, a rate by mapQ needs more windows.
- **N2. Germline/somatic denominators include unclassified reads.** The notebook
  section *Rate per bin* splits all of a sample's read positions over the bins.
  But it drops the calls on reads missing from the likelihood table (7.5% of
  calls, about 16% in ob007), so the bin rates come out low by about that share.
- **N3. A printed tissue table mixes sizes.** In *Group: tissue samples*, the
  summary table sums calls over the three sizes but shows the 1 bp
  opportunities (`first(opportunities)`). The plot is correct.
- **N4. README and notebook disagree.**
  - The README gives `max_hap_distance` a default of `Inf`
    ([README.md:704](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/README.md#L704)); the notebook sets 0.
  - The README lists a `unique_only` setting, default `TRUE`, for both
    notebooks ([README.md:703](https://github.com/OliverKjHansen/detecting_InDels/blob/542fa20/README.md#L703)), but
    nonrepetitive_indels.ipynb has no such setting: it keeps the two classes
    apart instead.
  - The README's step 4 of *Haploid concordance* describes the H1 walk as
    intended.
- **N5. human/tsaed is unexplained.** All 1,808 of its final calls are
  `non-unique`, and 98% have a soft clip in the haploid alignment (notebook,
  *Haploid check*). Its haploid assembly or alignment is worth a look.

## Tested on real reads

[scripts/export_local_subset.py](../scripts/export_local_subset.py) copies one
window of a sample from the cluster, every alignment overlapping it, along with
what the scan needs and the cluster's results for that window. The findings
were tested on one such window: human/chk `h1tg000001l:0-10000000`, one complete
scan job with 13,917 alignments, 10,487 calls before the haploid check and 113
final calls.

**The copy is faithful.** Rerunning workflow 3 on it locally (`run_local.sh`)
gave the cluster's calls and opportunities line for line. The counts after the
haploid check (10,487) and after final filtering (113) also matched.

| ID | Test | Result |
|---|---|---|
| S1 | Correlation of subread support with base quality, per read | forward 0.32; reverse 0.000 as the scan indexes the tags, 0.32 in the read's own order |
| S1 | The tags at the call sites | the scan's values reproduced for 100% of calls; on reverse-strand reads, full support 70.0% as recorded against 37.7% correctly (forward 39.3%) |
| H1 | The pipeline's walk against the walk with inserted bases counted | the recorded `haploid_check` reproduced for 100% of calls; 2–3 bp pass rates rise from 14–25% to 100% |
| H2, H3 | Rates on contigs whose haploid counterpart runs the other way (7 samples, from `results/`) | 2 bp at 0.14–0.25×, 3 bp at 0.07–0.15×, 1 bp at 0.66–0.85× of the other contigs |
| H4 | Passing calls without a matching event, once H1–H3 are fixed | 12 of 7,740 (0.2%); 7 of them on 2 reads aligned to other haploid contigs |
| L1 | The base after the anchor in the mask; pbmm2's gap placement | 17.4% of unique 2 bp calls and 93% of kmerQ ≥ 93 2 bp calls (chance 0.92%); no InDel could slide left |
| U1 | The regions rebuilt as coded and with the bounds fixed | as coded = the cluster's BED; fixed +10,566 bp (1.4%) |
| M1 | Exact repeats (k = 2–6, ≥ 3 copies) left unmasked | 0.7% of their bases; 9 calls, none at kmerQ ≥ 93 |
| D1 | Read bases passing the 1 bp opportunity rule | 50.6%: 1 bp insertion rates about 2× too high |
| D3 | Truncated k-mers against opportunities | at most 1.2×10⁻⁴ |
| D4 | All-Q93 windows of three and four bases | 34.9% and 30.7%: 1 bp deletions low by 12% |
| F1 | Calls on supplementary or secondary alignments | 13 of 10,487 supplementary, no secondary |
| N1 | Diploid mapQ of the kmerQ ≥ 93 calls | 122 of 136 below 60; too few unique calls for a rate |
| yardstick | Q93 read positions passing the subread filter, read correctly | 28.0%, against about 32% of the controls used in [rate_floor.md](rate_floor.md) and [candidate_filters.md](candidate_filters.md) |

Not testable in this window:
- **H2 and H3, read by read.** The window's haploid contig runs the same way, so
  only 7 calls are on flipped reads. They were tested across samples instead.
- **E1, D2 and N2–N5.** These need many samples, another sample, or the
  notebook.

To rerun every test on a window, from the repository root:

```bash
python documents/review_checks/check_mirror_window.py mirror/human/chk/h1tg000001l_0_10000000
```

## Checks to run on the cluster

These need `output/`, which stays on the cluster. Each runs for one sample;
set `S`, and `BAM` to that sample's diploid BAM. Checks 1 and 3 were answered
for one window under *Tested on real reads*; run genome-wide, they give the
full size of L1 and of H1–H3.

**1. Measure L1 directly.** Count the calls whose first affected base, the base
after the anchor, lies in the mask. `position` is the anchor's 1-based
coordinate, so as a 0-based start it is that base. Without the leak, the share
should be close to the number of mask intervals over the number of unmasked
bases, about 1%.

```bash
S=human/da1
awk -F'\t' -v OFS='\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    {print $c["contig"], $c["position"], $c["position"] + 1, length($c["indel_string"]), $c["unique"]}' \
    results/nonrepetitive_indels/$S/indels.tsv | sort -k1,1 -k2,2n > next_base.bed
cat output/nonrepetitive_indels/$S/mask/*.bed | sort -k1,1 -k2,2n > mask.bed
bedtools intersect -u -a next_base.bed -b mask.bed | cut -f4 | sort | uniq -c   # next base masked, by size
cut -f4 next_base.bed | sort | uniq -c                                           # all calls, by size
```

**2. Check the denominator over all kmerQ values.** Compare the size-3
opportunities in the unfiltered table with the read bases at unmasked positions
inside the coverage window. The two should agree to within about 10%: the
5-mer homopolymer rule removes about 5%, and `samtools depth` filters reads
slightly differently from the scan.

```bash
awk -F'\t' 'NR > 1 && length($2) == 3 {s += $6} END {printf "%.0f\n", s}' \
    output/nonrepetitive_indels/$S/merged/opportunities.tsv
bedtools subtract -a output/precomute/$S/long_contigs.bed -b mask.bed > unmasked.bed
C=$(cat output/precomute/$S/avg_coverage.txt)
samtools depth -b unmasked.bed $BAM | awk -v c=$C '$3 >= c/2 && $3 <= 2*c {s += $3} END {printf "%.0f\n", s}'
```

**3. Count the real calls the haploid check removes (H1–H3).** Re-run
`filtering_haploid_mappping.py` with a fixed walk on `calls/*.tsv` and
`haploid/locations.tsv`. Compare the `haploid_check` pass rate by size and
strand flip with `merged/indels.tsv`.

**4. Find the quality scale of la4 and ob003h (D2).** Print the highest base
quality in the first reads of the BAM:

```bash
samtools view $BAM | head -n 2000 | cut -f 11 \
    | python3 -c "import sys; print(max(ord(c) - 33 for line in sys.stdin for c in line.strip()))"
```

## Reproducing the numbers

Run from the repository root, with `results/` copied from the cluster:

```bash
.venv/bin/python documents/review_checks/check_denominator.py
```

| Script | Prints | Behind |
|---|---|---|
| [check_denominator.py](review_checks/check_denominator.py) | Q93 opportunities against the bases sequenced; 1 bp against 3 bp opportunities; the mapQ ≤ 1 share; samples with no Q93 position | the denominator section, D1, D2, N1 |
| [check_repeat_leak.py](review_checks/check_repeat_leak.py) | 2 bp rate per context; outside against inside M0 and M1 k = 2 tracts; how much the missing events explain; tracts per model | L1 |
| [check_error_floor.py](review_checks/check_error_floor.py) | rate per sample and size; the share kept by kmerQ ≥ 93 | E1 |
| [check_haploid_check.py](review_checks/check_haploid_check.py) | the CIGAR cases; strand flip by size; distance to the matched event; distance 0 along the read; rates on flipped and same-strand contigs | H1–H4 |
| [check_rate_floor.py](review_checks/check_rate_floor.py) | rate by read length, passes, subread support, SMRT Cell and instrument; rates after the cleanest filters | E1, [rate_floor.md](rate_floor.md) |
| [check_candidate_filters.py](review_checks/check_candidate_filters.py) | rate by further call properties; the strand split of the subread signal; rates with a working subread filter; clustering; contigs with an excess; supplementary alignments | S1, F1, [candidate_filters.md](candidate_filters.md) |
| [check_mirror_window.py](review_checks/check_mirror_window.py) | every test on real reads, from one mirrored window (pass the window's folder) | *Tested on real reads* |

[common.py](review_checks/common.py) holds the loading code and the filters
under *How the numbers were made*. Figures quoted from the notebook come from
the run saved in [nonrepetitive_indels.ipynb](../plotting/nonrepetitive_indels.ipynb).
