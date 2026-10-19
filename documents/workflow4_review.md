# Review: InDel rates inside repeats (workflow 4)

*19 October 2026. Covers workflow 4
([repetitive_indel_call.py](../workflows/repetitive_indel_call.py)): the repeat
scan (`tandem_repeat_scan()` in [utilities.py](../scripts/utilities.py)), the
haploid check and the final filters as workflow 4 applies them, and the tracts it
reads from workflow 2 ([resolve_repeat_overlaps.py](../scripts/resolve_repeat_overlaps.py)).
Code at commit a787976 on `false-positive-reduction`, which the links to lines
of code point at. Data: the tables in `results/` copied from the cluster on 30
September 2026. Rates come from
[check_workflow4.py](review_checks/check_workflow4.py). How the code treats
real reads comes from [check_workflow4_window.py](review_checks/check_workflow4_window.py)
on the human/chk window mirrored for ticket 07; the mirror is for software
testing only, so no rate is drawn from it.*

*This is the counterpart of [nonrepetitive_indel_review.md](nonrepetitive_indel_review.md)
(workflow 3), whose finding IDs (L1, S1, H1–H4, …) are used here. Findings
new to workflow 4 are R1–R9.*

## Summary

Workflow 4 has never been reviewed, and its repeat rates are not usable as
they stand. They are too low, and by different amounts in different tracts,
so neither the rates nor the comparison between M0, M1 and M2 can be read at
face value.

- **It barely sees slippage in pure tracts (R1).** pbmm2 places an InDel of
  whole units at its leftmost position, which in a pure tract is the base
  before the tract. Workflow 4 scans only bases inside tracts. In the human/chk
  window, 2,262 one-unit InDels of an M0 tract's motif sit on the base before
  the tract and 29 inside it. Genome-wide, within M1, one-unit InDels per base
  are 0.07–0.10× as frequent in exactly pure dinucleotide tracts as in
  90–99.9% pure ones of the same length.
- **What it does count is often the next tract's slippage (R2).** When two
  tracts touch, an InDel of the second slides onto the last base of the first
  and is credited to it. 171 of the window's 225 M0 calls are not the tract's
  own motif.
- **The haploid check removes most of what is left (R3).** It drops a median
  80–85% of repeat calls, against 31% outside repeats. Only 2–11% of the
  final repeat calls are on reads whose haploid alignment runs the other way,
  against 39% of 1 bp calls outside repeats.
- **The opportunity unit and the tract choice don't match the calls (R4–R8).**
  Opportunities are counted per read base, over every read, under a quality
  window sized by the period. Calls are counted per event, under a window sized
  by the InDel. Calls whose size isn't a multiple of the period are dropped
  silently.
- **Overlap resolution drops whole tracts (R9).** A tract that overlaps one of
  smaller period loses even the part that overlaps nothing, in both workflows.

Most of this is already planned. Tickets 17, 18, 22 and 23 fix R1–R3 and R5–R8,
and R9 is new as ticket 30. Ticket 22's opportunity unit has to settle R7 and
R8 when it is built. Those notes are added to the tickets.

| ID | Finding | Where | Effect on repeat rates | Ticket |
|---|---|---|---|---|
| R1 | Slippage in a pure tract is anchored before the tract and never scanned | [utilities.py:510](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/utilities.py#L510-L517) | too low, most in pure tracts: about 10× in pure M1 tracts | 17 |
| R2 | A tract is credited with the slippage of the tract after it | as R1 | wrong tract, wrong motif | 17, 23 |
| R3 | The haploid check removes 80–85% of repeat calls (H1–H3) | [filtering_haploid_mappping.py:106](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/filtering_haploid_mappping.py#L106), [:202](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/filtering_haploid_mappping.py#L202) | too low, several-fold | 18 |
| R4 | The subread tags are read at the wrong base on reverse-strand reads (S1) | [utilities.py:598](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/utilities.py#L598-L602) | none directly; no subread filter possible | 12 |
| R5 | InDels whose size isn't a multiple of the period are dropped silently | [utilities.py:574](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/utilities.py#L574-L582) | unmeasured events; multiples not checked for the motif | 23 |
| R6 | One opportunity per read base, partial reads included | [utilities.py:546](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/utilities.py#L546-L562) | wrong unit; rates by tract length approximate | 22 |
| R7 | Opportunities test base quality over a window sized by the period, calls by the InDel | [utilities.py:558](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/utilities.py#L558-L559) | multi-unit calls (12–31%) tested on other bases | 22 |
| R8 | The largest-k / smallest-k choice never applies | [utilities.py:534](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/utilities.py#L534), [:574](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/utilities.py#L574) | none: one tract per base per model | 22 |
| R9 | Overlap resolution drops a tract whole, not just its overlap | [resolve_repeat_overlaps.py:62](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/resolve_repeat_overlaps.py#L62) | repeat bases in no tract of the model; holes in workflow 3's mask (M1) | **30 (new)** |
| — | The 100 bp read-end cut is applied to calls only | [repetitive_indel_call.py:253](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/workflows/repetitive_indel_call.py#L253) | about 1.5% too low | 13, 24 |
| — | Unresolved haplotype regions, homopolymers ≥ 6 bp and periods ≥ 7 are not excluded | the mask is workflow 3's only | too high where a germline length difference mismaps | 20, 21 |
| — | F1 and N1 apply as in workflow 3 | | small | N1: 14; F1: none |

## R1. Slippage in a pure tract is anchored before the tract

**Mechanism.** The scan looks at a pileup column only if its position lies
inside a tract of some model ([utilities.py:510–517](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/utilities.py#L510-L517)),
and a call belongs to the tract containing its anchor, the base before the gap.
pbmm2 places an InDel that can slide at its leftmost position (verified on real
reads in the workflow 3 review, *L1*). Inside a pure tract a whole-unit InDel
can slide to the tract's first unit, so its anchor is the base before the
tract. That base is in no tract, so workflow 4 never sees the event. Workflow
3 sees it, because its mask is tested at the anchor only (L1). In an impure
tract, an impurity can stop the slide inside the tract, and the event is
counted.

**Direction.** Repeat rates too low, most of all in pure tracts and so in M0,
whose tracts are all pure. This is the other half of L1: the events workflow 4
misses are the ones workflow 3 counts as outside repeats.

**Evidence.**
- **On real reads** (check_workflow4_window.py, human/chk window). Of the
  one-unit InDels of an M0 tract's motif in the reads' alignments, 2,262 are
  anchored on the base before the tract and 29 inside it. Workflow 4 can see
  about 1% of them.
- **Genome-wide, within M1** (check_workflow4.py, section 4). For k = 2,
  one-unit calls per base of tract sequence in exactly pure tracts are 0.09×
  those in tracts 90% to < 100% pure. That's the median over 33 samples (range
  0.03–0.18). Slippage should be faster in purer tracts, not 11× slower.
  - **Tract length doesn't explain it.** Pooled over samples within bins of
    copy number, the ratio is 0.10 for up to 15 copies, 0.08 for 15–20 and
    0.07 for 20–30. Exactly pure M1 dinucleotide tracts of more than 30 copies
    carry essentially no calls at all.
  - **The denominator** is tract sequence from the summaries, because the
    opportunities carry no purity; depth cancels out of a ratio within a
    sample.
- **Genome-wide, models compared** (section 3). One-unit InDels per Q93 read
  base, M1 over M0, median over samples: 106× for k = 2, 263× for k = 3, 31×
  for k = 4, 68× for k = 5, 6× for k = 6. The medians for k = 2 are 6.2×10⁻⁶ in
  M1 and 4.7×10⁻⁸ in M0. M0 also has many short tracts, which slip slowly, so
  this comparison alone mixes R1 with tract length. It is consistent with the
  comparison within M1, which doesn't.

**Test.** Ticket 17's failing test: a pure (CA)n tract and a read with a
one-unit deletion placed leftmost, anchored on the base before the tract.
Workflow 4 must report it in that tract, and workflow 3 must not report it.

## R2. A tract is credited with the slippage of the tract after it

**Mechanism.** R1 with two tracts that touch. Where one tract runs straight
into another, for example (AC)n then (AT)n, the base before the second tract is
the last base of the first. A one-unit InDel of the second slides onto it and
is called in the first tract. If its size is a multiple of the first tract's
period, it passes the size test (R5) and is counted there under the first
tract's motif.

**Direction.** Calls in the wrong tract, under the wrong motif and period.
Rates per motif are mixed up; the total per model is too high by these calls
and too low by R1.

**Evidence** (check_workflow4_window.py, human/chk window, the local replay of
workflow 4). Of the window's 225 M0 calls, 172 are anchored in the last unit
of their tract and 5 on its first base. 171 of the 225 are not whole units of
the tract's motif. For example, three deletions of AT are credited to an (AC)18
tract that is followed directly by (TA)n.

**Test.** Two adjacent tracts, (AC)n then (AT)n, and a read with an AT deletion
anchored on the last base of the (AC)n tract. It belongs to the (AT)n tract,
as slippage. Today it is an (AC)n call. Tickets 17 (which tract) and 23
(slippage or other) settle it together.

## R3. The haploid check removes 80–85% of repeat calls

**Mechanism.** H1–H3 of the workflow 3 review, on InDels that can almost
always slide. Insertions before the call don't advance the CIGAR walk (H1). A
strand flip maps the anchor, not the event (H2). On flipped reads, the two
alignments put a repeat InDel at opposite ends of the stretch it can slide over
(H3), so the walk misses it whatever H1 and H2 do.

**Direction.** Repeat rates too low, several-fold, and unevenly: calls on
flipped reads are lost almost entirely.

**Evidence** (check_workflow4.py, sections 1 and 2).
- **Share of calls the check removes** (`filter_report.tsv`): a median of
  85.2% in M0, 80.0% in M1 and 83.1% in M2, against 30.7% in workflow 3.
- **Calls on flipped reads.** A read is flipped whatever InDel it carries, so
  without a strand-dependent loss every class would show the same share of
  flipped reads. Among the final calls it is 39% for 1 bp InDels outside
  repeats, which cannot slide, and 2–11% for every model and period inside
  repeats.

**Test.** Ticket 18's CIGAR cases, plus one for a repeat: a whole-unit
deletion in a tract on a flipped read, placed leftmost on each strand. It must
be found at distance 0.

## R4. S1 in the repeat scan

**Mechanism.** The repeat scan reads `sm`, `sx` and `sa` the same way as
workflow 3 ([utilities.py:598–602](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/utilities.py#L598-L602)),
so on reverse-strand reads they describe the mirrored base.

**Direction.** None on the rates; a subread filter can't be used in repeats
either.

**Evidence** (check_workflow4.py, section 7). Unique M1 calls with full support
(`sm` equal to `sa` and no `sx`): 29.4% on forward-strand reads and 40.6% on
reverse-strand ones, median over 27 Revio samples. Outside repeats, the same
inversion is what showed S1.

**Test.** Ticket 12's reverse-strand fixture, run through the repeat scan too.

## R5. InDels whose size isn't a multiple of the period are dropped

**Mechanism.** A call is kept only if its size is a multiple of the period of a
tract at the anchor ([utilities.py:574–582](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/utilities.py#L574-L582));
anything else is not written anywhere. A multiple is kept whatever its
sequence, so R2's calls pass.

**Direction.** Events in repeats that aren't whole-unit slippage are
unmeasured, and some whole-unit calls aren't the tract's motif.

**Evidence.** By construction; the dropped calls leave no trace in `results/`.
In the chk window, 171 of 225 M0 calls are multiples of the period but not the
motif (R2).

**Test.** Ticket 23's three cases.

## R6. One opportunity per read base

**Mechanism.** An opportunity is a read base inside a tract
([utilities.py:546–562](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/utilities.py#L546-L562)),
on every read that covers it, including reads that end inside the tract. A
call is an event in a tract. A read that ends inside the tract cannot show a
change in its length, yet it adds opportunities; and a long tract adds more
opportunities per read than a short one, though each read can show at most a
few events in it.

**Direction.** The rate's unit is wrong, and rates by tract length or purity
can only be approximated. The direction for a given class depends on its
tract lengths.

**Evidence.** The design itself; [test_read_shapes.py](../tests/test_read_shapes.py)
pins today's behaviour: a read ending inside a tract adds 49 opportunities
for the 49 tract bases it covers.

**Test.** Ticket 22's partial-read fixtures.

## R7. Opportunities and calls test quality over different bases

**Mechanism.** An opportunity's quality window is two bases on each side of a
window of k bases, with k the period
([utilities.py:558–559](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/utilities.py#L558-L559)).
A call's window is two bases on each side of its gap. For a call of one unit
the two are of the same size (D4 applies as in workflow 3); for a call of two
or more units they test different bases.

**Direction.** Small for one-unit calls. For multi-unit calls the
kmerQ ≥ 93 filter does not have a matching denominator.

**Evidence** (check_workflow4.py, section 6). Calls of two or more units are a
median of 30.6% of the final M0 calls and 11.9–12.5% of M1 and M2.

**Test.** Comes with ticket 22: a per-tract opportunity needs its own quality
rule, and the calls of that tract must pass the same rule over the same bases.

## R8. The largest-k / smallest-k choice never applies

**Mechanism.** When a position lies in several tracts of one model, calls take
the tract with the largest period that divides their size and opportunities
the one with the smallest period
([utilities.py:534](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/utilities.py#L534), [:574](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/utilities.py#L574)).
But workflow 2 resolves overlaps within each model before workflow 4 reads the
tracts, so no position is in two tracts of one model.

**Direction.** None today. The rule would bias rates if overlapping tracts ever
reached the scan, because calls and opportunities would use different tracts.

**Evidence.** The resolved BEDs of all four mirrored windows: 0 overlapping
pairs in M0, M1 or M2 (12 BEDs, up to 85,882 tracts each).
[test_read_shapes.py](../tests/test_read_shapes.py) shows the rule on
overlapping tracts made by hand.

**Test.** When ticket 22 replaces the unit, one tract per base per model can be
asserted, and the two rules collapse into one.

## R9. Overlap resolution drops a tract whole

**Mechanism.** [resolve_repeat_overlaps.py](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/scripts/resolve_repeat_overlaps.py#L62)
keeps tracts of the smallest period first and drops every later tract that
overlaps a kept one, including the bases of the dropped tract that overlap
nothing. On a constructed BED with an AC tract at 100–140 and an AAG tract at
130–190, it keeps only the AC tract: bases 140–190 are in no tract of that
model.

**Direction.** In workflow 4, repeat sequence of larger periods is missing from
the model, with its calls and opportunities: rates by period are biased in an
unknown direction. In workflow 3, those bases are outside the mask unless
another model covers them, so their slippage counts as non-repetitive. This is
the first kind of hole in the workflow 3 review's M1, which no ticket in the
plan closes; ticket 21 closes the second kind.

**Evidence.** The constructed example above. Genome-wide size: not in
`results/`, which holds only resolved tracts. On the cluster, per sample and
model:

```bash
S=human/chk; M=M0
cat output/repetitive_dna/$S/per_k/*_${M}_k*.bed | grep -v '^contig' | cut -f1-3 \
    | sort -k1,1 -k2,2n | bedtools merge -i - | awk '{s += $3 - $2} END {print s}'   # bp in any tract
awk 'NR > 1 {s += $3 - $2} END {print s}' output/repetitive_dna/$S/$M.bed            # bp kept
```

**Test.** Ticket 30: the constructed BED above; every base of a dropped tract
outside the kept one must end up in a defined group, not in no group.

## Other issues

- **The read-end cut is applied to calls only.** Workflow 4's final filters drop
  calls within 100 bp of a read end
  ([repetitive_indel_call.py:253](https://github.com/OliverKjHansen/detecting_InDels/blob/a787976/workflows/repetitive_indel_call.py#L253)),
  but the opportunities keep those bases. The cut removes a median 1.5–1.8% of
  the calls; with a median read length of 13 kb, about 1.5% of read bases lie
  within 100 bp of an end, so rates are about 1.5% too low. Workflow 3 doesn't
  apply the cut at all. Tickets 13 (read-end bins on both tables) and 24 (one
  filter set, applied to both) close it.
- **No exclusion of unresolved regions, homopolymers ≥ 6 bp or periods ≥ 7.**
  Only workflow 3 masks anything. Tickets 20 and 21.
- **F1 and N1.** Workflow 4 has no alignment-flag filter either, and its calls
  carry no diploid mapQ. N1 is ticket 14. F1 has no ticket in the plan; the
  workflow 3 review measured it as small (0.1% of calls on supplementary
  alignments), and nothing here suggests repeats differ.
- **The haploid mapping is shared** between the models and with nothing else,
  as intended: one remapping per sample for all three models.

## Reproducing the numbers

From the repository root:

```bash
.venv/bin/python documents/review_checks/check_workflow4.py
.venv/bin/python documents/review_checks/check_workflow4_window.py mirror/human/chk/h1tg000001l_0_10000000
```

The second needs the mirrored window, and its section on calls needs the
window replayed first (`mirror/human/chk/h1tg000001l_0_10000000/run_local.sh`).
