import pysam
import csv
import math
from bisect import bisect_right
from collections import defaultdict, namedtuple
from intervaltree import IntervalTree

from read_bins import load as load_read_bins
from read_tags import (Blocks, describe, events, indels, nearest_other_event, nearest_outside, subread_support,
                       subread_tags, tag)
from strata import Site, Strata, Tract

# What a call and the opportunities both need of one read, computed once per
# read: where its alignment ends, its sequence and base qualities, its read
# keys, its subread tags and support in the alignment's order, and the
# assembly positions of its own InDels and mismatches.
ReadEntry = namedtuple("ReadEntry", "end sequence qualities key tags support events")


def read_entry(alignment, tissue_table, classification_table, strata, fasta, contig):
    """The ReadEntry of an alignment."""
    tags = subread_tags(alignment)
    reference = fasta.fetch(contig, alignment.reference_start, alignment.reference_end).upper()
    name = alignment.query_name
    return ReadEntry(alignment.reference_end, alignment.query_sequence, alignment.query_qualities,
                     strata.of_read(describe(alignment, tissue_table.get(name, "NA"),
                                             classification_table.get(name, "NA"))),
                     tags, subread_support(tags), events(alignment, reference))


def read_tag_columns(alignment):
    """The calls' raw tag columns np, ec, rq, fn, rn and mapQ, NA where the read lacks a tag."""
    values = [tag(alignment, name) for name in ("np", "ec", "rq", "fn", "rn")] + [alignment.mapping_quality]
    return ["NA" if value is None else value for value in values]

class Mask:
    """
    Workflow 3's mask: the BED's intervals, merged, and how far the unmasked
    sequence runs from a position. With window_start/window_end only the
    intervals a scan of [window_start, window_end) can meet are kept: from the
    base before the window (an opportunity's anchor) to `reach` bases after it
    (the bases an event anchored in the window touches).
    """

    def __init__(self, bedfile, window_start=None, window_end=None, reach=0):
        intervals = []
        with open(bedfile) as handle:
            for line in handle:
                columns = line.strip().split("\t")
                start, end = int(columns[1]), int(columns[2])
                if window_start is not None and end < window_start:
                    continue
                if window_end is not None and start >= window_end + reach:
                    continue
                intervals.append((start, end))
        self.starts, self.ends = [], []
        for start, end in sorted(intervals):
            if self.ends and start <= self.ends[-1]:
                self.ends[-1] = max(self.ends[-1], end)
            else:
                self.starts.append(start)
                self.ends.append(end)

    def clear(self, position):
        """How many bases from `position` on are unmasked: 0 if it is masked, infinite if no masked base follows."""
        i = bisect_right(self.starts, position) - 1
        if i >= 0 and self.ends[i] > position:
            return 0
        return self.starts[i + 1] - position if i + 1 < len(self.starts) else math.inf


def load_tissue_classification(tissue_file):
    """
    readID -> tissue, from a two-column TSV.

    Samples without read labels pass "NA" instead of a file; they get an empty
    table, so every read is reported with tissue NA.
    """
    if tissue_file == "NA":
        return {}
    with open(tissue_file, "r") as f:
        return dict(csv.reader(f, delimiter="\t"))


def canonical_context(sequence):
    """
    Group a context with the same context read from the other strand: the
    smaller of the sequence and its reverse complement.

    A GT deletion on one strand is an AC deletion on the other, so both are
    reported as AC.
    """
    return min(sequence, reverse_complement(sequence))

def canonical_repeat_indel(sequence):
    """
    Canonical form of an InDel inside a tandem repeat: the smallest rotation of
    the sequence or of its reverse complement.

    Inside a (TC)n tract, gaining one unit reads TC or CT depending on where the
    aligner placed the gap, and GA or AG read from the other strand. All four
    are the same event. This is also how parse_trf.py canonicalises the tract
    motif, so an InDel and the tract it sits in carry the same label.

    Rotation only makes sense when the sequence is a whole number of repeat
    units; contexts outside repeats use canonical_context() instead.
    """
    strands = (sequence, reverse_complement(sequence))
    return min(strand[i:] + strand[:i] for strand in strands for i in range(len(strand)))

def homopolymer_check(seq):
    """
    Checks if a sequence is mostly a homopolymer.
    Returns True if the sequence has at most one break in homopolymer continuity.
    """
    # if len(set(seq)) == 1:
    #     return True  # entire sequence is a homopolymer

    breaks = 0

    for i in range(len(seq) - 1):
        if seq[i] != seq[i + 1]:
            breaks += 1
            if breaks == 2:
                return False

    return True

def has_flanks(before, after, read_length):
    """
    Whether an event has two read bases on each side: `before` is the read
    index of the last base before it, `after` of the first base after it.
    Calls and opportunities without them are neither called nor counted (D3).
    """
    return before >= 1 and after + 2 <= read_length

def flanking(values, before, after):
    """
    Per-base values (base qualities or subread support) of the two read bases
    before an event and the two after it, never of the event's own bases: the
    bases calls and opportunities alike are tested on (D4).
    """
    return values[before - 1:before + 1] + values[after:after + 2]

def nonrepetitive_scan(filters, contig, bamfile, fastafile, average_coverage_file, tissue_file,
                       window_start=None, window_end=None,
                       calls_out=None, opportunities_out=None, fastq_out=None,
                       *, coverage_band, indel_sizes, key_parameters, flank_width, read_bins_file="NA"):
    """
    Call InDels of the sizes in indel_sizes, (smallest, largest) in bp, and
    count opportunities outside the mask in one pass over the BAM, and write
    the reads carrying a call as FASTQ.

    calls_out, opportunities_out and fastq_out are open files; whichever is
    None is not produced.

    Only calls and opportunities whose anchor has a depth within
    coverage_band, (low, high) multiples of the sample's mean coverage, are
    called and counted. key_parameters are the `strata` block of config.yaml.
    Each call carries flank_width bases of the assembly on each side of its
    anchor. read_bins_file is the compact table of the sample's read bins
    (read_bins.py, ticket 44), or "NA".

    Calls: InDels of those sizes at positions outside the mask (the filters
    BED), except those whose sequence plus one flanking base on each side is a
    homopolymer (changes base at most once).

    Opportunities: at every read position outside the mask, a deletion of each
    size, counted under the canonical k-mer it would delete (the smaller of
    the k-mer and its reverse complement) with the same homopolymer rule, and
    an insertion of each size in the gap before the position, counted under
    the canonical pair of bases around the gap. Counted per type, size and
    value of the keys in strata.py; each call carries the same keys as
    columns, and its `context` is the opportunity context it is counted
    against.

    Calls and opportunities alike are tested (kmerQ, subread support) on the
    two read bases on each side of the event, never on the event's own bases,
    and are neither called nor counted without them (D3, D4).

    fastq_out gets every read carrying a call, once, as forward-strand FASTQ
    for haploid remapping.

    With window_start/window_end only positions in [window_start, window_end)
    are scanned; each position belongs to one window, so windows never report
    the same call twice and their opportunity counts add up to the whole
    contig. Without them the whole contig is scanned.
    """
    with open(average_coverage_file, 'r') as f: # load in the average coverage smaple
        average_coverage_of_species = float(f.read().strip())
    low_depth, high_depth = (average_coverage_of_species * factor for factor in coverage_band)
    smallest, largest = indel_sizes
    sizes = range(smallest, largest + 1)

    fasta = pysam.FastaFile(fastafile)
    mask = Mask(filters, window_start, window_end, reach=largest)
    d = load_tissue_classification(tissue_file)
    classified = load_read_bins(read_bins_file)
    strata = Strata(**key_parameters)

    min_mapping_qual = 0 # to have more contorl over pysam
    min_base_qual = 0 #--||---

    # (context, *read keys, *site keys) -> count
    data = defaultdict(int)

    if calls_out is not None:
        print("#contig", "start", "end", "ref", "alt",
        "indel_type","unique","indel_string","canonical_indel","context","read_length",
        "np","ec","rq","fn","rn","mapQ","sm","sx","sa","sm_pre","sx_pre","ToreadEnd","kmerQ",
        "position","readID", "diploid_direction", "flank_sequence", *strata.columns,
        sep= "\t", file=calls_out) #header for output

    # Each read's ReadEntry, built once rather than at every position it
    # covers: copying a whole HiFi read at each position is most of the cost
    # of counting opportunities. Keyed by read, start and flag, with the read's
    # end so reads the scan has passed can be dropped every 10 kb.
    reads = {}
    next_eviction = 0

    def entry_of(alignment):
        key = (alignment.query_name, alignment.reference_start, alignment.flag)
        entry = reads.get(key)
        if entry is None:
            entry = reads[key] = read_entry(alignment, d, classified, strata, fasta, contig)
        return entry

    # The last column's position and depth: the depth at the anchor of the
    # opportunities of the column after it.
    previous_position, previous_depth = None, 0

    # Canonical k-mer and homopolymer check of a run (k-mer plus its flanking
    # bases), and canonical pair of bases around a gap, computed once each.
    runs = {}
    gaps = {}

    fastq_reads = set()

    with pysam.AlignmentFile(bamfile, "rb") as bam: # print the header
        # From the base before the window, for its depth.
        for pileupcolumn in bam.pileup(contig=contig, #pileup
                                       start=None if window_start is None else max(0, window_start - 1),
                                       stop=window_end,
                                       fastafile=fasta,
                                       truncate=True,
                                       stepper="nofilter",
                                       min_base_quality=min_base_qual,#quite light
                                       min_mapping_quality=min_mapping_qual,#conservative
                                       ):

            position = pileupcolumn.reference_pos # right. the position in the reference

            # Coverage is tested at the anchor, for calls and opportunities
            # alike: a call is anchored here, an opportunity on the base before
            # (ticket 31). A base without reads has depth 0.
            depth = pileupcolumn.nsegments
            anchor_depth = previous_depth if previous_position == position - 1 else 0
            previous_position, previous_depth = position, depth
            if window_start is not None and position < window_start:
                continue # the base before the window, read for its depth only
            # Whether this column's calls, and its opportunities, are made.
            calling = calls_out is not None and low_depth <= depth <= high_depth
            counting = opportunities_out is not None and low_depth <= anchor_depth <= high_depth
            if not (calling or counting):
                continue

            # An event of k bases anchored at base a, call or opportunity, is
            # counted only if a to a+k are unmasked: InDels that slide to the
            # base before a tract are the tract's (L1). A call is anchored
            # here, an opportunity on the base before.
            clear = mask.clear(position)
            if not clear:
                continue
            anchor_clear = mask.clear(position - 1)

            if position >= next_eviction:
                reads = {key: entry for key, entry in reads.items() if entry.end > position}
                next_eviction = position + 10_000

            pileups = pileupcolumn.pileups

            if calling:
                indel_count = sum(1 for pileupread in pileups if (pileupread.indel != 0) or (pileupread.is_del == True))

                seen = "unique"
                if indel_count > 1:
                    seen = "non-unique"

            for pileupread in pileups:
                if pileupread.query_position is None: # we are in an indel and the position is none
                    continue

                read_pos   = pileupread.query_position
                alignment  = pileupread.alignment
                query_name = alignment.query_name

                # ---- opportunities ----
                # The gap before this read base, for an insertion, and the
                # k-mers starting at it, for deletions, each tested on two read
                # bases on each side, as the calls are.
                if counting and read_pos >= 2 and anchor_clear > smallest:
                    entry = entry_of(alignment)
                    read_seq, base_qualities, support = entry.sequence, entry.qualities, entry.support
                    read_length = len(read_seq)
                    # An opportunity's anchor is the base before its gap or k-mer, as a call's is the base before its gap.
                    anchor = position - 1
                    to_read_end = min(anchor - alignment.reference_start, entry.end - anchor)
                    nearest = nearest_other_event(entry.events, anchor)

                    def keys(before, after):
                        site = Site(flanking(base_qualities, before, after),
                                    None if support is None else flanking(support, before, after), to_read_end, nearest)
                        return entry.key + strata.of_site(site)

                    #quality, insertion   deletion of 1    deletion of 3
                    #NT|AC                 NT|A|CN          NT|ACT|CN
                    #00|00                 00|x|00          00|xxx|00
                    if has_flanks(read_pos - 1, read_pos, read_length):
                        pair = read_seq[read_pos-1:read_pos+1]
                        context = gaps.get(pair)
                        if context is None:
                            context = gaps[pair] = canonical_context(pair) # grouped with the other strand
                        gap = keys(read_pos - 1, read_pos) # the same bases tested for every size
                        for indel_length in sizes:
                            if indel_length >= anchor_clear:
                                break # it would touch the mask, as would any longer one
                            data[("insertion", indel_length, context) + gap] += 1

                    for indel_length in sizes: #indel lengths here means kmer length
                        if indel_length >= anchor_clear or not has_flanks(read_pos - 1, read_pos + indel_length, read_length):
                            break # nor does any longer k-mer
                        run = read_seq[read_pos-1:read_pos+indel_length+1]   #run dependent on indel length # includes kmer + flanking bases
                        info = runs.get(run)
                        if info is None:
                            #context of kmer without the flanking bases, grouped with the other strand
                            info = runs[run] = (canonical_context(run[1:-1]), homopolymer_check(run))
                        lexi, homopolymer = info
                        if homopolymer == True: # filter out homopolymers
                            continue
                        data[("deletion", indel_length, lexi) + keys(read_pos - 1, read_pos + indel_length)] += 1

                # ---- calls ----
                if not calling:
                    continue

                indel = pileupread.indel

                if not smallest <= abs(indel) <= largest:
                    continue

                indel_length = abs(indel)
                if indel_length >= clear:
                    continue # a base it touches is masked (L1)
                # The last read base before the gap and the first after it.
                before, after = read_pos, read_pos + 1 + (indel_length if indel > 0 else 0)
                if not has_flanks(before, after, alignment.query_length):
                    continue

                sm = sx = sa = sm_pre = sx_pre = "NA" # initialise tags
                entry = entry_of(alignment)
                readstart = alignment.reference_start
                ToreadEnd = min([(position-readstart),(entry.end-position)]) #tjek
                nearest = nearest_other_event(entry.events, position)
                flank = fasta.fetch(contig, max(0, position - flank_width), position + flank_width + 1)

                read_seq = entry.sequence
                base_qualities = entry.qualities

                # The subread tags in the alignment's order (S1).
                tags = entry.tags
                call_support = entry.support
                if tags is not None:
                    sm_all, sx_all, sa_all = tags
                    sm, sx = sm_all[read_pos+1], sx_all[read_pos+1] # the base after the anchor
                    sa, sm_pre, sx_pre = sa_all[read_pos], sm_all[read_pos], sx_all[read_pos] # the anchor

                diploid_strand = "-" if alignment.is_reverse == True else "+"

                #quality check, as for the opportunities
                #AA|ACG|TT
                #00|xxx|00
                site = Site(flanking(base_qualities, before, after),
                            None if call_support is None else flanking(call_support, before, after), ToreadEnd, nearest)

                if indel < 0:
                # -------------------
                # CASE 1: DELETIONS
                # -------------------
                    hcheck = fasta.fetch(contig, position, position + indel_length+2) #for homopolymer check
                    deletion_context = hcheck[1:-1] #the deletion context, as observed

                    homopolymer = homopolymer_check(hcheck)
                    if homopolymer == True: # filter out homopolymers
                        continue

                    row = (contig, read_pos+1,read_pos+1,hcheck[0:-1],hcheck[0], #start and end are the same for delettion, becaus its is in the read
                           "deletion", seen, deletion_context, canonical_context(deletion_context),
                           canonical_context(deletion_context), len(read_seq),
                           *read_tag_columns(alignment),sm,sx,sa,sm_pre,sx_pre,ToreadEnd,min(site.qualities),
                           position+1,query_name, diploid_strand, flank)
                else:
                # -------------------
                # CASE 2: INSERTION
                # -------------------
                    insertion_context = read_seq[read_pos+1:read_pos + indel_length +1] # as observed

                    #homopolymer check
                    homopolymer = homopolymer_check(read_seq[read_pos:read_pos + indel_length +2])
                    if homopolymer == True: # filter out homopolymers
                        continue

                    # Counted against the opportunities of the gap between the anchor and the base after.
                    gap = canonical_context(read_seq[before] + read_seq[after])
                    row = (contig, read_pos+1,read_pos + indel_length +1,read_seq[read_pos],read_seq[read_pos:read_pos + indel_length +1],
                           "insertion", seen, insertion_context, canonical_context(insertion_context), gap, len(read_seq),
                           *read_tag_columns(alignment),sm,sx,sa,sm_pre,sx_pre,ToreadEnd,min(site.qualities),
                           position+1,query_name, diploid_strand, flank)

                print(*row, *entry.key, *strata.of_site(site), sep="\t", file=calls_out)

                # The read, for remapping to the haploid assembly
                if fastq_out is not None and query_name not in fastq_reads:
                    fastq_reads.add(query_name)
                    qual = ''.join(chr(x + 33) for x in alignment.get_forward_qualities())
                    fastq_out.write(f"@{query_name}\n{alignment.get_forward_sequence()}\n+\n{qual}\n")

    fasta.close()

    if opportunities_out is not None:
        print("contig", "indel_type", "size", "context", *strata.columns, "count", sep="\t", file=opportunities_out) #header for output
        for key, count in data.items():
            print(contig, *key, count, sep="\t", file=opportunities_out)
    return


def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    complement_map = str.maketrans("ATGC", "TACG")
    complement = sequence.translate(complement_map)
    return complement[::-1]

# ---------------------------------------------------------------------------
# Tandem-repeat InDel pipeline helpers
# ---------------------------------------------------------------------------

def load_tracts(bedfile, window_start, window_end, reach):
    """
    A model's tracts (an M0, M1 or M2 BED, 11 columns with a header) for a
    scan of [window_start, window_end): the tracts that begin in it, which
    are the window's, and the tracts a call in them can be counted in instead,
    those within `reach` bases. Returns (the window's tracts, an IntervalTree
    of all of them), each tract's data a dict of canonical_motif, k, copies,
    trf_purity, model_purity and model.

    Workflow 2 leaves one tract per base per model; tracts that overlap raise
    a ValueError, since calls and opportunities would then need a rule for
    choosing between them (R8).
    """
    rows = []
    with open(bedfile) as handle:
        next(handle)   # header
        for line in handle:
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            start, end = int(cols[1]), int(cols[2])
            if end > start:
                rows.append((start, end, cols))
    in_window = [row for row in rows if window_start <= row[0] < window_end]
    reach_end = max([window_end] + [end for _, end, _ in in_window]) + reach
    tree = IntervalTree()
    for start, end, cols in sorted(rows):
        if end <= window_start - reach or start >= reach_end:
            continue
        if tree.overlap(start, end):
            raise ValueError(f"{bedfile}: tracts overlap at {start}-{end}; workflow 4 needs one tract per base per model")
        tree[start:end] = {
            'canonical_motif': cols[3],
            'k':               int(cols[4]),
            'copies':          float(cols[5]),
            'trf_purity':      float(cols[6]),
            'model_purity':    float(cols[7]),
            'model':           cols[8],
        }
    window_tracts = sorted(iv for iv in tree if window_start <= iv.begin < window_end)
    return window_tracts, tree


def indel_class(indel_string, motif, k):
    """
    A repeat InDel's class (Q23): "slippage" if it is whole copies of the
    tract's motif, in any rotation, on either strand, else "other". Case is
    ignored: assemblies can be soft-masked.
    """
    copies, rest = divmod(len(indel_string), k)
    if rest or not copies:
        return "other"
    same = canonical_repeat_indel(indel_string.upper()) == canonical_repeat_indel(motif.upper() * copies)
    return "slippage" if same else "other"


def assigned_tract(tree, anchor, size, reach):
    """
    The tract an InDel of `size` bases anchored at `anchor` is counted in, as
    an interval of `tree`, or None.

    pbmm2 places an InDel that can slide at its leftmost position, so
    slippage in a tract is anchored on the base before it, or before a
    partial unit in front of it (R1), and slippage of a tract right after
    another on the first tract's last base (R2). So the InDel is counted in
    the tract of the first base after the anchor, up to `size` and at most
    `reach` bases on, that lies in a tract, or else the anchor's: the bases
    workflow 3 tests against its mask (L1).
    """
    for base in [*range(anchor + 1, anchor + min(size, reach) + 1), anchor]:
        here = tree[base]
        if here:
            return next(iter(here))
    return None


def tandem_repeat_scan(tr_beds, contig, bamfile, fastafile,
                       average_coverage_file, tissue_file,
                       window_start=None, window_end=None,
                       calls_out=None, opportunities_out=None, fastq_out=None,
                       *, coverage_band, key_parameters, flank_width, reach, min_flank,
                       read_bins_file="NA"):
    """
    Call InDels and count opportunities inside the tandem repeat tracts of
    several repeat models (M0, M1, M2) in one pass over the BAM.

    tr_beds maps a model name to its BED file. calls_out and opportunities_out
    map a model name to an open file for its calls and for its opportunity
    counts; a model left out of one gets no such output. fastq_out, when given,
    gets every read carrying a call in any model, once, as forward-strand FASTQ
    for haploid remapping.

    The unit is one tract in one read that spans it and at least min_flank
    bases on each side (Q22): a read ending inside a tract or its flank can't
    show a change in the tract's length, so it counts for neither.

    Opportunities: one per tract per spanning read, counted per motif, k, the
    tract keys (the tract's copies and purity, the read's flank on each side)
    and the read and site keys of strata.py.

    Calls: every InDel of a spanning read that is counted in the tract
    (assigned_tract(): anchored in it, or up to `reach` bases before it), of
    any size, with its class (indel_class(): slippage or other) and size, the
    tract's columns and the same keys as its opportunity. No homopolymer rule:
    k = 1 slippage is the signal.

    A tract's call and opportunity on a read are tested on the same bases
    (R7): kmerQ and subread support over the read's bases from two before
    the tract to two after it. The site's distance to the read end is the
    nearer flank; its nearest event is the read's nearest InDel or mismatch
    other than the tract's own InDels.

    Only tracts whose first base has a depth within coverage_band, (low, high)
    multiples of the sample's mean coverage, are used. key_parameters are the
    `strata` block of config.yaml. Each call carries flank_width bases of the
    assembly on each side of its anchor. read_bins_file is the compact table
    of the sample's read bins (read_bins.py, ticket 44), or "NA".

    With window_start/window_end only the tracts that begin in [window_start,
    window_end) are used; each tract belongs to one window, so windows never
    report the same call twice and their opportunity counts add up to the
    whole contig. Without them the whole contig is scanned.
    """
    calls_out = calls_out or {}
    opportunities_out = opportunities_out or {}

    with open(average_coverage_file, 'r') as f:
        average_coverage_of_species = float(f.read().strip())
    low_depth, high_depth = (average_coverage_of_species * factor for factor in coverage_band)

    fasta = pysam.FastaFile(fastafile)
    window_start = 0 if window_start is None else window_start
    window_end = fasta.get_reference_length(contig) if window_end is None else window_end
    tracts = {model: load_tracts(bedfile, window_start, window_end, reach) for model, bedfile in tr_beds.items()}

    d = load_tissue_classification(tissue_file)
    classified = load_read_bins(read_bins_file)
    strata = Strata(**key_parameters)

    # counts[model][(canonical_motif, k, *tract keys, *read keys, *site keys)] = count
    counts = {model: defaultdict(int) for model in opportunities_out}

    for out in calls_out.values():
        print("#contig", "start", "end", "ref", "alt",
              "indel_type", "unique", "indel_string", "canonical_indel", "read_length",
              "np", "ec", "rq", "fn", "rn", "mapQ", "sm", "sx", "sa", "sm_pre", "sx_pre", "ToreadEnd", "kmerQ",
              "position", "readID", "diploid_direction", "flank_sequence",
              "canonical_motif", "k", "copies", "trf_purity", "model_purity", "model", "indel_class", "size",
              *strata.tract_columns, *strata.columns,
              sep="\t", file=out)

    # The reads that can span a window's tract, or carry an InDel counted in
    # one, or cover a base where such an InDel is anchored.
    fetch_start = max(0, window_start - reach)
    fetch_end = max([window_end] + [iv.end for window_tracts, _ in tracts.values() for iv in window_tracts])

    def alignments(bam):
        for alignment in bam.fetch(contig, fetch_start, fetch_end):
            if not alignment.is_unmapped and alignment.cigartuples:
                yield alignment

    # First pass: each window tract's depth at its first base, as pysam's
    # pileup counts it (every read whose alignment spans the base), and, per
    # base, the reads with an InDel anchored there or a deletion over it, for
    # the calls' unique flag.
    starts, ends = [], []
    at_base = defaultdict(int)
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for alignment in alignments(bam):
            starts.append(alignment.reference_start)
            ends.append(alignment.reference_end)
            for indel in indels(alignment):
                at_base[indel.anchor] += 1
                if indel.deletion:
                    for base in range(indel.anchor + 1, indel.anchor + 1 + indel.length):
                        at_base[base] += 1
    starts.sort()
    ends.sort()

    def depth(base):
        return bisect_right(starts, base) - bisect_right(ends, base)

    usable = {model: {iv.begin for iv in window_tracts if low_depth <= depth(iv.begin) <= high_depth}
              for model, (window_tracts, _) in tracts.items()}

    fastq_reads = set()

    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for alignment in alignments(bam):
            read_start, read_end = alignment.reference_start, alignment.reference_end
            spanned = {}
            inner_start, inner_end = read_start + min_flank, read_end - min_flank
            for model, (_, tree) in tracts.items():
                if inner_start >= inner_end:
                    break
                # The tracts the read spans with min_flank bases on each side.
                found = sorted(iv for iv in tree.envelop(inner_start, inner_end) if iv.begin in usable[model])
                if found:
                    spanned[model] = found
            if not spanned:
                continue

            entry = read_entry(alignment, d, classified, strata, fasta, contig)
            blocks = Blocks(alignment)
            read_indels = indels(alignment)
            query_name = alignment.query_name
            read_seq = entry.sequence
            diploid_strand = "-" if alignment.is_reverse else "+"
            carries_call = False

            for model, found in spanned.items():
                tree = tracts[model][1]
                # The read's InDels by the tract they are counted in.
                own = defaultdict(list)
                for indel in read_indels:
                    tract = assigned_tract(tree, indel.anchor, indel.length, reach)
                    if tract is not None:
                        own[tract.begin].append(indel)

                for iv in found:
                    tr = iv.data
                    k = tr['k']
                    left, right = iv.begin - read_start, read_end - iv.end
                    first, last = blocks.span(iv.begin - 2, iv.end + 2)
                    support = None if entry.support is None else entry.support[first:last]
                    nearest = nearest_outside(entry.events, [indel.anchor for indel in own[iv.begin]],
                                              iv.begin, iv.end)
                    site = Site(entry.qualities[first:last], support, min(left, right), nearest)
                    keys = (strata.of_tract(Tract(tr['copies'], tr['model_purity'], left, right))
                            + entry.key + strata.of_site(site))
                    if model in counts:
                        counts[model][(tr['canonical_motif'], k) + keys] += 1

                    if model not in calls_out:
                        continue
                    for indel in own[iv.begin]:
                        position, read_pos, indel_length = indel.anchor, indel.read_pos, indel.length
                        seen = "unique" if at_base[position] <= 1 else "non-unique"
                        ToreadEnd = min(position - read_start, read_end - position)
                        flank = fasta.fetch(contig, max(0, position - flank_width), position + flank_width + 1)

                        # The subread tags in the alignment's order (S1).
                        sm = sx = sa = sm_pre = sx_pre = "NA"
                        if entry.tags is not None:
                            sm_all, sx_all, sa_all = entry.tags
                            sm, sx = sm_all[read_pos + 1], sx_all[read_pos + 1]   # the base after the anchor
                            sa, sm_pre, sx_pre = sa_all[read_pos], sm_all[read_pos], sx_all[read_pos]   # the anchor

                        if indel.deletion:
                            hcheck = fasta.fetch(contig, position, position + indel_length + 2)
                            indel_string = hcheck[1:-1]
                            row = (contig, read_pos + 1, read_pos + 1, hcheck[0:-1], hcheck[0],
                                   "deletion", seen, indel_string, canonical_repeat_indel(indel_string))
                        else:
                            indel_string = read_seq[read_pos + 1:read_pos + indel_length + 1]
                            row = (contig, read_pos + 1, read_pos + indel_length + 1,
                                   read_seq[read_pos], read_seq[read_pos:read_pos + indel_length + 1],
                                   "insertion", seen, indel_string, canonical_repeat_indel(indel_string))
                        print(*row, len(read_seq),
                              *read_tag_columns(alignment), sm, sx, sa, sm_pre, sx_pre, ToreadEnd, min(site.qualities),
                              position + 1, query_name, diploid_strand, flank,
                              tr['canonical_motif'], tr['k'], tr['copies'],
                              tr['trf_purity'], tr['model_purity'], tr['model'],
                              indel_class(indel_string, tr['canonical_motif'], k), indel_length,
                              *keys, sep="\t", file=calls_out[model])
                        carries_call = True

            # The read, for remapping to the haploid assembly
            if carries_call and fastq_out is not None and query_name not in fastq_reads:
                fastq_reads.add(query_name)
                qual = ''.join(chr(x + 33) for x in alignment.get_forward_qualities())
                fastq_out.write(f"@{query_name}\n{alignment.get_forward_sequence()}\n+\n{qual}\n")

    for model, out in opportunities_out.items():
        print("contig", "canonical_motif", "k", *strata.tract_columns, *strata.columns, "count",
              sep="\t", file=out)
        for key, count in counts[model].items():
            print(contig, *key, count, sep="\t", file=out)

    fasta.close()
    return


def resolve_haplotypes(contig, window_start, window_end, bamfile, fastafile, average_coverage_file,
                       *, min_base_quality, min_reads):
    # A region is the span of the reads that share a mismatch with the
    # assembly: at least min_reads reads with the same base, each at base
    # quality min_base_quality or more. A region that overlaps or touches the
    # one being built joins it, which becomes their union.
    # Only positions in [window_start, window_end) are scanned. Regions are
    # reported with the full span of their reads, so they can reach outside
    # the window; the workflow merges regions from neighbouring windows.
    fasta = pysam.FastaFile(fastafile)
    #average_coverage_of_species = extract_avg_cov(average_coverage_file)
    #print("works?")
    min_mapping_qual = 0
    min_base_qual = 0
    # The region being built, None until the first; one at the contig's first
    # base starts at 0.
    previous_start = None
    previous_end = None
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        #print("contig","start","end", sep= "\t")

        for pileupcolumn in bam.pileup(contig=contig,
                                       start=window_start,
                                       stop=window_end,
                                       fastafile=fasta,
                                       truncate=True,
                                       stepper="nofilter",
                                       min_base_quality=min_base_qual,#quite light
                                       min_mapping_quality=min_mapping_qual#conservative
        ):

            position = pileupcolumn.reference_pos
            #reference_base = fasta.fetch(contig, position, position + 1).upper()
            #print(position+1)
            # if not (average_coverage_of_species/2 <= pileupcolumn.nsegments <= average_coverage_of_species*2): #tjek on avergae depth of position
            #     continue
            bases = pileupcolumn.get_query_sequences(mark_matches=True, mark_ends=True,add_indels=False)
            bases = [base.upper() for base in bases]
            #print(bases)
            if len(set(bases)) != 1:
                reference_base = fasta.fetch(contig, position, position + 1).upper()
                #print(bases, reference_base, position, sep= "\t") 
                mismatch_positions = {}
                for pileupread in pileupcolumn.pileups:
                    query_position = pileupread.query_position
                    if query_position is None:
                        continue  # Skip if query position is None #potentiel indel
                    
                    read_base = pileupread.alignment.query_sequence[query_position].upper()
                    baseq = pileupread.alignment.query_qualities[query_position]
                    if baseq < min_base_quality:
                        continue

                    if read_base != reference_base: # if the read base is not the same as the reference base
                        readstart = pileupread.alignment.reference_start # keep zero indexing
                        readend = pileupread.alignment.reference_end
                        if read_base not in mismatch_positions:
                            mismatch_positions[read_base] = [[],[]]
                        mismatch_positions[read_base][0].append(readstart)
                        mismatch_positions[read_base][1].append(readend)
                        #print(bases, reference_base, position+1, sep= "\t") 
                # if position+1 == 172309:
                #     print(mismatch_positions)
                for base, intervals in mismatch_positions.items():
                    if len(intervals[0]) >= min_reads: # number of identical bases i need to see to say there is a unresolved haplotype
                        #print(intervals)
                        start = min(intervals[0])
                        end = max(intervals[1])
                        #print(intervals,position+1,start,end)
                        if previous_start is None:
                            previous_start = start
                            previous_end = end
                        elif start <= previous_end: # overlapping or touching: the union, never shrinking
                            previous_start = min(previous_start, start)
                            previous_end = max(previous_end, end)
                        else: # new region
                            print(contig, previous_start, previous_end, sep="\t")
                            previous_start = start
                            previous_end = end
    if previous_start is not None:
        print(contig, previous_start, previous_end, sep="\t")