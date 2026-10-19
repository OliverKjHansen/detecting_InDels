import pysam
import csv
from collections import defaultdict, namedtuple
from intervaltree import IntervalTree

from read_tags import describe, events, nearest_other_event, subread_support, subread_tags, tag
from strata import Site, Strata

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

def load_filters_with_tree(bedfile, window_start=None, window_end=None): #for intervals interval trees are quicker then the previous sweep-line approach Ol(og(n)) = O(n)
    """
    Reads bed file and creates an IntervalTree for O(log n) lookups.
    Much faster than sweep-line when checking if a position is in any interval.
    With window_start/window_end only the intervals overlapping
    [window_start, window_end) are kept, which are all the intervals a position
    in that window can fall in.

    Returns:
        IntervalTree: Tree for fast interval overlap queries
    """
    tree = IntervalTree()

    with open(bedfile, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            start, end = int(columns[1]), int(columns[2])
            if window_start is not None and end <= window_start:
                continue
            if window_end is not None and start >= window_end:
                continue
            # Create interval [start, end) - end is exclusive
            tree[start:end] = True
    
    return tree


def load_read_classification(classification_file):
    """
    readID -> germline/somatic bin, from read_classification.sh's table.
    Samples without methylation scores pass "NA": every read is NA.
    """
    if classification_file == "NA":
        return {}
    with open(classification_file) as f:
        header = f.readline().rstrip("\n").split("\t")
        read, label = header.index("readID"), header.index("classification_bin")
        return {fields[read]: fields[label] for fields in (line.rstrip("\n").split("\t") for line in f)}


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

def subread_array(array,readposition):
    """
    takes the array for subreads coverage and find the position of the read in the array
    this script might change with time when we better understand the sa tag
    """
    arrayposition = 0
    for i in range(0,len(array),2):
        previous_position = arrayposition
        arrayposition += array[i]
        coverage = array[i + 1]
        if previous_position <= readposition < arrayposition:
            #return [previous_position,readposition,arrayposition,coverage]
            return coverage

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
                       *, coverage_band, indel_sizes, key_parameters, flank_width, classification_file="NA"):
    """
    Call InDels of the sizes in indel_sizes, (smallest, largest) in bp, and
    count opportunities outside the mask in one pass over the BAM, and write
    the reads carrying a call as FASTQ.

    calls_out, opportunities_out and fastq_out are open files; whichever is
    None is not produced.

    Only positions whose depth is within coverage_band, (low, high) multiples
    of the sample's mean coverage, are scanned. key_parameters are the
    `strata` block of config.yaml. Each call carries flank_width bases of the
    assembly on each side of its anchor. classification_file is
    read_classification.sh's table of the sample's reads, or "NA".

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
    interval_tree = load_filters_with_tree(filters, window_start, window_end) # IntervalTree
    d = load_tissue_classification(tissue_file)
    classified = load_read_classification(classification_file)
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

    # Canonical k-mer and homopolymer check of a run (k-mer plus its flanking
    # bases), and canonical pair of bases around a gap, computed once each.
    runs = {}
    gaps = {}

    fastq_reads = set()

    with pysam.AlignmentFile(bamfile, "rb") as bam: # print the header
        for pileupcolumn in bam.pileup(contig=contig, #pileup
                                       start=window_start,
                                       stop=window_end,
                                       fastafile=fasta,
                                       truncate=True,
                                       stepper="nofilter",
                                       min_base_quality=min_base_qual,#quite light
                                       min_mapping_quality=min_mapping_qual,#conservative
                                       ):

            # quality fitler on reference position before determining indel type
            if not (low_depth <= pileupcolumn.nsegments <= high_depth):
                continue

            position = pileupcolumn.reference_pos # right. the position in the reference

            if interval_tree[position]: #saved me many lines of code
                continue

            if position >= next_eviction:
                reads = {key: entry for key, entry in reads.items() if entry.end > position}
                next_eviction = position + 10_000

            pileups = pileupcolumn.pileups

            if calls_out is not None:
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
                if opportunities_out is not None and read_pos >= 2:
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
                            data[("insertion", indel_length, context) + gap] += 1

                    for indel_length in sizes: #indel lengths here means kmer length
                        if not has_flanks(read_pos - 1, read_pos + indel_length, read_length):
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
                if calls_out is None:
                    continue

                indel = pileupread.indel

                if not smallest <= abs(indel) <= largest:
                    continue

                indel_length = abs(indel)
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

def load_tr_bed_with_tree(bedfile, window_start=None, window_end=None):
    """
    Load an M0 or M1 BED file (11-column, with header) into an IntervalTree.

    Each interval's data payload is a dict with keys:
        canonical_motif, k, copies, trf_purity, model_purity, model

    The first line (header) is skipped automatically.
    Intervals where end <= start are silently skipped.
    With window_start/window_end only the intervals overlapping
    [window_start, window_end) are kept, which are all the intervals a position
    in that window can fall in.
    """
    tree = IntervalTree()
    with open(bedfile, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                continue  # skip header
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            start, end = int(cols[1]), int(cols[2])
            if end <= start:
                continue
            if window_start is not None and end <= window_start:
                continue
            if window_end is not None and start >= window_end:
                continue
            tree[start:end] = {
                'canonical_motif': cols[3],
                'k':               int(cols[4]),
                'copies':          float(cols[5]),
                'trf_purity':      float(cols[6]),
                'model_purity':    float(cols[7]),
                'model':           cols[8],
            }
    return tree


def tandem_repeat_scan(tr_beds, contig, bamfile, fastafile,
                       average_coverage_file, tissue_file,
                       window_start=None, window_end=None,
                       calls_out=None, opportunities_out=None, fastq_out=None,
                       *, coverage_band, key_parameters, flank_width, classification_file="NA"):
    """
    Call InDels and count opportunities inside tandem repeat tracts for several
    repeat models (M0, M1, M2) in one pass over the BAM.

    tr_beds maps a model name to its BED file. calls_out and opportunities_out
    map a model name to an open file for its calls and for its opportunity
    counts; a model left out of one gets no such output. fastq_out, when given,
    gets every read carrying a call in any model, once, as forward-strand FASTQ
    for haploid remapping.

    For each model the output is exactly what that model gives on its own; the
    pileup, the coverage filter and the work per read are shared between the
    models.

    Calls mirror nonrepetitive_scan() but:
    - Use the TR BED as an *inclusion* filter (positions must be inside a TR interval).
    - InDel size filter: abs(indel) % k == 0 (any size, not just 1-3 bp).
    - Homopolymer filter is NOT applied (k=1 homopolymer InDels are the signal).
    - Extra output columns: canonical_motif, k, copies, trf_purity, model_purity, model.
    - When a position overlaps multiple TR intervals matching the InDel size,
      the interval with the largest k is used (most specific period).

    Opportunities mirror nonrepetitive_scan(): reads passing the quality
    filters at each position inside a TR interval, counted per (canonical_motif,
    k, tissue, mapQ, kmerQ). The kmerQ window uses the period k as the InDel
    length proxy. When a position overlaps multiple TR intervals, the one with
    the smallest k is used.

    Only positions whose depth is within coverage_band, (low, high) multiples
    of the sample's mean coverage, are scanned. key_parameters are the
    `strata` block of config.yaml. Each call carries flank_width bases of the
    assembly on each side of its anchor. classification_file is
    read_classification.sh's table of the sample's reads, or "NA".

    With window_start/window_end, only positions in [window_start, window_end)
    are scanned; each position belongs to one window, so windows never report
    the same call twice and their opportunity counts add up to the whole
    contig. Without them the whole contig is scanned.
    """
    calls_out = calls_out or {}
    opportunities_out = opportunities_out or {}

    with open(average_coverage_file, 'r') as f:
        average_coverage_of_species = float(f.read().strip())
    low_depth, high_depth = (average_coverage_of_species * factor for factor in coverage_band)

    fasta = pysam.FastaFile(fastafile)
    trees = {model: load_tr_bed_with_tree(bedfile, window_start, window_end)
             for model, bedfile in tr_beds.items()}

    d = load_tissue_classification(tissue_file)
    classified = load_read_classification(classification_file)

    min_mapping_qual = 0
    min_base_qual = 0

    strata = Strata(**key_parameters)

    # counts[model][(canonical_motif, k, *read keys, *site keys)] = count
    counts = {model: defaultdict(int) for model in opportunities_out}

    for out in calls_out.values():
        print("#contig", "start", "end", "ref", "alt",
              "indel_type", "unique", "indel_string", "canonical_indel", "read_length",
              "np", "ec", "rq", "fn", "rn", "mapQ", "sm", "sx", "sa", "sm_pre", "sx_pre", "ToreadEnd", "kmerQ",
              "position", "readID", "diploid_direction", "flank_sequence",
              "canonical_motif", "k", "copies", "trf_purity", "model_purity", "model",
              *strata.columns,
              sep="\t", file=out)

    # Each read's ReadEntry, built once rather than at every tract position it
    # covers: copying a whole HiFi read's qualities at each position is most of
    # the cost of counting opportunities. Keyed by read, start and flag, with
    # the read's end so reads the scan has passed can be dropped every 10 kb.
    reads = {}
    next_eviction = 0

    def entry_of(alignment):
        key = (alignment.query_name, alignment.reference_start, alignment.flag)
        entry = reads.get(key)
        if entry is None:
            entry = reads[key] = read_entry(alignment, d, classified, strata, fasta, contig)
        return entry

    fastq_reads = set()

    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for pileupcolumn in bam.pileup(contig=contig,
                                       start=window_start, stop=window_end,
                                       region=None, end=None,
                                       fastafile=fasta, truncate=True,
                                       stepper="nofilter",
                                       min_base_quality=min_base_qual,
                                       min_mapping_quality=min_mapping_qual):

            if not (low_depth <= pileupcolumn.nsegments <= high_depth):
                continue

            position = pileupcolumn.reference_pos

            # Only process positions INSIDE a TR interval, per model
            hits = {}
            for model, tree in trees.items():
                model_hits = tree[position]
                if model_hits:
                    hits[model] = model_hits
            if not hits:
                continue

            if position >= next_eviction:
                reads = {key: entry for key, entry in reads.items() if entry.end > position}
                next_eviction = position + 10_000

            pileups = pileupcolumn.pileups

            call_models = [model for model in hits if model in calls_out]
            if call_models:
                indel_count = sum(
                    1 for pr in pileups
                    if (pr.indel != 0) or pr.is_del
                )
                seen = "unique" if indel_count <= 1 else "non-unique"

            # Opportunities use the interval with the smallest k at this position
            opportunity_tracts = [(model, min(hits[model], key=lambda iv: iv.data['k']).data)
                                  for model in hits if model in opportunities_out]

            for pileupread in pileups:
                if pileupread.query_position is None:
                    continue

                read_pos   = pileupread.query_position
                alignment  = pileupread.alignment
                query_name = alignment.query_name

                # ---- opportunities ----
                if opportunity_tracts and 1 <= read_pos <= alignment.query_length - 2:
                    entry = entry_of(alignment)
                    base_qualities, support = entry.qualities, entry.support
                    # An opportunity's anchor is the base before its window, as a call's is the base before its gap.
                    anchor = position - 1
                    to_read_end = min(anchor - alignment.reference_start, entry.end - anchor)
                    nearest = nearest_other_event(entry.events, anchor)

                    for model, tr in opportunity_tracts:
                        k = tr['k']
                        # kmerQ window: 2 bases flanking a k-length window
                        kmerq = (base_qualities[read_pos - 1:read_pos + 1]
                                 + base_qualities[read_pos + k - 1:read_pos + k + 1])
                        if not kmerq:
                            continue
                        tested_support = None if support is None else (
                            support[read_pos - 1:read_pos + 1] + support[read_pos + k - 1:read_pos + k + 1])
                        counts[model][(tr['canonical_motif'], k) + entry.key
                                      + strata.of_site(Site(kmerq, tested_support, to_read_end, nearest))] += 1

                # ---- calls ----
                if not call_models:
                    continue

                indel = pileupread.indel
                if indel == 0:
                    continue

                indel_length = abs(indel)

                # Per model, pick the matching interval with the largest k
                # (most specific period for this InDel size)
                call_tracts = []
                for model in call_models:
                    matching = [iv for iv in hits[model]
                                if indel_length % iv.data['k'] == 0]
                    if matching:
                        call_tracts.append((model, max(matching, key=lambda iv: iv.data['k']).data))
                if not call_tracts:
                    continue

                sm = sx = sa = sm_pre = sx_pre = "NA"
                read_pos_after = pileupread.query_position_or_next
                entry          = entry_of(alignment)
                readstart      = alignment.reference_start
                ToreadEnd      = min(position - readstart, entry.end - position)
                nearest        = nearest_other_event(entry.events, position)
                flank          = fasta.fetch(contig, max(0, position - flank_width), position + flank_width + 1)

                read_seq       = entry.sequence
                base_qualities = entry.qualities

                # The subread tags in the alignment's order (S1).
                tags = entry.tags
                call_support = entry.support
                if tags is not None:
                    sm_all, sx_all, sa_all = tags
                    sm, sx = sm_all[read_pos + 1], sx_all[read_pos + 1]   # the base after the anchor
                    sa, sm_pre, sx_pre = sa_all[read_pos], sm_all[read_pos], sx_all[read_pos]   # the anchor

                diploid_strand = "-" if alignment.is_reverse else "+"

                if indel < 0:
                    # DELETION
                    hcheck = fasta.fetch(contig, position,
                                         position + indel_length + 2)
                    deletion_context = hcheck[1:-1]
                    def tested(values):   # the bases tested, for qualities and subread support alike
                        return (values[max(0, read_pos - 1):read_pos + 1]
                                + values[read_pos_after + 1:read_pos_after + 3])
                    del_q = tested(base_qualities)
                    # Quality 0, and no support, when no base is left to test.
                    site = Site(del_q or (0,), (call_support and tested(call_support)) or None, ToreadEnd, nearest)

                    row = (contig,
                           read_pos + 1, read_pos + 1,
                           hcheck[0:-1], hcheck[0],
                           "deletion", seen, deletion_context, canonical_repeat_indel(deletion_context),
                           len(read_seq),
                           *read_tag_columns(alignment), sm, sx, sa, sm_pre, sx_pre, ToreadEnd, min(site.qualities),
                           position + 1, query_name, diploid_strand, flank)
                else:
                    # INSERTION
                    insertion_context = read_seq[read_pos + 1:read_pos + indel_length + 1]
                    def tested(values):   # the bases tested, for qualities and subread support alike
                        return (values[max(0, read_pos):read_pos + 2]
                                + values[max(0, read_pos) + indel_length:read_pos + indel_length + 2])
                    ins_q = tested(base_qualities)
                    # Quality 0, and no support, when no base is left to test.
                    site = Site(ins_q or (0,), (call_support and tested(call_support)) or None, ToreadEnd, nearest)

                    row = (contig,
                           read_pos + 1, read_pos + indel_length + 1,
                           read_seq[read_pos],
                           read_seq[read_pos:read_pos + indel_length + 1],
                           "insertion", seen, insertion_context, canonical_repeat_indel(insertion_context),
                           len(read_seq),
                           *read_tag_columns(alignment), sm, sx, sa, sm_pre, sx_pre, ToreadEnd, min(site.qualities),
                           position + 1, query_name, diploid_strand, flank)

                keys = entry.key + strata.of_site(site)
                for model, tr in call_tracts:
                    print(*row,
                          tr['canonical_motif'], tr['k'], tr['copies'],
                          tr['trf_purity'], tr['model_purity'], tr['model'],
                          *keys,
                          sep="\t", file=calls_out[model])

                # The read, for remapping to the haploid assembly
                if fastq_out is not None and query_name not in fastq_reads:
                    fastq_reads.add(query_name)
                    qual = ''.join(chr(x + 33) for x in alignment.get_forward_qualities())
                    fastq_out.write(f"@{query_name}\n{alignment.get_forward_sequence()}\n+\n{qual}\n")

    for model, out in opportunities_out.items():
        print("contig", "canonical_motif", "k", *strata.columns, "count",
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