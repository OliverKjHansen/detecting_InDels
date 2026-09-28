import pysam
import sys
import csv
from collections import defaultdict
from intervaltree import IntervalTree, Interval

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


def small_indels_dict(microsatfile, contig, bamfile, fastafile, average_coverage_file, tissue_file,
                      window_start=None, window_end=None):

    with open(average_coverage_file, 'r') as f: # load in the average coverage smaple
        average_coverage_of_species = float(f.read().strip())

    fasta = pysam.FastaFile(fastafile)
    intervals, metadata = load_intervals(microsatfile) # artifact from microsat script, but might be useful later
    
    #####
    d = load_tissue_classification(tissue_file)
    #print(d)

    #####
    min_mapping_qual = 0 # to have more contorl over pysam
    min_base_qual = 0 #--||---

    def counts_dict():
        return {"filter": 0, "nofilter": 0, "3mer_del":0}

    data = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(counts_dict))))
    # data = {
    # "insertion": {
    #   "unique" : {
    #     "AT": {
    #         "liver": {"filter": 10, "nofilter": 4, 3mer: 3},
    #         "brain": {"count_1": 7, "count_2": 2}
    #     }
    #   }
    # } 
    # "deletion": {}
    # }
    #
    # multiple reads small indels
    #
    # we'll iterate using a sweep-line over intervals to skip unnecessary positions from microsats, because they overlap. re-thinkt his a bit
    next_iv = 0                              # CHANGED: index of next interval to activate
    active = []                              # CHANGED: list of (start,end) intervals currently covering position

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
            if not (average_coverage_of_species/2 <= pileupcolumn.nsegments <= average_coverage_of_species*2):
                continue

            position = pileupcolumn.reference_pos # right. the position in the reference

            ##skipping position filters
            # activate any intervals that start at/before this position
            while next_iv < len(intervals) and intervals[next_iv][0] - 2 <= position:
                active.append(intervals[next_iv])
                next_iv += 1

            # drop intervals that ended before/at this position (end is exclusive)
            if active:
                active = [(s, e) for (s, e) in active if e > position]
                #print(f"DEBUG position {position+1}, active intervals: {active}")
            if active:
                continue
            # do nothing if we are in  microsat intervals
            # if active[0][0] <= position <= active[0][1]:
            #     continue

            indel_count = sum(1 for pileupread in pileupcolumn.pileups if (pileupread.indel != 0) or (pileupread.is_del == True))

            seen = "unique"
            if indel_count > 1:
                seen = "non-unique"
            
            # fitler on positions on reads
            for pileupread in pileupcolumn.pileups: # looking for reads without indels and have the tags we are interested i
                    if pileupread.query_position is None: # we are in an indel and the position is none
                        continue

                    indel = pileupread.indel

                    if (abs(indel) > 3) or (abs(indel) == 0): # only interested in 1,2,3 sizes
                        continue

                    np = ec = sm = sx = sa = sm_pre = sx_pre = "NA" # initialise tags
                    read_pos = pileupread.query_position
                    read_pos_after = pileupread.query_position_or_next
                    readstart = pileupread.alignment.reference_start
                    readend = pileupread.alignment.reference_end
                    ToreadEnd = min([(position-readstart),(readend-position)]) #tjek

                    # if ToreadEnd < 1000: # if we are too close to the end of the read
                    #     continue

                    # for (start, end) in active: # dont think i need this ????
                    #     #print(f"DEBUG loop start-end {start+1}-{end}, position {position+1}, read {pileupread.alignment.query_name}")# shift for left-alignment
                    #     if not ((start-2) <= position < end): # because indels are left aligned
                    #         continue
                    #     if position >= end:
                    #         continue

                    # microsat_buffer = 2 #artifact ?!?!?

                    # read fully spans interval with buffer?
                    # if not ((readstart < (start - microsat_buffer)) and (end < (readend - microsat_buffer))):
                    #     continue
                        # make easier

                    read_seq = pileupread.alignment.query_sequence
                    base_qualities = pileupread.alignment.query_qualities

                    #readoverlap = min(start - readstart, readend - end)
                    #print("working?", position, active, read_pos, len(base_qualities),pileupread.alignment.query_name,pileupread.indel)
                    if pileupread.alignment.has_tag("np"):
                        np = pileupread.alignment.get_tag("np") # full read
                        ec = pileupread.alignment.get_tag("ec")
                    if pileupread.alignment.has_tag("sa"): # er et tjek nok?
                        sm = pileupread.alignment.get_tag("sm")[read_pos+1] # first base in the insertion
                        sx = pileupread.alignment.get_tag("sx")[read_pos+1]
                        sa = pileupread.alignment.get_tag("sa")#[pileupread.query_position_or_next] # this is an array so I index later in the function subread_array()
                        sa = subread_array(sa,read_pos)
                        sm_pre = pileupread.alignment.get_tag("sm")[read_pos] #base just before the insertion                                
                        sx_pre = pileupread.alignment.get_tag("sx")[read_pos]

                    indel_length = abs(indel)
                    upstream,downstream = closest_indel(pileupread.alignment.cigartuples,read_pos,indel_length)

                    #Make a tjek that sees if the indel is unique or not
                    
                    if indel < 0:
                    # -------------------
                    # CASE 1: DELETIONS
                    # -------------------
                        hcheck = fasta.fetch(contig, position, position + indel_length+2) #for homopolymer check
                        deletion_context = hcheck[1:-1] #context of the deletion
                        #deletion_context = fasta.fetch(contig, position+1, position + 1 + indel_length)
                        lexi = lexicographically(deletion_context)

                        #homopolymer_check

                        homopolymer = homopolymer_check(hcheck)
                        if homopolymer == True: # filter out homopolymers
                            continue

                        tissue_type = d.get(pileupread.alignment.query_name, "NA")
                        data["deletion"][seen][lexi][tissue_type]["nofilter"] += 1

                        #quality check
                        #AA|ACG|TT
                        #00|xxx|00
                        deletion_qualities_4mer = base_qualities[max(0, read_pos-1):read_pos+1]+base_qualities[read_pos_after+1:read_pos_after+3] # think this is correct
                        kmerQ = min(deletion_qualities_4mer)

                        if (len(deletion_context) == 1) and (min(deletion_qualities_4mer[0:3])) == 93:
                            data["deletion"][seen][lexi][tissue_type]["3mer_del"] += 1

                        if kmerQ < 93:
                            continue
                
                        data["deletion"][seen][lexi][tissue_type]["filter"] += 1

                    if indel > 0:
                    # -------------------
                    # CASE 2: INSERTION
                    # -------------------
                        insertion_context = read_seq[read_pos+1:read_pos + indel_length +1] 
                        lexi = lexicographically(insertion_context)

                        #homopolymer check
                        homopolymer = homopolymer_check(read_seq[read_pos:read_pos + indel_length +2])
                        if homopolymer == True: # filter out homopolymers
                            continue

                        tissue_type = d.get(pileupread.alignment.query_name, "NA")
                        data["insertion"][seen][lexi][tissue_type]["nofilter"] += 1

                        #quality check
                        #AA|ACG|TT
                        #x0|0x0|0x
                        insertion_qualities_4mer = base_qualities[max(0, read_pos):read_pos+2]+base_qualities[max(0, read_pos)+indel_length:read_pos+indel_length+2] # will count the middle base twice in 1 bp, but it doesnt matter
                        kmerQ = min(insertion_qualities_4mer)

                        if kmerQ < 93:
                            continue
                        
                        data["insertion"][seen][lexi][tissue_type]["filter"] += 1

                        lexi = lexicographically(insertion_context) # add lowest substring to sat_type
    fasta.close()

    # print the data
    print("contig","indel_type","unique","indel_string","tissue","nofilter","3mer","filter", sep= "\t", file=sys.stdout) #header for output
    for indel_type, uniques in data.items():
        for unique, indel_strings in uniques.items():
            for indel_string, tissue in indel_strings.items():
                for tissue, counts in tissue.items():
                    print(contig,indel_type,unique,indel_string,tissue,counts['nofilter'],counts["3mer_del"],counts['filter'], sep="\t", file=sys.stdout)
    return

def lexicographically(pattern: str) -> str:
    """
    Find the lexicographically of string and reverse complement
    """
    complement_map = str.maketrans("ATGC", "TACG")
    complement = pattern.translate(complement_map)
    all_strings = [pattern,complement]
    return min(all_strings)

def canonical_context(sequence):
    """
    Group a context with the same context read from the other strand: the
    smaller of the sequence and its reverse complement.

    A GT deletion on one strand is an AC deletion on the other, so both are
    reported as AC. lexicographically() above complements without reversing,
    which pairs GT with CA instead; it is left as it is because only the
    disabled code paths still call it.
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

def closest_indel(cigartuple, position, indel):   

    """
    Finds the closest indel up- and downstream in a read from a given position or indel (insertions and deletions) for a specified contig in a BAM file

    Parameters:
        cigartuple of the read in which the indel is from 
        postion (int): The position (0-index) in the read of the position of interest
    
    Returns:
        distance to the closest downstream indel
        distance to the closests upstream indel
    """
    #somewhere this function beomse 1 indexed for 
    upstream_indel = None
    downstream_indel = None
    readposition = 0
    for ctuple in cigartuple:
        operation = ctuple[0]
        length = ctuple[1]
        if downstream_indel != None and readposition > position:
            break
        #print(readposition, position, length, operation, downstream_indel,upstream_indel, operation,length)
        #print(readposition, position)
        if operation in (1,2) and readposition < position:
            upstream_indel = position-readposition # so it is zero indexed to the end of the indel, if i want start i should +len. furthermore if i want to 1 index i should plus with 1 in the end
            if operation == 1:
                readposition += length
        if operation in (1,2) and readposition > position+indel: #somewhere here
            downstream_indel = readposition-position
            if operation == 1:
                readposition += length
        if operation not in (1,2) and downstream_indel == None:
                readposition += length
    return upstream_indel, downstream_indel

def load_intervals(microsatfile):

    """
    Reads microsatellite intervals into a sorted list.
    """
    intervals = []
    metadata = {}
    
    with open(microsatfile, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            start, end = int(columns[1]), int(columns[2])
            sat_purity, sat_type, sat_length, sample, species = columns[6],columns[14], len(columns[15]), columns[18], columns[19] # could add purity here # col 7(col 6)
            
            intervals.append((start, end))
            metadata[(start, end)] = (sat_purity, sat_type, sat_length, sample, species)

    return sorted(intervals), metadata


def load_filters(bedfile):
    """
    Reads bed file that contains filters ww use to mask out the genome and creates a sorted list.
    """
    intervals = []
    
    with open(bedfile, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            start, end = int(columns[1]), int(columns[2])
            #sat_purity, sat_type, sat_length, sample, species = columns[6],columns[14], len(columns[15]), columns[18], columns[19] # could add purity here # col 7(col 6)
            
            intervals.append((start, end))
            #metadata[(start, end)] = (sat_purity, sat_type, sat_length, sample, species)

    return sorted(intervals)#, metadata

def large_indels(contig, bamfile, fastafile, average_coverage_file, tissue_file, microsatfile,
                 window_start=None, window_end=None):

    with open(average_coverage_file, 'r') as f: # load in the average coverage smaple
        average_coverage_of_species = float(f.read().strip())

    fasta = pysam.FastaFile(fastafile)
    intervals, metadata = load_intervals(microsatfile) # artifact from microsat script, but might be useful later
    
    #####
    d = load_tissue_classification(tissue_file)
    #print(d)
    ### ADD mircosat filter in this analysis and in tmotif flag
    print("contig","position","ref","alt",
    "indeltype","indel_context",
    "tissue_type","read_position","read_length",
    "np","ec","sm","sx","sa","sm_pre","sx_pre","upstream","downstream",
    "2mer_before","2mer_after",
    "microsat_inteval", "microsat_type", sep="\t",file=sys.stdout)


    #####
    min_mapping_qual = 0 # to have more contorl over pysam
    min_base_qual = 0 #--||---

    next_iv = 0                             
    active = []                              

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

            microsat_interval = False
            microsat_type = "NA"

            if not (average_coverage_of_species/2 <= pileupcolumn.nsegments <= average_coverage_of_species*2):
                continue

            position = pileupcolumn.reference_pos # right. the position in the reference

            ##skipping position filters
            # activate any intervals that start at/before this position
            while next_iv < len(intervals) and intervals[next_iv][0] - 2 <= position:
                active.append(intervals[next_iv])
                next_iv += 1

            # drop intervals that ended before/at this position (end is exclusive)
            if active:
                active = [(s, e) for (s, e) in active if e > position]
                #print(f"DEBUG position {position+1}, active intervals: {active}")


            for (start, end) in active:
                        #print(f"DEBUG loop start-end {start+1}-{end}, position {position+1}, read {pileupread.alignment.query_name}")# shift for left-alignment
                #print(f"DEBUG loop start-end {start}-{end}, position {position}, read {pileupread.alignment.query_name}")# shift for left-alignment

                # if not ((start-2) <= position < end): # because indels are left aligned
                #      continue
                # if position >= end:
                #     continue

                # microsat_buffer = 2
                    # read fully spans interval with buffer?
                # if not ((readstart < (start - microsat_buffer)) and (end < (readend - microsat_buffer))):
                #     continue

            # if active:
                microsat_interval = True
                sat_purity, sat_type, sat_length, sample, species = metadata[(start, end)]
                microsat_type = find_lowest_substring(sat_type)

            for pileupread in pileupcolumn.pileups: # looking for reads without indels and have the tags we are interested i
                    if pileupread.query_position is None: # we are in an indel and the position is none
                        continue
                    indel = pileupread.indel

                    if (abs(indel) < 10): # only interested larger then 10 indels
                        continue

                    np = ec = sm = sx = sa = sm_pre = sx_pre = "NA" # initialise tags
                    read_pos = pileupread.query_position
                    read_pos_after = pileupread.query_position_or_next
                    readstart = pileupread.alignment.reference_start
                    readend = pileupread.alignment.reference_end
                   # ToreadEnd = min([(position-readstart),(readend-position)]) #tjek

                    read_seq = pileupread.alignment.query_sequence
                    base_qualities = pileupread.alignment.query_qualities

                    #readoverlap = min(start - readstart, readend - end)
                    # print("working?", position, active, read_pos, len(base_qualities),pileupread.alignment.query_name,pileupread.indel)
                    if pileupread.alignment.has_tag("np"):
                        np = pileupread.alignment.get_tag("np") # full read
                        ec = pileupread.alignment.get_tag("ec")
                    if pileupread.alignment.has_tag("sa"): # er et tjek nok?
                        sm = pileupread.alignment.get_tag("sm")[read_pos+1] # first base in the insertion
                        sx = pileupread.alignment.get_tag("sx")[read_pos+1]
                        sa = pileupread.alignment.get_tag("sa")#[pileupread.query_position_or_next] # this is an array so I index later in the function subread_array()
                        sa = subread_array(sa,read_pos)
                        sm_pre = pileupread.alignment.get_tag("sm")[read_pos] #base just before the insertion                                
                        sx_pre = pileupread.alignment.get_tag("sx")[read_pos]

                    indel_length = abs(indel)
                    upstream,downstream = closest_indel(pileupread.alignment.cigartuples,read_pos,indel_length)
                    base_qualities = pileupread.alignment.query_qualities

                    if indel < 0:
                    # -------------------
                    # CASE 1: DELETIONS
                    # -------------------
                        deletion_context = fasta.fetch(contig, position+1, position + 1 + indel_length) #the bases deleted
                        #lexi = lexicographically(insertion_context)

                        ref = read_seq[read_pos]+deletion_context
                        alt = read_seq[read_pos]

                        tissue_type = d.get(pileupread.alignment.query_name, "NA")

                        mer2_before = base_qualities[max(0, read_pos-1):read_pos+1]
                        mer2_after = base_qualities[read_pos_after+1:read_pos_after+3]

                        print(contig,position+1,ref,alt,"deletion",deletion_context,tissue_type,read_pos,len(read_seq),
                        np,ec,sm,sx,sa,sm_pre,sx_pre,upstream,downstream,
                        min(mer2_before),min(mer2_after),
                        microsat_interval,microsat_type, sep="\t",file=sys.stdout)


                    if indel > 0:
                    # -------------------
                    # CASE 2: INSERTION
                    # -------------------
                        insertion_context = read_seq[read_pos+1:read_pos + indel_length +1] 
                        #lexi = lexicographically(insertion_context)
                        ref = read_seq[read_pos]
                        alt = read_seq[read_pos]+insertion_context

                        tissue_type = d.get(pileupread.alignment.query_name, "NA")

                        mer2_before = base_qualities[max(0, read_pos):read_pos+2]
                        mer2_after = base_qualities[max(0, read_pos)+indel_length:read_pos+indel_length+2]

                        print(contig,position+1,ref,alt,"insertion",insertion_context,tissue_type,read_pos,len(read_seq),
                        np,ec,sm,sx,sa,sm_pre,sx_pre,upstream,downstream,
                        min(mer2_before),min(mer2_after),
                        microsat_interval,microsat_type, sep="\t",file=sys.stdout)

    fasta.close()
    return

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

def nonrepetitive_scan(filters, contig, bamfile, fastafile, average_coverage_file, tissue_file,
                       window_start=None, window_end=None,
                       calls_out=None, opportunities_out=None, fastq_out=None):
    """
    Call 1-3 bp InDels and count opportunities outside the mask in one pass
    over the BAM, and write the reads carrying a call as FASTQ.

    calls_out, opportunities_out and fastq_out are open files; whichever is
    None is not produced. small_indels() and indel_opportunities() are this
    function with only calls_out or only opportunities_out, on stdout.

    Calls: 1-3 bp InDels at positions outside the mask (the filters BED),
    except those whose sequence plus one flanking base on each side is a
    homopolymer (changes base at most once).

    Opportunities: every read position outside the mask, counted once per
    InDel size 1, 2 and 3 under the canonical k-mer of that size (the smaller
    of the k-mer and its reverse complement), per tissue, mapQ and kmerQ, with
    the same homopolymer rule applied to the k-mer and its flanking bases.

    fastq_out gets every read carrying a call, once, as forward-strand FASTQ
    for haploid remapping.

    With window_start/window_end only positions in [window_start, window_end)
    are scanned; each position belongs to one window, so windows never report
    the same call twice and their opportunity counts add up to the whole
    contig. Without them the whole contig is scanned.
    """
    with open(average_coverage_file, 'r') as f: # load in the average coverage smaple
        average_coverage_of_species = float(f.read().strip())

    fasta = pysam.FastaFile(fastafile)
    interval_tree = load_filters_with_tree(filters, window_start, window_end) # IntervalTree
    d = load_tissue_classification(tissue_file)

    min_mapping_qual = 0 # to have more contorl over pysam
    min_base_qual = 0 #--||---

    def counts_dict():
        return {str(i): 0 for i in range(94)}

    # data[context][tissue][mapQ][kmerQ] = count
    data = defaultdict(
    lambda: defaultdict(
        lambda: defaultdict(counts_dict)
        )
    )

    if calls_out is not None:
        print("#contig", "start", "end", "ref", "alt",
        "indel_type","unique","indel_string","canonical_indel","tissue","read_length",
        "np","ec","sm","sx","sa","sm_pre","sx_pre","ToreadEnd","kmerQ",
        "position","readID", "diploid_direction",
        sep= "\t", file=calls_out) #header for output

    # Sequence and base qualities of each read, fetched once rather than at
    # every position it covers: copying a whole HiFi read at each position is
    # most of the cost of counting opportunities. Keyed by read, start and
    # flag, with the read's end so reads the scan has passed can be dropped
    # every 10 kb.
    reads = {}
    next_eviction = 0

    # Canonical k-mer and homopolymer check of a run (k-mer plus its flanking
    # bases), computed once per distinct run.
    runs = {}

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
            if not (average_coverage_of_species/2 <= pileupcolumn.nsegments <= average_coverage_of_species*2):
                continue

            position = pileupcolumn.reference_pos # right. the position in the reference

            if interval_tree[position]: #saved me many lines of code
                continue

            if position >= next_eviction:
                reads = {key: value for key, value in reads.items() if value[0] > position}
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
                # we need at least one base on each side of the kmer
                if opportunities_out is not None and 1 <= read_pos <= alignment.query_length - 2:
                    key = (query_name, alignment.reference_start, alignment.flag)
                    cached = reads.get(key)
                    if cached is None:
                        cached = reads[key] = (alignment.reference_end, alignment.query_sequence, alignment.query_qualities)
                    _, read_seq, base_qualities = cached
                    mapQ = str(alignment.mapping_quality)
                    tissue_type = d.get(query_name, "NA")

                    for indel_length in (1, 2, 3): #indel lengths here means kmer length
                        run = read_seq[read_pos-1:read_pos+indel_length+1]   #run dependent on indel length # includes kmer + flanking bases
                        kmerq = base_qualities[read_pos-1:read_pos+1]+base_qualities[read_pos+indel_length-1:read_pos+indel_length+1] #look weird but it makes sure that the 3 indel does not inlude the central base to insure a 4merQ

                        #quality 3size
                        #NT|ACT|CN
                        #x0|0x0|0x

                        #quality 2size
                        #NTA|CT|CN
                        #xx0|00|0x

                        #quality 1size
                        #TA|C|TC
                        #x0|0|0x

                        info = runs.get(run)
                        if info is None:
                            #context of kmer without the flanking bases, grouped with the other strand
                            info = runs[run] = (canonical_context(run[1:-1]), homopolymer_check(run))
                        lexi, homopolymer = info
                        if homopolymer == True: # filter out homopolymers
                            continue

                        data[lexi][tissue_type][mapQ][str(min(kmerq))] += 1

                # ---- calls ----
                if calls_out is None:
                    continue

                indel = pileupread.indel

                if (abs(indel) > 3) or (abs(indel) == 0): # only interested in 1,2,3 sizes
                    continue

                np = ec = sm = sx = sa = sm_pre = sx_pre = "NA" # initialise tags
                read_pos_after = pileupread.query_position_or_next
                readstart = alignment.reference_start
                readend = alignment.reference_end
                ToreadEnd = min([(position-readstart),(readend-position)]) #tjek

                read_seq = alignment.query_sequence
                base_qualities = alignment.query_qualities

                if alignment.has_tag("np"):
                    np = alignment.get_tag("np") # full read
                    ec = alignment.get_tag("ec")
                if alignment.has_tag("sa"): # er et tjek nok?
                    sm = alignment.get_tag("sm")[read_pos+1] # first base in the insertion
                    sx = alignment.get_tag("sx")[read_pos+1]
                    sa = alignment.get_tag("sa") # this is an array so I index later in the function subread_array()
                    sa = subread_array(sa,read_pos)
                    sm_pre = alignment.get_tag("sm")[read_pos] #base just before the insertion
                    sx_pre = alignment.get_tag("sx")[read_pos]

                indel_length = abs(indel)
                diploid_strand = "-" if alignment.is_reverse == True else "+"

                if indel < 0:
                # -------------------
                # CASE 1: DELETIONS
                # -------------------
                    hcheck = fasta.fetch(contig, position, position + indel_length+2) #for homopolymer check
                    deletion_context = hcheck[1:-1] #the deletion context, as observed

                    homopolymer = homopolymer_check(hcheck)
                    if homopolymer == True: # filter out homopolymers
                        continue

                    tissue_type = d.get(query_name, "NA")

                    #quality check
                    #AA|ACG|TT
                    #00|xxx|00
                    deletion_qualities_4mer = base_qualities[max(0, read_pos-1):read_pos+1]+base_qualities[read_pos_after+1:read_pos_after+3] # need to dobblet check this
                    kmerQ = min(deletion_qualities_4mer)

                    row = (contig, read_pos+1,read_pos+1,hcheck[0:-1],hcheck[0], #start and end are the same for delettion, becaus its is in the read
                           "deletion", seen, deletion_context, canonical_context(deletion_context), tissue_type, len(read_seq),
                           np,ec,sm,sx,sa,sm_pre,sx_pre,ToreadEnd,kmerQ,
                           position+1,query_name, diploid_strand)
                else:
                # -------------------
                # CASE 2: INSERTION
                # -------------------
                    insertion_context = read_seq[read_pos+1:read_pos + indel_length +1] # as observed

                    #homopolymer check
                    homopolymer = homopolymer_check(read_seq[read_pos:read_pos + indel_length +2])
                    if homopolymer == True: # filter out homopolymers
                        continue

                    tissue_type = d.get(query_name, "NA")

                    #quality check
                    #AA|ACG|TT
                    #x0|0x0|0x
                    insertion_qualities_4mer = base_qualities[max(0, read_pos):read_pos+2]+base_qualities[max(0, read_pos)+indel_length:read_pos+indel_length+2] # will count the middle base twice in 1 bp, but it doesnt matter
                    kmerQ = min(insertion_qualities_4mer)

                    row = (contig, read_pos+1,read_pos + indel_length +1,read_seq[read_pos],read_seq[read_pos:read_pos + indel_length +1],
                           "insertion", seen, insertion_context, canonical_context(insertion_context), tissue_type, len(read_seq),
                           np,ec,sm,sx,sa,sm_pre,sx_pre,ToreadEnd,kmerQ,
                           position+1,query_name, diploid_strand)

                print(*row, sep="\t", file=calls_out)

                # The read, for remapping to the haploid assembly
                if fastq_out is not None and query_name not in fastq_reads:
                    fastq_reads.add(query_name)
                    qual = ''.join(chr(x + 33) for x in alignment.get_forward_qualities())
                    fastq_out.write(f"@{query_name}\n{alignment.get_forward_sequence()}\n+\n{qual}\n")

    fasta.close()

    if opportunities_out is not None:
        print("contig","context","tissue","mapQ","kmerQ","count", sep= "\t", file=opportunities_out) #header for output
        for context, tissues in data.items():
            for tissue, mapQs in tissues.items():
                for mapQ, kmerQs in mapQs.items():
                    for kmerQ, count in kmerQs.items():
                        print(contig,context,tissue,mapQ,kmerQ,count, sep="\t", file=opportunities_out)
    return


def indel_opportunities(filters, contig, bamfile, fastafile, average_coverage_file, tissue_file,
                        window_start=None, window_end=None):
    """
    Count opportunities outside the mask and print them to stdout.
    nonrepetitive_scan() without the calls.
    """
    nonrepetitive_scan(filters, contig, bamfile, fastafile, average_coverage_file, tissue_file,
                       window_start, window_end, opportunities_out=sys.stdout)


def small_indels(filters, contig, bamfile, fastafile, average_coverage_file, tissue_file,
                 window_start=None, window_end=None):
    """
    Call 1-3 bp InDels outside the mask and print them to stdout.
    nonrepetitive_scan() without the opportunities.
    """
    nonrepetitive_scan(filters, contig, bamfile, fastafile, average_coverage_file, tissue_file,
                       window_start, window_end, calls_out=sys.stdout)

def find_lowest_substring(pattern: str) -> str:
    """Find the lexicographically smallest substring or its reverse complement
    in the cyclic permutations of the pattern."""
    size = len(pattern)
    doubled = pattern + pattern
    substrings = [doubled[i:i + size] for i in range(len(pattern))]
    reversed_substrings = [reverse_complement(s) for s in substrings]
    all_substrings = substrings + reversed_substrings
    return min(all_substrings)

def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    complement_map = str.maketrans("ATGC", "TACG")
    complement = sequence.translate(complement_map)
    return complement[::-1]

def haploid_check():
    """
    Placeholder for haploid check function.
    """
    pass


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
                       calls_out=None, opportunities_out=None, fastq_out=None):
    """
    Call InDels and count opportunities inside tandem repeat tracts for several
    repeat models (M0, M1, M2) in one pass over the BAM.

    tr_beds maps a model name to its BED file. calls_out and opportunities_out
    map a model name to an open file for its calls and for its opportunity
    counts; a model left out of one gets no such output. fastq_out, when given,
    gets every read carrying a call in any model, once, as forward-strand FASTQ
    for haploid remapping.

    For each model the output is exactly what that model gives on its own
    (tandem_repeat_indels() and tandem_repeat_opportunities() are this function
    with one model); the pileup, the coverage filter and the work per read are
    shared between the models.

    Calls mirror small_indels() but:
    - Use the TR BED as an *inclusion* filter (positions must be inside a TR interval).
    - InDel size filter: abs(indel) % k == 0 (any size, not just 1-3 bp).
    - Homopolymer filter is NOT applied (k=1 homopolymer InDels are the signal).
    - Extra output columns: canonical_motif, k, copies, trf_purity, model_purity, model.
    - When a position overlaps multiple TR intervals matching the InDel size,
      the interval with the largest k is used (most specific period).

    Opportunities mirror indel_opportunities(): reads passing the quality
    filters at each position inside a TR interval, counted per (canonical_motif,
    k, tissue, mapQ, kmerQ). The kmerQ window uses the period k as the InDel
    length proxy. When a position overlaps multiple TR intervals, the one with
    the smallest k is used.

    With window_start/window_end, only positions in [window_start, window_end)
    are scanned; each position belongs to one window, so windows never report
    the same call twice and their opportunity counts add up to the whole
    contig. Without them the whole contig is scanned.
    """
    calls_out = calls_out or {}
    opportunities_out = opportunities_out or {}

    with open(average_coverage_file, 'r') as f:
        average_coverage_of_species = float(f.read().strip())

    fasta = pysam.FastaFile(fastafile)
    trees = {model: load_tr_bed_with_tree(bedfile, window_start, window_end)
             for model, bedfile in tr_beds.items()}

    d = load_tissue_classification(tissue_file)

    min_mapping_qual = 0
    min_base_qual = 0

    def counts_dict():
        return {str(i): 0 for i in range(94)}

    # counts[model][canonical_motif][k][tissue][mapQ][kmerQ] = count
    counts = {
        model: defaultdict(
            lambda: defaultdict(
                lambda: defaultdict(
                    lambda: defaultdict(counts_dict)
                )
            )
        )
        for model in opportunities_out
    }

    for out in calls_out.values():
        print("#contig", "start", "end", "ref", "alt",
              "indel_type", "unique", "indel_string", "canonical_indel", "tissue", "read_length",
              "np", "ec", "sm", "sx", "sa", "sm_pre", "sx_pre", "ToreadEnd", "kmerQ",
              "position", "readID", "diploid_direction",
              "canonical_motif", "k", "copies", "trf_purity", "model_purity", "model",
              sep="\t", file=out)

    # Base qualities of each read, fetched once rather than at every tract
    # position it covers: copying a whole HiFi read's qualities at each
    # position is most of the cost of counting opportunities. Keyed by read,
    # start and flag, with the read's end so reads the scan has passed can be
    # dropped every 10 kb.
    qualities = {}
    next_eviction = 0

    fastq_reads = set()

    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for pileupcolumn in bam.pileup(contig=contig,
                                       start=window_start, stop=window_end,
                                       region=None, end=None,
                                       fastafile=fasta, truncate=True,
                                       stepper="nofilter",
                                       min_base_quality=min_base_qual,
                                       min_mapping_quality=min_mapping_qual):

            if not (average_coverage_of_species / 2
                    <= pileupcolumn.nsegments
                    <= average_coverage_of_species * 2):
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
                qualities = {key: value for key, value in qualities.items() if value[0] > position}
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
                    key = (query_name, alignment.reference_start, alignment.flag)
                    cached = qualities.get(key)
                    if cached is None:
                        cached = qualities[key] = (alignment.reference_end, alignment.query_qualities)
                    base_qualities = cached[1]
                    mapQ = alignment.mapping_quality
                    tissue_type = d.get(query_name, "NA")

                    for model, tr in opportunity_tracts:
                        k = tr['k']
                        # kmerQ window: 2 bases flanking a k-length window
                        kmerq = (base_qualities[read_pos - 1:read_pos + 1]
                                 + base_qualities[read_pos + k - 1:read_pos + k + 1])
                        if not kmerq:
                            continue
                        counts[model][tr['canonical_motif']][k][tissue_type][str(mapQ)][str(min(kmerq))] += 1

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

                np = ec = sm = sx = sa = sm_pre = sx_pre = "NA"
                read_pos_after = pileupread.query_position_or_next
                readstart      = alignment.reference_start
                readend        = alignment.reference_end
                ToreadEnd      = min(position - readstart, readend - position)

                read_seq       = alignment.query_sequence
                base_qualities = alignment.query_qualities

                if alignment.has_tag("np"):
                    np = alignment.get_tag("np")
                    ec = alignment.get_tag("ec")
                if alignment.has_tag("sa"):
                    sm     = alignment.get_tag("sm")[read_pos + 1]
                    sx     = alignment.get_tag("sx")[read_pos + 1]
                    sa     = subread_array(alignment.get_tag("sa"), read_pos)
                    sm_pre = alignment.get_tag("sm")[read_pos]
                    sx_pre = alignment.get_tag("sx")[read_pos]

                diploid_strand = "-" if alignment.is_reverse else "+"
                tissue_type = d.get(query_name, "NA")

                if indel < 0:
                    # DELETION
                    hcheck = fasta.fetch(contig, position,
                                         position + indel_length + 2)
                    deletion_context = hcheck[1:-1]
                    del_q = (base_qualities[max(0, read_pos - 1):read_pos + 1]
                             + base_qualities[read_pos_after + 1:read_pos_after + 3])
                    kmerQ = min(del_q) if del_q else 0

                    row = (contig,
                           read_pos + 1, read_pos + 1,
                           hcheck[0:-1], hcheck[0],
                           "deletion", seen, deletion_context, canonical_repeat_indel(deletion_context),
                           tissue_type, len(read_seq),
                           np, ec, sm, sx, sa, sm_pre, sx_pre, ToreadEnd, kmerQ,
                           position + 1, query_name, diploid_strand)
                else:
                    # INSERTION
                    insertion_context = read_seq[read_pos + 1:read_pos + indel_length + 1]
                    ins_q = (base_qualities[max(0, read_pos):read_pos + 2]
                             + base_qualities[max(0, read_pos) + indel_length:
                                              read_pos + indel_length + 2])
                    kmerQ = min(ins_q) if ins_q else 0

                    row = (contig,
                           read_pos + 1, read_pos + indel_length + 1,
                           read_seq[read_pos],
                           read_seq[read_pos:read_pos + indel_length + 1],
                           "insertion", seen, insertion_context, canonical_repeat_indel(insertion_context),
                           tissue_type, len(read_seq),
                           np, ec, sm, sx, sa, sm_pre, sx_pre, ToreadEnd, kmerQ,
                           position + 1, query_name, diploid_strand)

                for model, tr in call_tracts:
                    print(*row,
                          tr['canonical_motif'], tr['k'], tr['copies'],
                          tr['trf_purity'], tr['model_purity'], tr['model'],
                          sep="\t", file=calls_out[model])

                # The read, for remapping to the haploid assembly
                if fastq_out is not None and query_name not in fastq_reads:
                    fastq_reads.add(query_name)
                    qual = ''.join(chr(x + 33) for x in alignment.get_forward_qualities())
                    fastq_out.write(f"@{query_name}\n{alignment.get_forward_sequence()}\n+\n{qual}\n")

    for model, out in opportunities_out.items():
        print("contig", "canonical_motif", "k", "tissue", "mapQ", "kmerQ", "count",
              sep="\t", file=out)
        for motif, ks in counts[model].items():
            for period, tissues in ks.items():
                for tissue, mapQs in tissues.items():
                    for mapQ, kmerQs in mapQs.items():
                        for kmerQ, count in kmerQs.items():
                            print(contig, motif, period, tissue, mapQ, kmerQ, count,
                                  sep="\t", file=out)

    fasta.close()
    return


def tandem_repeat_indels(tr_bedfile, contig, bamfile, fastafile,
                         average_coverage_file, tissue_file,
                         window_start=None, window_end=None):
    """
    Call InDels whose size is a multiple of the repeat period, scanning only
    positions *inside* tandem repeat regions defined by an M0, M1 or M2 BED
    file, and print them to stdout. tandem_repeat_scan() for one model.
    """
    tandem_repeat_scan({"model": tr_bedfile}, contig, bamfile, fastafile,
                       average_coverage_file, tissue_file, window_start, window_end,
                       calls_out={"model": sys.stdout})


def tandem_repeat_opportunities(tr_bedfile, contig, bamfile, fastafile,
                                 average_coverage_file, tissue_file,
                                 window_start=None, window_end=None):
    """
    Count opportunities (non-InDel read positions) inside tandem repeat
    regions and print them to stdout. tandem_repeat_scan() for one model.
    """
    tandem_repeat_scan({"model": tr_bedfile}, contig, bamfile, fastafile,
                       average_coverage_file, tissue_file, window_start, window_end,
                       opportunities_out={"model": sys.stdout})


def resolve_haplotypes(contig, window_start, window_end, bamfile, fastafile, average_coverage_file):
    # Only positions in [window_start, window_end) are scanned. Regions are
    # reported with the full span of their reads, so they can reach outside
    # the window; the workflow merges regions from neighbouring windows.
    fasta = pysam.FastaFile(fastafile)
    #average_coverage_of_species = extract_avg_cov(average_coverage_file)
    #print("works?")
    min_mapping_qual = 0
    min_base_qual = 0
    previous_start = 0
    previous_end = 0
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
                    if baseq < 40:
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
                    if len(intervals[0]) > 1: # number of identical bases i need to see to say there is a unresolved haplotype
                        #print(intervals)
                        start = min(intervals[0])
                        end = max(intervals[1])
                        #print(intervals,position+1,start,end)
                        if (previous_start == 0) and (previous_end==0):
                            previous_start = start
                            previous_end = end
                        elif start < previous_end: # extending region
                            previous_end = end
                        elif start > previous_end: # new region
                            print(contig, previous_start, previous_end, sep="\t")
                            previous_start = start
                            previous_end = end
    if previous_start != 0 and previous_end != 0:
        print(contig, previous_start, previous_end, sep="\t")