import pysam
import sys
import csv
from collections import defaultdict
from intervaltree import IntervalTree, Interval

def load_filters_with_tree(bedfile): #for intervals interval trees are quicker then the previous sweep-line approach Ol(og(n)) = O(n)
    """
    Reads bed file and creates an IntervalTree for O(log n) lookups.
    Much faster than sweep-line when checking if a position is in any interval.
    
    Returns:
        IntervalTree: Tree for fast interval overlap queries
    """
    tree = IntervalTree()
    
    with open(bedfile, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            start, end = int(columns[1]), int(columns[2])
            # Create interval [start, end) - end is exclusive
            tree[start:end] = True
    
    return tree



def small_indels_dict(microsatfile, contig, bamfile, fastafile, average_coverage_file,tissue_file):

    with open(average_coverage_file, 'r') as f: # load in the average coverage smaple
        average_coverage_of_species = float(f.read().strip())

    fasta = pysam.FastaFile(fastafile)
    intervals, metadata = load_intervals(microsatfile) # artifact from microsat script, but might be useful later
    
    #####
    with open(tissue_file, "r") as f:
        tissue_classification = csv.reader(f, delimiter="\t")
        d = dict(tissue_classification)
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
                                       start=None,
                                       stop=None,
                                       region=None,
                                       end=None,
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

def large_indels(contig, bamfile, fastafile, average_coverage_file,tissue_file, microsatfile):

    with open(average_coverage_file, 'r') as f: # load in the average coverage smaple
        average_coverage_of_species = float(f.read().strip())

    fasta = pysam.FastaFile(fastafile)
    intervals, metadata = load_intervals(microsatfile) # artifact from microsat script, but might be useful later
    
    #####
    with open(tissue_file, "r") as f:
        tissue_classification = csv.reader(f, delimiter="\t")
        d = dict(tissue_classification)
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
                                       start=None,
                                       stop=None,
                                       region=None,
                                       end=None,
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

def indel_opportunities(filters, contig, bamfile, fastafile, average_coverage_file,tissue_file):

    with open(average_coverage_file, 'r') as f: # load in the average coverage smaple
        average_coverage_of_species = float(f.read().strip())

    fasta = pysam.FastaFile(fastafile)
    #intervals, metadata = load_intervals(microsatfile) # artifact from microsat script, but might be useful later
    interval_tree = load_filters_with_tree(filters) # IntervalTree
    #####
    with open(tissue_file, "r") as f:
        tissue_classification = csv.reader(f, delimiter="\t")
        d = dict(tissue_classification)
    #print(d)

    #####
    min_mapping_qual = 0 # to have more contorl over pysam
    min_base_qual = 0 #--||---

    def counts_dict():
        return {str(i): 0 for i in range(94)}

    data = defaultdict(
    lambda: defaultdict(
        lambda: defaultdict(counts_dict)
        )
    )

    #next_iv = 0                              # CHANGED: index of next interval to activate
    #active = []                              # CHANGED: list of (start,end) intervals currently covering position

    with pysam.AlignmentFile(bamfile, "rb") as bam: # print the header
        for pileupcolumn in bam.pileup(contig=contig, #pileup 
                                       start=None,
                                       stop=None,
                                       region=None,
                                       end=None,
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

            ##skipping position filters
            # activate any intervals that start at/before this position
            # while next_iv < len(intervals) and intervals[next_iv][0] - 2 <= position:
            #     active.append(intervals[next_iv])
            #     next_iv += 1

            # # drop intervals that ended before/at this position (end is exclusive)
            # if active:
            #     active = [(s, e) for (s, e) in active if e > position]
            #     #print(f"DEBUG position {position+1}, active intervals: {active}")
            # if active:
            #     continue
            # # do nothing if we are in  microsat intervals
            # # if active[0][0] <= position <= active[0][1]:
            # #     continue
            
            # fitler on positions on reads
            for pileupread in pileupcolumn.pileups: # looking for reads without indels and have the tags we are interested i
                    if pileupread.query_position is None: # we are in an indel and the position is none
                        continue

                    #non indel or point mutation positions

                    #np = ec = sm = sx = sa = sm_pre = sx_pre = "NA" # initialise tags
                    read_pos = pileupread.query_position
                    readstart = pileupread.alignment.reference_start
                    readend = pileupread.alignment.reference_end
                    #ToreadEnd = min([(position-readstart),(readend-position)]) #tjek

                    # if ToreadEnd < 1000: # if we are too close to the end of the read
                    #     continue

                    read_seq = pileupread.alignment.query_sequence
                    base_qualities = pileupread.alignment.query_qualities
                    mapQ = pileupread.alignment.mapping_quality

                    #readoverlap = min(start - readstart, readend - end)

                    # if pileupread.alignment.has_tag("np"):
                    #     np = pileupread.alignment.get_tag("np") # full read
                    #     ec = pileupread.alignment.get_tag("ec")
                    # if pileupread.alignment.has_tag("sa"): # er et tjek nok?
                    #     sm = pileupread.alignment.get_tag("sm")[read_pos+1] # first base in the insertion
                    #     sx = pileupread.alignment.get_tag("sx")[read_pos+1]
                    #     sa = pileupread.alignment.get_tag("sa")#[pileupread.query_position_or_next] # this is an array so I index later in the function subread_array()
                    #     sa = subread_array(sa,read_pos)
                    #     sm_pre = pileupread.alignment.get_tag("sm")[read_pos] #base just before the insertion                                
                    #     sx_pre = pileupread.alignment.get_tag("sx")[read_pos]

                    
                    #check if we are close to the end of read
                    if (read_pos < 1) or (read_pos > (len(read_seq))-2): # we need at least one base on each side of the kmer
                        continue
                    #make kmer quality check 
                    for indel_length in [1,2,3]: #indel lengths here means kmer length
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

                        tissue_type = d.get(pileupread.alignment.query_name, "NA")

                        lexi = lexicographically(run[1:-1]) #context of kmer without the flanking bases

                        homopolymer = homopolymer_check(run)
                        if homopolymer == True: # filter out homopolymers
                            continue

                        data[lexi][tissue_type][str(mapQ)][str(min(kmerq))] += 1

                        # if min(kmerq) < 93: 
                        #     continue
                    
                        # data[lexi][tissue_type]["filter"] += 1

    print("contig","context","tissue","mapQ","kmerQ","count", sep= "\t", file=sys.stdout) #header for output
    for context, tissues in data.items():
        for tissue, mapQs in tissues.items():
            for mapQ, kmerQs in mapQs.items():
                for kmerQ, count in kmerQs.items():
                    print(contig,context,tissue,mapQ,kmerQ,count, sep="\t", file=sys.stdout)
                    
    return
 
def small_indels(filters, contig, bamfile, fastafile, average_coverage_file,tissue_file):

    with open(average_coverage_file, 'r') as f: # load in the average coverage smaple
        average_coverage_of_species = float(f.read().strip())

    fasta = pysam.FastaFile(fastafile)
    # intervals = load_filters(filters)  #sweep-line approach
    # next_iv = 0                           
    # active = []    

    interval_tree = load_filters_with_tree(filters) # IntervalTree
    
    #####
    with open(tissue_file, "r") as f:
        tissue_classification = csv.reader(f, delimiter="\t")
        d = dict(tissue_classification)
    #print(d)

    #####
    min_mapping_qual = 0 # to have more contorl over pysam
    min_base_qual = 0 #--||---

    # def counts_dict():
    #     return {"filter": 0, "nofilter": 0, "3mer_del":0}
    print("#contig", "start", "end", "ref", "alt",
    "indel_type","unique","indel_string","tissue","read_length",
    "np","ec","sm","sx","sa","sm_pre","sx_pre","ToreadEnd","kmerQ",
    "position","readID", "diploid_direction",
    sep= "\t", file=sys.stdout) #header for output

    # print("contig","indel_type","unique","indel_string","tissue","read_position","read_length",
    # "np","ec","sm","sx","sa","sm_pre","sx_pre", 
    # "position","readID", "diploid_direction",
    # sep= "\t", file=sys.stdout) #header for output
    #
    # multiple reads small indels
    #
    with pysam.AlignmentFile(bamfile, "rb") as bam: # print the header
        for pileupcolumn in bam.pileup(contig=contig, #pileup 
                                       start=None,
                                       stop=None,
                                       region=None,
                                       end=None,
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

            ##skipping position filters
            # activate any intervals that start at/before this position
            # while next_iv < len(intervals) and intervals[next_iv][0] - 2 <= position:
            #     active.append(intervals[next_iv])
            #     next_iv += 1

            # # drop intervals that ended before/at this position (end is exclusive)
            # if active:
            #     active = [(s, e) for (s, e) in active if e > position]
            #     #print(f"DEBUG position {position+1}, active intervals: {active}")
            # if active:
            #     continue
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
                        deletion_context = hcheck[1:-1] #the deletion context
                        #deletion_context = fasta.fetch(contig, position+1, position + 1 + indel_length)
                        lexi = lexicographically(deletion_context)

                        #homopolymer_check

                        homopolymer = homopolymer_check(hcheck)
                        if homopolymer == True: # filter out homopolymers
                            continue

                        tissue_type = d.get(pileupread.alignment.query_name, "NA")
                        #data["deletion"][seen][lexi][tissue_type]["nofilter"] += 1

                        #quality check
                        #AA|ACG|TT
                        #00|xxx|00
                        #print(read_pos,read_pos_after+1)
                        deletion_qualities_4mer = base_qualities[max(0, read_pos-1):read_pos+1]+base_qualities[read_pos_after+1:read_pos_after+3] # need to dobblet check this 
                        #print(deletion_qualities_4mer)
                        kmerQ = min(deletion_qualities_4mer)

                        # if (len(deletion_context) == 1) and (min(deletion_qualities_4mer[0:3])) == 93:
                        #     data["deletion"][seen][lexi][tissue_type]["3mer_del"] += 1

                        # if kmerQ < 93:
                        #     continue
                        
                        if pileupread.alignment.is_reverse == True:
                            diploid_strand = "-"
                        else:
                            diploid_strand = "+"

                        #data["deletion"][seen][lexi][tissue_type]["filter"] += 1
                        print(contig, read_pos+1,read_pos+1,hcheck[0:-1],hcheck[0], #start and end are the same for delettion, becaus its is in the read
                        "deletion", seen, lexi, tissue_type, len(read_seq),
                        np,ec,sm,sx,sa,sm_pre,sx_pre,ToreadEnd,kmerQ,
                        position+1,pileupread.alignment.query_name, diploid_strand,
                        sep="\t", file=sys.stdout)

                        # print("contig", "start", "end", "ref", "alt",
                        # "indel_type","unique","indel_string","tissue","read_length",
                        # "np","ec","sm","sx","sa","sm_pre","sx_pre",
                        # "position","readID", "diploid_direction",
                        # sep= "\t", file=sys.stdout)

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
                        #data["insertion"][seen][lexi][tissue_type]["nofilter"] += 1

                        #quality check
                        #AA|ACG|TT
                        #x0|0x0|0x
                        insertion_qualities_4mer = base_qualities[max(0, read_pos):read_pos+2]+base_qualities[max(0, read_pos)+indel_length:read_pos+indel_length+2] # will count the middle base twice in 1 bp, but it doesnt matter
                        kmerQ = min(insertion_qualities_4mer)

                        # if kmerQ < 93:
                        #     continue

                        if pileupread.alignment.is_reverse == True:
                            diploid_strand = "-"
                        else:
                            diploid_strand = "+"                        
                        #data["insertion"][seen][lexi][tissue_type]["filter"] += 1

                        lexi = lexicographically(insertion_context) 

                        print(contig, read_pos+1,read_pos + indel_length +1,read_seq[read_pos],read_seq[read_pos:read_pos + indel_length +1], #start and end are the same for delettion, becaus its is in the read
                        "insertion", seen, lexi, tissue_type, len(read_seq),
                        np,ec,sm,sx,sa,sm_pre,sx_pre,ToreadEnd,kmerQ,
                        position+1,pileupread.alignment.query_name, diploid_strand,
                        sep="\t", file=sys.stdout)
    fasta.close()

    # print the data
    # print("contig","indel_type","unique","indel_string","tissue","nofilter","3mer","filter", sep= "\t", file=sys.stdout) #header for output
    # for indel_type, uniques in data.items():
    #     for unique, indel_strings in uniques.items():
    #         for indel_string, tissue in indel_strings.items():
    #             for tissue, counts in tissue.items():
    #                 print(contig,indel_type,unique,indel_string,tissue,counts['nofilter'],counts["3mer_del"],counts['filter'], sep="\t", file=sys.stdout)
    return

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

def load_tr_bed_with_tree(bedfile):
    """
    Load an M0 or M1 BED file (11-column, with header) into an IntervalTree.

    Each interval's data payload is a dict with keys:
        canonical_motif, k, copies, trf_purity, model_purity, model

    The first line (header) is skipped automatically.
    Intervals where end <= start are silently skipped.
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
            tree[start:end] = {
                'canonical_motif': cols[3],
                'k':               int(cols[4]),
                'copies':          float(cols[5]),
                'trf_purity':      float(cols[6]),
                'model_purity':    float(cols[7]),
                'model':           cols[8],
            }
    return tree


def tandem_repeat_indels(tr_bedfile, contig, bamfile, fastafile,
                         average_coverage_file, tissue_file):
    """
    Call InDels whose size is a multiple of the repeat period, scanning only
    positions *inside* tandem repeat regions defined by an M0 or M1 BED file.

    Mirrors small_indels() but:
    - Uses TR BED as an *inclusion* filter (positions must be inside a TR interval).
    - InDel size filter: abs(indel) % k == 0 (any size, not just 1-3 bp).
    - Homopolymer filter is NOT applied (k=1 homopolymer InDels are the signal).
    - Extra output columns: canonical_motif, k, copies, trf_purity, model_purity, model.
    - When a position overlaps multiple TR intervals matching the InDel size,
      the interval with the largest k is used (most specific period).
    """
    with open(average_coverage_file, 'r') as f:
        average_coverage_of_species = float(f.read().strip())

    fasta = pysam.FastaFile(fastafile)
    interval_tree = load_tr_bed_with_tree(tr_bedfile)

    with open(tissue_file, "r") as f:
        tissue_classification = csv.reader(f, delimiter="\t")
        d = dict(tissue_classification)

    min_mapping_qual = 0
    min_base_qual = 0

    print("#contig", "start", "end", "ref", "alt",
          "indel_type", "unique", "indel_string", "tissue", "read_length",
          "np", "ec", "sm", "sx", "sa", "sm_pre", "sx_pre", "ToreadEnd", "kmerQ",
          "position", "readID", "diploid_direction",
          "canonical_motif", "k", "copies", "trf_purity", "model_purity", "model",
          sep="\t", file=sys.stdout)

    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for pileupcolumn in bam.pileup(contig=contig,
                                       start=None, stop=None,
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

            # Only process positions INSIDE a TR interval
            hits = interval_tree[position]
            if not hits:
                continue

            indel_count = sum(
                1 for pr in pileupcolumn.pileups
                if (pr.indel != 0) or pr.is_del
            )
            seen = "unique" if indel_count <= 1 else "non-unique"

            for pileupread in pileupcolumn.pileups:
                if pileupread.query_position is None:
                    continue

                indel = pileupread.indel
                if indel == 0:
                    continue

                indel_length = abs(indel)

                # Pick the matching interval with the largest k
                # (most specific period for this InDel size)
                matching = [iv for iv in hits
                            if indel_length % iv.data['k'] == 0]
                if not matching:
                    continue
                best_iv = max(matching, key=lambda iv: iv.data['k'])
                tr = best_iv.data

                np = ec = sm = sx = sa = sm_pre = sx_pre = "NA"
                read_pos      = pileupread.query_position
                read_pos_after = pileupread.query_position_or_next
                readstart     = pileupread.alignment.reference_start
                readend       = pileupread.alignment.reference_end
                ToreadEnd     = min(position - readstart, readend - position)

                read_seq       = pileupread.alignment.query_sequence
                base_qualities = pileupread.alignment.query_qualities

                if pileupread.alignment.has_tag("np"):
                    np = pileupread.alignment.get_tag("np")
                    ec = pileupread.alignment.get_tag("ec")
                if pileupread.alignment.has_tag("sa"):
                    sm     = pileupread.alignment.get_tag("sm")[read_pos + 1]
                    sx     = pileupread.alignment.get_tag("sx")[read_pos + 1]
                    sa     = subread_array(pileupread.alignment.get_tag("sa"), read_pos)
                    sm_pre = pileupread.alignment.get_tag("sm")[read_pos]
                    sx_pre = pileupread.alignment.get_tag("sx")[read_pos]

                upstream, downstream = closest_indel(
                    pileupread.alignment.cigartuples, read_pos, indel_length
                )
                diploid_strand = "-" if pileupread.alignment.is_reverse else "+"
                tissue_type = d.get(pileupread.alignment.query_name, "NA")

                if indel < 0:
                    # DELETION
                    hcheck = fasta.fetch(contig, position,
                                         position + indel_length + 2)
                    deletion_context = hcheck[1:-1]
                    lexi = lexicographically(deletion_context)
                    del_q = (base_qualities[max(0, read_pos - 1):read_pos + 1]
                             + base_qualities[read_pos_after + 1:read_pos_after + 3])
                    kmerQ = min(del_q) if del_q else 0

                    print(contig,
                          read_pos + 1, read_pos + 1,
                          hcheck[0:-1], hcheck[0],
                          "deletion", seen, lexi, tissue_type, len(read_seq),
                          np, ec, sm, sx, sa, sm_pre, sx_pre, ToreadEnd, kmerQ,
                          position + 1, pileupread.alignment.query_name, diploid_strand,
                          tr['canonical_motif'], tr['k'], tr['copies'],
                          tr['trf_purity'], tr['model_purity'], tr['model'],
                          sep="\t", file=sys.stdout)

                if indel > 0:
                    # INSERTION
                    insertion_context = read_seq[read_pos + 1:read_pos + indel_length + 1]
                    lexi = lexicographically(insertion_context)
                    ins_q = (base_qualities[max(0, read_pos):read_pos + 2]
                             + base_qualities[max(0, read_pos) + indel_length:
                                              read_pos + indel_length + 2])
                    kmerQ = min(ins_q) if ins_q else 0

                    print(contig,
                          read_pos + 1, read_pos + indel_length + 1,
                          read_seq[read_pos],
                          read_seq[read_pos:read_pos + indel_length + 1],
                          "insertion", seen, lexi, tissue_type, len(read_seq),
                          np, ec, sm, sx, sa, sm_pre, sx_pre, ToreadEnd, kmerQ,
                          position + 1, pileupread.alignment.query_name, diploid_strand,
                          tr['canonical_motif'], tr['k'], tr['copies'],
                          tr['trf_purity'], tr['model_purity'], tr['model'],
                          sep="\t", file=sys.stdout)

    fasta.close()
    return


def tandem_repeat_opportunities(tr_bedfile, contig, bamfile, fastafile,
                                 average_coverage_file, tissue_file):
    """
    Count opportunities (non-InDel read positions) inside tandem repeat regions.

    For each position inside a TR interval, counts reads that pass quality filters,
    stratified by (canonical_motif, k, tissue, mapQ, kmerQ).
    The kmerQ window uses the period k as the InDel length proxy.
    When a position overlaps multiple TR intervals, the one with the smallest k is used.

    Mirrors indel_opportunities() but uses TR BED as an inclusion filter.
    """
    with open(average_coverage_file, 'r') as f:
        average_coverage_of_species = float(f.read().strip())

    fasta = pysam.FastaFile(fastafile)
    interval_tree = load_tr_bed_with_tree(tr_bedfile)

    with open(tissue_file, "r") as f:
        tissue_classification = csv.reader(f, delimiter="\t")
        d = dict(tissue_classification)

    min_mapping_qual = 0
    min_base_qual = 0

    def counts_dict():
        return {str(i): 0 for i in range(94)}

    # data[canonical_motif][k][tissue][mapQ][kmerQ] = count
    data = defaultdict(
        lambda: defaultdict(
            lambda: defaultdict(
                lambda: defaultdict(counts_dict)
            )
        )
    )

    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for pileupcolumn in bam.pileup(contig=contig,
                                       start=None, stop=None,
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

            # Only process positions INSIDE a TR interval
            hits = interval_tree[position]
            if not hits:
                continue

            # Use the interval with the smallest k at this position
            tr = min(hits, key=lambda iv: iv.data['k']).data
            canonical_motif = tr['canonical_motif']
            k = tr['k']

            for pileupread in pileupcolumn.pileups:
                if pileupread.query_position is None:
                    continue

                read_pos       = pileupread.query_position
                read_seq       = pileupread.alignment.query_sequence
                base_qualities = pileupread.alignment.query_qualities
                mapQ           = pileupread.alignment.mapping_quality

                if (read_pos < 1) or (read_pos > (len(read_seq)) - 2):
                    continue

                # kmerQ window: 2 bases flanking a k-length window
                kmerq = (base_qualities[read_pos - 1:read_pos + 1]
                         + base_qualities[read_pos + k - 1:read_pos + k + 1])
                if not kmerq:
                    continue

                tissue_type = d.get(pileupread.alignment.query_name, "NA")
                data[canonical_motif][k][tissue_type][str(mapQ)][str(min(kmerq))] += 1

    print("contig", "canonical_motif", "k", "tissue", "mapQ", "kmerQ", "count",
          sep="\t", file=sys.stdout)
    for motif, ks in data.items():
        for period, tissues in ks.items():
            for tissue, mapQs in tissues.items():
                for mapQ, kmerQs in mapQs.items():
                    for kmerQ, count in kmerQs.items():
                        print(contig, motif, period, tissue, mapQ, kmerQ, count,
                              sep="\t", file=sys.stdout)

    fasta.close()
    return

def resolve_haplotypes(contig, bamfile, fastafile, average_coverage_file):
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
                                       start=None,
                                       stop=None,
                                       region=None,
                                       end=None,
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