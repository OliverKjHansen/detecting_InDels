import sys
import pysam
from collections import defaultdict
import pandas
import re


#input should be contig mutations call and the haploid mapping file
#join on read name
#if 16==sequence have been reverse complemented
#use read_position fra small indel
# LOOP OVER CIGAR BLOCKS


# position_in_read = ifelse(
#       haploid_strand == diploid_strand,
#       haploid_pos_in_read,
#       haploid_read_length - haploid_pos_in_read - 1
#     )

haploid_header = [
    "readID",
    "bitflag_hap",
    "contig_hap",
    "pos_1based_hap",
    "map_mapq",
    "cigar_hap"
]

def find_previous_op_with_distance(cigar, start_index, target_op):
    """
    Search backwards in a CIGAR tuple list for target_op,
    accumulating distance while excluding deletions ('D').
    """
    distance = 0
    for idx in range(start_index, -1, -1):
        length, op = cigar[idx]
        length = int(length)
        
        if op != 'D':
            distance += length
        
        if op == target_op:
            return distance
    
    return None

def find_next_op_with_distance(cigar, start_index, target_op):
    """
    Search forwards in a CIGAR tuple list for target_op,
    accumulating distance while excluding deletions ('D').
    """
    distance = 0
    for idx in range(start_index, len(cigar)):
        length, op = cigar[idx]
        length = int(length)
        
        if op != 'D':
            distance += length
        
        if op == target_op:
            return distance
    
    return None


def find_indel_position_in_read(cigar, position_in_read, indel_type, indel_string):
    """
    Given a CIGAR string and a position in the read (1-based),
    determine if there is an insertion or deletion at that position in a haploid mapping CIGAR string.
    everything should be 1-based
    """
    softclipped_haploid = False
    others = None
    current_pos = 0
    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar)    #print(cigar_tuples)
    target_op = indel_type[0].upper()
    target_size = len(indel_string)
    for idx in range(0,len(cigar_tuples)): # check om previus op er den samme som den efterfølgende
        length,op = cigar_tuples[idx]
        length = int(length)
        if op == 'M' or op == '=' or op == 'X':
          current_pos += length

        elif op == 'I':
            if (current_pos >= position_in_read) and (op == target_op):
                # Found an insertion at or past target position
                forward_distance = abs(current_pos - position_in_read)
                
                # Search backwards for another insertion
                backward_distance = find_previous_op_with_distance(cigar_tuples, idx - 1, op)
                
                # Search forwards for the next insertion # redundant
                #forward_found_distance = find_next_op_with_distance(cigar_tuples, idx + 1, op)
                
                # Find the minimum distance among all occurrences
                distances = [d for d in [forward_distance, backward_distance] if d is not None]
                distance_to_event = min(distances) if distances else forward_distance
                
                # check if the Indel we have found is of the same size as the target indel
                if target_size == length:
                  return True, distance_to_event, softclipped_haploid, others
                else:
                  return False, distance_to_event, softclipped_haploid, others
            #current_pos += length

        elif op == 'D':
            if (current_pos >= position_in_read) and (op == target_op):
                # Found a deletion at or past target position
                forward_distance = abs(current_pos - position_in_read)
                
                # Search backwards for another deletion
                backward_distance = find_previous_op_with_distance(cigar_tuples, idx - 1, op)
                
                # Search forwards for the next deletion
                #forward_found_distance = find_next_op_with_distance(cigar_tuples, idx + 1, op)
                
                # Find the minimum distance among all occurrences
                distances = [d for d in [forward_distance, backward_distance] if d is not None]
                distance_to_event = min(distances) if distances else forward_distance

                # check if the Indel we have found is of the same size as the target indel
                if target_size == length:
                  return True, distance_to_event, softclipped_haploid, others
                else:
                  return False, distance_to_event, softclipped_haploid, others
            # Deletions do not advance the read position
        elif op == "S":
          softclipped_haploid = True
          current_pos += length
         # current_pos =+ length i dont this should advance the position in read for softclip
        elif op == "H" or op == "P" or op == "N":
          others = "HPN"
         # current_pos =+ length i dont this should advance the position in read for hardclip, padding or N
        else:
          others = "unknown"
      
    
    nearest_backward = find_previous_op_with_distance(cigar_tuples, len(cigar_tuples) - 1, target_op)
    nearest_forward = find_next_op_with_distance(cigar_tuples, 0, target_op)
    
    distances = [d for d in [nearest_backward, nearest_forward] if d is not None]
    if distances:
        distance_to_event = min(distances)
        return False, distance_to_event, softclipped_haploid, others
    else:
        return False, None, softclipped_haploid, others
  
def check_call(read, flag, cigar, start, indel_type, indel_string, diploid_direction):
  """
  Whether a call is also in its read's alignment to the haploid assembly, as
  (haploid_check, distance_to_indel_in_hap, softclip_haploid, other_ops).

  read is the read's sequence as sequenced, from the FASTQ that was remapped;
  flag and cigar are the read's haploid alignment record; start, indel_type,
  indel_string and diploid_direction are the call's columns.
  """
  haploid_direction = "-" if flag == 16 else "+"
  position = start
  if diploid_direction != haploid_direction:
    position = len(read) - start - 1
  return find_indel_position_in_read(cigar, position, indel_type, indel_string)


# The columns of a call that check_call() needs.
CALL_COLUMNS = ["readID", "start", "indel_type", "indel_string", "diploid_direction", "read_length"]


def check_tables(calls_files, records, reads_file):
  """
  check_call() for every call of every table, one list per table in row order.

  The FASTQ holds every read carrying a call, which can be most of a sample's
  reads, so it is streamed once and each read's calls are checked as it
  passes. A call is None when its read has no haploid record of flag 0 or 16;
  check_calls() leaves it out.
  """
  checks, waiting = [], defaultdict(list)
  for number, calls_file in enumerate(calls_files):
    calls = pandas.read_csv(calls_file, sep="\t", usecols=CALL_COLUMNS, dtype=str, keep_default_na=False)
    checks.append([None] * len(calls))
    for row, call in enumerate(calls.itertuples(index=False)):
      waiting[call.readID].append((number, row, call))

  flags = dict(zip(records.index, records.bitflag_hap))
  cigars = dict(zip(records.index, records.cigar_hap))
  with pysam.FastxFile(reads_file) as reads:
    for read in reads:
      calls = waiting.pop(read.name, None)
      if calls is None or flags.get(read.name) not in (0, 16):
        continue
      for number, row, call in calls:
        checks[number][row] = check_call(read.sequence, flags[read.name], cigars[read.name], int(call.start),
                                         call.indel_type, call.indel_string, call.diploid_direction)
  return checks


def main():
  # Usage: python filtering_haploid_mappping.py <calls.tsv> [<calls.tsv> ...] <locations.tsv> <reads.fastq>
  # Several call tables (e.g. one per contig) are checked against the locations
  # and the reads, which are read only once, and written out as one table with
  # one header. The reads are the FASTQ that was remapped to the haploid
  # assembly.
  if len(sys.argv) < 4:
    sys.exit("Usage: filtering_haploid_mappping.py <calls.tsv> [<calls.tsv> ...] <locations.tsv> <reads.fastq>")
  *calls_files, locations_file, reads_file = sys.argv[1:]

  haploid_alignment = pandas.read_csv(locations_file, sep='\t', names=haploid_header)
  # One record per read: the first, which is the primary alignment as it sorts first.
  haploid_alignment = haploid_alignment.drop_duplicates(subset="readID").set_index("readID")

  checks = check_tables(calls_files, haploid_alignment, reads_file)
  for file_number, calls_file in enumerate(calls_files):
    check_calls(calls_file, haploid_alignment, checks[file_number], print_header=(file_number == 0))


def check_calls(calls_file, haploid_alignment, checks, print_header):
  """Join one call table to the haploid alignments and print the checked calls."""
  indels = pandas.read_csv(calls_file, sep='\t')

  merged = indels.join(haploid_alignment, on="readID")
  merged = merged.rename(columns=lambda c: c.lstrip("#"))

  new_columns = ["haploid_check","distance_to_indel_in_hap","softclip_haploid", "other_ops"]
  out_cols = [col for col in merged.columns if col != "cigar_hap"]
  if print_header:
    print("\t".join(out_cols + new_columns))

  for row, check in zip(merged.itertuples(index=False), checks):
    if check is None:
      continue
    values = [getattr(row, col) for col in out_cols]
    print(*values, *check, sep="\t")


if __name__ == "__main__":
  main()
