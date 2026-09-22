import sys
import pysam
import csv
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
  
def main():
  # Usage: python filtering_haploid_mappping.py <calls.tsv> [<calls.tsv> ...] <locations.tsv>
  # Several call tables (e.g. one per contig) are checked against the locations,
  # which are read only once, and written out as one table with one header.
  if len(sys.argv) < 3:
    sys.exit("Usage: filtering_haploid_mappping.py <calls.tsv> [<calls.tsv> ...] <locations.tsv>")
  *calls_files, locations_file = sys.argv[1:]

  haploid_alignment = pandas.read_csv(locations_file, sep='\t', names=haploid_header)
  haploid_alignment = haploid_alignment.set_index("readID")

  for file_number, calls_file in enumerate(calls_files):
    check_calls(calls_file, haploid_alignment, print_header=(file_number == 0))


def check_calls(calls_file, haploid_alignment, print_header):
  """Join one call table to the haploid alignments and print the checked calls."""
  indels = pandas.read_csv(calls_file, sep='\t')

  merged = indels.join(haploid_alignment, on="readID")
  merged = merged.drop_duplicates(subset=["readID","start"]) # filter when a read maps to multiple locations in the haploid assembly
  merged = merged.rename(columns=lambda c: c.lstrip("#"))

  new_columns = ["haploid_check","distance_to_indel_in_hap","softclip_haploid", "other_ops"]
  header_cols = [col for col in merged.columns if col != "cigar_hap"] + new_columns
  if print_header:
    print("\t".join(header_cols))

  for row in merged.itertuples(index=False):
    read_id = row.readID
    read_position = row.start # prettysure i make it 1 base when i print it: check #1-based
    cigar = row.cigar_hap
    indel_type = row.indel_type
    indel_string = row.indel_string
    diploid_direction = row.diploid_direction #change
    bitflag_hap = row.bitflag_hap
    read_length = row.read_length
    # quality = int(row.kmerQ)
    # if quality < 93:
    #   continue 

    if bitflag_hap == 16:
      haploid_direction = "-"
    elif bitflag_hap == 0:
      haploid_direction = "+"
    else:
      #print("low quality mapping, skipping read:", read_id, bitflag_hap)
      continue
    #if pandas.isna(row.cigar_hap): #should not exist
      #print(row)
    #if read_id == "m84108_251007_135540_s2/186516316/ccs":
    if diploid_direction != haploid_direction:
      read_position = read_length - read_position - 1

    res,target,softclip,thers = find_indel_position_in_read(cigar, read_position, indel_type, indel_string)
  
    out_cols = [col for col in merged.columns if col != "cigar_hap"]
    values = [getattr(row, col) for col in out_cols]

    print(*values,res,target,softclip,thers, sep="\t")
    #print(res,target,softclip,thers,read_id,row.start, read_position, diploid_direction, haploid_direction, sep=" ")
    # if res != "no_indel_found":
    #   print(res,target,softclip,thers,read_id,diploid_direction, haploid_direction,sep=" ")
      
      
    # headers = merged.headers
    # new = ["distance_to_diploid_indel", softclip_haploid]
    #print(row,res,target,softclip,thers, sep="\t")



  #print(merged)




if __name__ == "__main__":
  main()

#   cigar = [(10, 'M'), (2, 'D'), (5, 'M'), (3, 'I'), (8, 'D')]
# # Search backwards from index 3 (the 'D')
#   result = find_previous_op_with_distance(cigar, len(cigar)-2, 'D')
#   print(result)
