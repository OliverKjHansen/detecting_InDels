"""
The haploid check's CIGAR walk as the review found it, at commit 542fa20, with
findings H1-H4 of ../nonrepetitive_indel_review.md in it. The pipeline's check
(scripts/filtering_haploid_mapping.py, then filtering_haploid_mappping.py) has since been fixed; this copy keeps
the review's numbers reproducible, since they are what the cluster's tables
recorded.

On a flipped read the pipeline handed it `read_length - start - 1`, and
`start` otherwise.
"""
import re


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
