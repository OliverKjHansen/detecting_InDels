"""
The haploid check: is each call also in its read's alignment to the haploid
assembly?

    python filtering_haploid_mappping.py <calls.tsv> [<calls.tsv> ...] <locations.tsv> <reads.fastq>

locations.tsv holds the reads' haploid alignment records (read, flag, contig,
position, mapQ, CIGAR) and reads.fastq the reads that were remapped, as
sequenced. Prints the calls with the read's haploid record and the columns
haploid_check, distance_to_indel_in_hap, softclip_haploid and other_ops.
"""
import re
import sys
from collections import defaultdict

import numpy
import pandas
import pysam

from utilities import reverse_complement

haploid_header = [
    "readID",
    "bitflag_hap",
    "contig_hap",
    "pos_1based_hap",
    "map_mapq",
    "cigar_hap"
]

def cigar_events(cigar):
  """
  Walk a haploid CIGAR along its read: (read length, InDels, aligned part,
  softclip_haploid, other_ops). Each InDel is (op, length, start, end), an
  interval of the read in the alignment's orientation: an insertion covers
  its inserted bases, a deletion is the empty interval at its gap. The aligned
  part is the interval between the clips. Clipped bases, hard ones too, are
  read bases, so positions are those of the whole read.
  """
  position, events = 0, []
  first = last = None
  softclip, other_ops = False, None
  for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
    length = int(length)
    if op in "SH":
      softclip = softclip or op == "S"
      other_ops = "HPN" if op == "H" else other_ops
      position += length
      continue
    first = position if first is None else first
    if op in "M=X":
      position += length
    elif op == "I":
      events.append(("I", length, position, position + length))
      position += length
    elif op == "D":
      events.append(("D", length, position, position))
    else:
      # P or N
      other_ops = "HPN"
    last = position
  aligned = (position, position) if first is None else (first, last)
  return position, events, aligned, softclip, other_ops


def mirrored(interval, length):
  """An interval of a read of this length, seen from the other strand."""
  start, end = interval
  return length - end, length - start


def slide(read, op, start, size, deleted):
  """
  The stretch of the read an InDel can slide over, (start, end): the same
  InDel can sit anywhere in it. An insertion of read bases start..start+size
  slides while the base it gives up equals the one it takes on. A deletion of
  `deleted` at the gap before read base `start` slides while the read base it
  passes equals the deleted base at the far end, which then moves round.
  """
  if op == "I":
    first, end = start, start + size
    while first > 0 and read[first - 1] == read[first - 1 + size]:
      first -= 1
    while end < len(read) and read[end] == read[end - size]:
      end += 1
    return first, end
  first, bases = start, deleted
  while first > 0 and read[first - 1] == bases[-1]:
    first, bases = first - 1, bases[-1] + bases[:-1]
  end, bases = start, deleted
  while end < len(read) and read[end] == bases[0]:
    end, bases = end + 1, bases[1:] + bases[0]
  return first, end


def check_call(read, flag, cigar, start, indel_type, indel_string, diploid_direction):
  """
  Whether a call is also in its read's alignment to the haploid assembly, as
  (haploid_check, distance_to_indel_in_hap, softclip_haploid, other_ops).

  read is the read's sequence as sequenced, from the FASTQ that was remapped;
  flag and cigar are the read's haploid alignment record; start, indel_type,
  indel_string and diploid_direction are the call's columns.

  Both alignments are walked in read bases (H1) and their InDels compared as
  intervals of the read as sequenced, so a flipped strand mirrors the InDel,
  not its anchor (H2). Two aligners can place an InDel in a repeat at
  opposite ends of the stretch it can slide over (H3), so the call matches an
  InDel of the same type and size anywhere in that stretch, and only there
  (H4). The distance is from the call to the nearest such InDel, in read
  bases, and there is none without a match.

  Every value is None when the record is not a plain forward or reverse
  alignment (flag 0 or 16), and haploid_check and the distance are None when
  the haploid alignment clips the call off: the check cannot tell either way.
  """
  if flag not in (0, 16):
    return None, None, None, None
  op = "I" if indel_type == "insertion" else "D"
  size = len(indel_string)
  # `start` read bases lead up to the anchor; the InDel follows it. The
  # deleted bases are the diploid assembly's, on the alignment's strand.
  call = (start, start + size) if op == "I" else (start, start)
  deleted = indel_string.upper()
  if diploid_direction == "-":
    call, deleted = mirrored(call, len(read)), reverse_complement(deleted)
  first, end = slide(read, op, call[0], size, deleted)

  length, events, aligned, softclip, other_ops = cigar_events(cigar)
  if flag == 16:
    events = [(event_op, n, *mirrored(interval, length)) for event_op, n, *interval in events]
    aligned = mirrored(aligned, length)
  if not aligned[0] <= call[0] <= call[1] <= aligned[1]:
    return None, None, softclip, other_ops
  distances = [abs(event_start - call[0]) for event_op, n, event_start, event_end in events
               if event_op == op and n == size and first <= event_start and event_end <= end]
  if not distances:
    return False, None, softclip, other_ops
  return True, min(distances), softclip, other_ops


# The columns of a call that check_call() needs.
CALL_COLUMNS = ["readID", "start", "indel_type", "indel_string", "diploid_direction", "read_length"]

# What check_call() gives a call it cannot check.
UNCHECKED = (None, None, None, None)


def haploid_records(locations_file):
  """
  One haploid alignment record per read, indexed by read: its plain forward
  or reverse alignment (flag 0 or 16), the primary, if it has one, and
  otherwise its first record, which check_call() cannot use. Every value is
  kept as text.
  """
  records = pandas.read_csv(locations_file, sep="\t", names=haploid_header, dtype=str, keep_default_na=False)
  usable = records.bitflag_hap.isin(["0", "16"])
  records = records.iloc[numpy.argsort(~usable.to_numpy(), kind="stable")]
  return records.drop_duplicates(subset="readID").set_index("readID")


def check_tables(calls_files, records, reads_file):
  """
  check_call() for every call of every table, one list per table in row order.

  The FASTQ holds every read carrying a call, which can be most of a sample's
  reads, so it is streamed once and each read's calls are checked as it
  passes. A call stays unchecked when its read is not in the FASTQ or has
  another length there than in the call, since its position would then be
  in another read.
  """
  checks, waiting = [], defaultdict(list)
  for number, calls_file in enumerate(calls_files):
    calls = pandas.read_csv(calls_file, sep="\t", usecols=CALL_COLUMNS, dtype=str, keep_default_na=False)
    checks.append([UNCHECKED] * len(calls))
    for row, call in enumerate(calls.itertuples(index=False)):
      waiting[call.readID].append((number, row, call))

  flags = dict(zip(records.index, records.bitflag_hap.astype(int)))
  cigars = dict(zip(records.index, records.cigar_hap))
  unusable = clipped = 0
  with pysam.FastxFile(reads_file) as reads:
    for read in reads:
      for number, row, call in waiting.pop(read.name, []):
        if int(call.read_length) != len(read.sequence):
          continue
        check = check_call(read.sequence, flags.get(read.name), cigars.get(read.name), int(call.start),
                           call.indel_type, call.indel_string, call.diploid_direction)
        unusable += check[2] is None
        clipped += check[2] is not None and check[0] is None
        checks[number][row] = check

  total = sum(len(table) for table in checks)
  unchecked = sum(check[0] is None for table in checks for check in table)
  print(f"{unchecked:,} of {total:,} calls without a haploid check: {unusable:,} on reads without a "
        f"flag 0 or 16 haploid record, {clipped:,} clipped off by the haploid alignment, "
        f"{unchecked - unusable - clipped:,} on reads missing from the FASTQ or of another length there",
        file=sys.stderr)
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

  records = haploid_records(locations_file)
  checks = check_tables(calls_files, records, reads_file)
  for file_number, calls_file in enumerate(calls_files):
    check_calls(calls_file, records, checks[file_number], print_header=(file_number == 0))


def check_calls(calls_file, records, checks, print_header):
  """
  Print one call table with its read's haploid record and its checks. Every
  call is printed: the ones the check could not make have NA in its columns.
  The call's own columns are printed as written.
  """
  calls = pandas.read_csv(calls_file, sep="\t", dtype=str, keep_default_na=False)
  merged = calls.join(records, on="readID").drop(columns="cigar_hap")
  merged = merged.rename(columns=lambda c: c.lstrip("#"))
  new_columns = ["haploid_check", "distance_to_indel_in_hap", "softclip_haploid", "other_ops"]
  for column, values in zip(new_columns, zip(*checks) if checks else [[]] * len(new_columns)):
    merged[column] = pandas.Series(values, index=merged.index, dtype=object)
  merged.to_csv(sys.stdout, sep="\t", index=False, header=print_header, na_rep="NA")


if __name__ == "__main__":
  main()
