#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import collections
import itertools
import logging
import os.path
import re
import string
import sys

import Bio.SeqIO

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

class Error(Exception):
  pass

class SampleError(Error):
  pass


DNA_PAIRINGS = string.maketrans('atcgATCG', 'tagcTAGC')

def revcomp(x):
  return x.translate(DNA_PAIRINGS)[::-1]


def get_comparison_regions(parent_plasmid):
  """Annotate targets according to overlaps with gff entries.
  Args:
    parent_plasmid: parsed genbank SeqIO object for parent plasmid.
  Returns:
    (roi_front, target, roi_back): The three parts of the must_have region.
  """
  logging.info(
      'Extracting regions from parent plasmid.'.format(**vars()))
  items = [x for x in parent_plasmid.features if x.qualifiers.has_key('label')]
  rois = [x for x in items if x.qualifiers['label'][0] == 'region_of_interest']
  if len(rois) > 1:
    logging.error("More than one region_of_interest label in plasmid.")
    sys.exit(1)
  if len(rois) == 0:
    logging.error("No region_of_interest label in plasmid.")
    sys.exit(1)
  roi = rois[0].extract(parent_plasmid)
  # (Yeah, recreating items is redundant, I'm just redoing it for clarity.)
  items = [x for x in parent_plasmid.features if x.qualifiers.has_key('label')]
  targets = [x for x in items if x.qualifiers['label'][0] == 'target']
  if len(targets) > 1:
    logging.error("More than one target label in plasmid.")
    sys.exit(1)
  if len(targets) == 0:
    logging.error("No target label in plasmid.")
    sys.exit(1)
  target = targets[0].extract(parent_plasmid)
  target_relstart = roi.seq.find(target.seq)
  if target_relstart < 0:
    logging.error("'target' label not inside region_of_interest.")
    sys.exit(1)
  roi_front = roi.seq[:target_relstart]
  roi_back = roi.seq[target_relstart + len(target.seq):]
  return (str(roi_front).upper(),
          str(target.seq).upper(),
          str(roi_back).upper())


def has_items_in_series(sequence, roi_front, target, roi_back):
  """Returns True IFF sequence contains roi_front and roi_back in order."""
  roi_front_start = sequence.find(roi_front)
  roi_back_start = sequence.find(roi_back)
  return (roi_front_start > 0 and
          roi_back_start > 0 and
          (roi_front_start + len(roi_front)) < roi_back_start)


def read_target_library(target_library_file):
  target_library = collections.defaultdict(str)
  for line in target_library_file:
    if line.startswith('#'):
      continue  # Allow comments.
    fields = line.strip().split('\t')
    if len(fields) is not 2:
      logging.error('Expected 2 fields in line:\n{line}'.format(**vars()))
      sys.exit(1)
    target_id, target_seq = fields
    target_seq = target_seq.upper()
    if target_seq in target_library:
      # Bail!  We saw the same target multiple times, maybe with different ids.
      oid = target_library[target_seq]
      nid = target_id
      seq = target_seq
      logging.error(
          'Saw {seq} twice, first as {oid} and then as {nid}'.format(**vars()))
      sys.exit(1)
    target_library[target_seq] = target_id
  return target_library


def read_sequence_records(input_sequence_file):
  input_records = collections.defaultdict(str)
  for record in Bio.SeqIO.parse(input_sequence_file, "fasta"):
    if record.id in input_records:
      # Bail!  We saw the same input multiple times, maybe with different seqs.
      oseq = input_records[record.id]
      nseq = record.seq
      id = record.id
      logging.error(
          'Saw {id} twice, first as\n{oseq}\nand then as\n{nseq}'.format(**vars()))
      sys.exit(1)
    seq = str(record.seq).upper()
    input_records[record.id] = seq
  return input_records


def desperately_search_for_targets(target_library, input_seq):
  for seq in (input_seq, revcomp(input_seq)):
    for target in target_library:
      if input_seq.find(target) > 0:
        return target
  return ''


def parse_args():
  """Read in the arguments for the sgrna sequencing result analysis code."""
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--input_target_lib',
                      type=argparse.FileType('r'), help=':')
  parser.add_argument('--input_sequences',
                      type=argparse.FileType('r'), help=':')
  parser.add_argument('--parent_plasmid', type=str, help=':')
  parser.add_argument('--output_file',
                      type=argparse.FileType('w'), help=':')
  args = parser.parse_args()
  # optional space for post-processing args.
  return args


def main():
  args = parse_args()
  # Read in the target librar(y|ies).
  logging.info(
      'Parsing target lib file: {input_target_lib}.'.format(**vars(args)))
  target_library = read_target_library(args.input_target_lib)
  # Read in the sequence records.
  logging.info(
      'Parsing sequencing output: {input_sequences}.'.format(**vars(args)))
  input_records = read_sequence_records(args.input_sequences)
  # Read in the plasmid.
  logging.info(
      'Parsing parent plasmid file: {parent_plasmid}.'.format(**vars(args)))
  parent_vector = Bio.SeqIO.read(args.parent_plasmid, 'genbank')
  # Determine the Region Of Interest and target region.
  (roi_front, parent_target, roi_back) = get_comparison_regions(parent_vector)
  target_library[parent_target] = 'PARENT'
  # Annotate the input.
  for input_id, input_seq in input_records.iteritems():
    # Try to find the region of interest somewhere.
    roi_matched = has_items_in_series(input_seq,
                                      roi_front,
                                      parent_target,
                                      roi_back)
    if not roi_matched:
      # Maybe this was a reverse primer.
      input_seq = revcomp(input_seq)
      roi_matched = has_items_in_series(input_seq,
                                        roi_front,
                                        parent_target,
                                        roi_back)
    if roi_matched:
      target_start = input_seq.find(roi_front) + len(roi_front)
      roi_back_start = input_seq.find(roi_back)
      target_seq = input_seq[target_start:roi_back_start]
    else:
      target_seq = desperately_search_for_targets(target_library, input_seq)
    # Annotate the rest with the target...
    known_target = target_seq in target_library
    if known_target:
      if not roi_matched:
        prefix = 'BAD_'
      else:
        prefix = ''
      target_id = target_library[target_seq]
    else:
      # ...where possible.
      target_id = ''
    # Write record to output
    fields = (input_id,
              target_id,
              target_seq,
              input_seq,
              roi_matched,
              known_target)
    strfields = (str(x) for x in fields)
    args.output_file.write('\t'.join(strfields) + '\n')


##############################################
if __name__ == "__main__":
  sys.exit(main())
