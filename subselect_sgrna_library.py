#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import bisect
import collections
import copy
import itertools
import logging
import os.path
import random
import re
import subprocess
import string
import sys

from Bio import SeqIO

from sgrna_target import sgrna_target


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

class Error(Exception):
  pass

class SampleError(Error):
  pass


DNA_PAIRINGS = string.maketrans('atcgATCG', 'tagcTAGC')

def revcomp(x):
  return x.translate(DNA_PAIRINGS)[::-1]

SUBSELECTOR_REGISTRY = dict()


def one_antisense_weissman(gene, target_group):
  chosen_subset = list()
  target_group = list(target_group)
  anti_sense_targets = [x for x in list(target_group) if not x.sense_strand]
  count = len(anti_sense_targets)
  if count == 0:
    logging.warn('No targets for gene {gene}.'.format(**vars()))
    return chosen_subset
  specific_targets = [x for x in anti_sense_targets if x.min_unique]
  count = len(specific_targets)
  if count == 0:
    logging.warn('No specific targets for gene {gene}.'.format(**vars()))
    return chosen_subset
  chosen_subset.extend(specific_targets[:1])
  return chosen_subset
SUBSELECTOR_REGISTRY['one_antisense_weissman'] = one_antisense_weissman


def select_four_lsq(gene, target_group):
  chosen_subset = list()
  # Choose one sense-strand target.
  target_group = list(target_group)
  sense_targets = [x for x in list(target_group) if x.sense_strand]
  count = len(sense_targets)
  if count == 0:
    logging.warn('No sense-strand targets for gene {gene}.'.format(**vars()))
  chosen_subset.extend(sense_targets[:1])
  # Choose three anti-sense targets, distributed as far as possible.
  anti_sense_targets = [x for x in list(target_group) if not x.sense_strand]
  count = len(anti_sense_targets)
  if count < 3:
    if count > 0:
      logging.warn(
          'Only {count} anti-sense targets for gene {gene}.'.format(**vars()))
    else:
      logging.warn('No anti-sense targets for gene {gene}.'.format(**vars()))
    chosen_subset.extend(anti_sense_targets)
    return chosen_subset
  def ideal_candidate(x):
    return (x.min_unique <= 12)
  def first_backoff(x):
    return (x.min_unique <= 15)
  def second_backoff(x):
    return (x.min_unique is not None)
  working_set = [x for x in anti_sense_targets if ideal_candidate(x)]
  if len(working_set) < 3:
    # Prefer seed-specific matches.
    working_set = [x for x in anti_sense_targets if first_backoff(x)]
  if len(working_set) < 3:
    # Tolerate some weaker matches.
    working_set = [x for x in anti_sense_targets if second_backoff(x)]
  if len(working_set) < 3:
    # Screw it, just take everything.
    working_set = anti_sense_targets
  # Take the ends.
  working_set.sort(key=lambda x: x.offset)
  chosen_subset.append(working_set[0])
  chosen_subset.append(working_set[-1])
  # And the closest thing to the middle.
  front = working_set[0].offset
  back = working_set[-1].offset
  center = (front + back) / 2
  chosen_subset.append(min(working_set[1:-1],
                           key=lambda x: abs(x.offset - center)))
  return chosen_subset
SUBSELECTOR_REGISTRY['select_four_lsq'] = select_four_lsq


def select_one_sense(gene, target_group):
  targets = [x for x in list(target_group) if x.sense_strand]
  count = len(targets)
  if count == 0:
    logging.warn('No targets for gene {gene}.'.format(**vars()))
    return targets
  return targets[:1]


def three_at_random(gene, target_group):
  """Just pick 3 arbitrarily from the input and return them."""
  return random.sample(list(target_group), 3)
SUBSELECTOR_REGISTRY['three_at_random'] = three_at_random


def parse_args():
  """Read in the arguments for the sgrna library construction code."""
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--input_tsv_file', type=file, help=':')
  parser.add_argument(
      '--output_tsv_file_name', type=str, help=':', default=None)
  parser.add_argument('--subselector', type=str, help=':',
                      default='three_at_random')
  parser.add_argument('--gene_list', type=file, help=':', default=None)
  parser.add_argument('--exclude_listed_genes', help=':',
                      action='store_true', default=False)
  args = parser.parse_args()
  if args.output_tsv_file_name is None:
    base, ext = os.path.splitext(args.input_tsv_file.name)
    args.output_tsv_file_name =  base + '.sub' + ext
  if args.gene_list is not None:
    args.gene_list = set([x.strip() for x in args.gene_list])
  return args


def main():
  args = parse_args()
  subselector_func = SUBSELECTOR_REGISTRY[args.subselector]
  logging.info('Reading in targets from {0}'.format(args.input_tsv_file.name))
  targets = (sgrna_target.from_tsv(x) for x in args.input_tsv_file)
  chosen_targets = list()
  counter = 0
  for gene, cluster in itertools.groupby(targets, lambda x: x.gene):
    if args.gene_list is not None:
      if args.exclude_listed_genes:
        if gene in args.gene_list:
          continue
      else:
        if gene not in args.gene_list:
          continue
    if gene is None:
      continue
    if counter % 500 == 0:
      logging.info('Selecting targets for gene {gene}'.format(**vars()))
    counter += 1
    chosen_targets.extend(subselector_func(gene, cluster))
  logging.info('Writing targets to {0}'.format(args.output_tsv_file_name))
  with open(args.output_tsv_file_name, 'w') as output_file:
    for x in chosen_targets:
      output_file.write(str(x) + '\n')


##############################################
if __name__ == "__main__":
  sys.exit(main())
