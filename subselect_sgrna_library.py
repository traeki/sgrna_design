#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import bisect
import collections
import copy
import itertools
import logging
import os.path
import pdb
import pprint
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


def two_antisense_deglib(gene, target_group):
  chosen_subset = list()
  target_group = list(target_group)
  specific_targets = [x for x in target_group if x.specificity > 30]
  anti_sense_targets = [x for x in specific_targets if not x.sense_strand]
  anti_sense_targets.sort(key=lambda x:x.offset)
  try:
    first = itertools.dropwhile(lambda x: False, anti_sense_targets).next()
    chosen_subset.append(first)
    second = itertools.dropwhile(lambda x: x.offset < (first.offset + 20),
                                 anti_sense_targets).next()
    chosen_subset.append(second)
  except StopIteration:
    n = len(chosen_subset)
    logging.warn(
        '{n} non-overlapping antisense guides for {gene}'.format(**vars()))
    return chosen_subset
  return chosen_subset
SUBSELECTOR_REGISTRY['two_antisense_deglib'] = two_antisense_deglib


def three_template_deglib(gene, target_group):
  chosen_subset = list()
  target_group = list(target_group)
  specific_targets = [x for x in target_group if x.specificity > 30]
  sense_targets = [x for x in specific_targets if x.sense_strand]
  sense_targets.sort(key=lambda x:x.offset)
  try:
    first = itertools.dropwhile(lambda x: False, sense_targets).next()
    chosen_subset.append(first)
    second = itertools.dropwhile(lambda x: x.offset < (first.offset + 20),
                                 sense_targets).next()
    chosen_subset.append(second)
    third = itertools.dropwhile(lambda x: x.offset < (second.offset + 20),
                                 sense_targets).next()
    chosen_subset.append(third)
  except StopIteration:
    n = len(chosen_subset)
    logging.warn(
        '{n} non-overlapping template guides for {gene}'.format(**vars()))
    return chosen_subset
  return chosen_subset
SUBSELECTOR_REGISTRY['three_template_deglib'] = three_template_deglib


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
  targets = (sgrna_target.from_tsv(x)
                  for x in args.input_tsv_file
                  if not x.startswith('#'))
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
