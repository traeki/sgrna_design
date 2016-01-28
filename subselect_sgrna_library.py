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


def partition_overlapping(inlist):
  good = list()
  bad = list()
  if not inlist:
    return (good, bad)
  good.append(inlist[0])
  for x in inlist[1:]:
    if x.offset < (good[-1].offset + 20):
      bad.append(x)
    else:
      good.append(x)
  return (good, bad)


def antisense(gene, target_group, wanted):
  target_group = list(target_group)
  target_group.sort(key=lambda x:x.offset)
  # for each specificity threshold, take n_non_overlapping, or just n
  sense_targets = [x for x in target_group if x.sense_strand]
  anti_sense_targets = [x for x in target_group if not x.sense_strand]
  for threshold in [30, 20, 10, 0]:
    specific = lambda x: x.specificity > threshold
    specific_targets = [x for x in anti_sense_targets if specific(x)]
    nonspecific_targets = [x for x in anti_sense_targets if not specific(x)]
    (spaced, overlapped) = partition_overlapping(specific_targets)
    # Try to get non_overlapping, specific guides
    if len(spaced) >= wanted:
      return spaced[:wanted]
    # Okay, can we get enough if we allow overlap?
    if len(specific_targets) >= wanted:
      needed = wanted - len(spaced)
      return spaced + overlapped[:needed]
  # if that never works, dig into nonspecific_targets
  if len(anti_sense_targets) >= wanted:
    needed = wanted - len(specific_targets)
    return specific_targets + nonspecific_targets[:needed]
  # if that doesn't even work, dip into the other strand, and then give up.
  needed = wanted - len(anti_sense_targets)
  return anti_sense_targets + sense_targets[:needed]
SUBSELECTOR_REGISTRY['antisense'] = antisense


def template(gene, target_group, wanted):
  target_group = list(target_group)
  target_group.sort(key=lambda x:x.offset)
  # for each specificity threshold, take n_non_overlapping, or just n
  sense_targets = [x for x in target_group if x.sense_strand]
  anti_sense_targets = [x for x in target_group if not x.sense_strand]
  for threshold in [30, 20, 10, 0]:
    specific = lambda x: x.specificity > threshold
    specific_targets = [x for x in sense_targets if specific(x)]
    nonspecific_targets = [x for x in sense_targets if not specific(x)]
    (spaced, overlapped) = partition_overlapping(specific_targets)
    # Try to get non_overlapping, specific guides
    if len(spaced) >= wanted:
      return spaced[:wanted]
    # Okay, can we get enough if we allow overlap?
    if len(specific_targets) >= wanted:
      needed = wanted - len(spaced)
      return spaced + overlapped[:needed]
  # if that never works, dig into nonspecific_targets and then give up
  needed = wanted - len(specific_targets)
  return specific_targets + nonspecific_targets[:needed]
SUBSELECTOR_REGISTRY['template'] = template


def parse_args():
  """Read in the arguments for the sgrna library construction code."""
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--input_tsv_file', type=file, help=':')
  parser.add_argument(
      '--output_tsv_file_name', type=str, help=':', default=None)
  parser.add_argument('--subselector', type=str, help=':', required=True)
  parser.add_argument('--wanted', type=int, help=':', default=3)
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
    chosen_targets.extend(subselector_func(gene, cluster, args.wanted))
  logging.info('Writing targets to {0}'.format(args.output_tsv_file_name))
  with open(args.output_tsv_file_name, 'w') as output_file:
    for x in chosen_targets:
      output_file.write(str(x) + '\n')


##############################################
if __name__ == "__main__":
  sys.exit(main())
