#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import bisect
import collections
import contextlib
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
import tempfile

from Bio import SeqIO

from .sgrna_target import sgrna_target


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

class Error(Exception):
  pass

class SampleError(Error):
  pass


DNA_PAIRINGS = string.maketrans('atcgATCG', 'tagcTAGC')

def revcomp(x):
  return x.translate(DNA_PAIRINGS)[::-1]


def parse_args():
  """Read in the arguments for the sgrna library construction code."""
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--input_tsv_file', type=file, help=':')
  parser.add_argument('--output_tsv_file', type=str, help=':', default=None)
  args = parser.parse_args()
  if args.output_tsv_file is None:
    base, ext = os.path.splitext(args.input_tsv_file.name)
    args.output_tsv_file =  base + '.uncut' + ext
  return args


def contains_site(target):
  bsai = 'GGTCTC'
  spacer = 'T'
  recaps = revcomp(spacer)
  iasb = revcomp(bsai)
  front_overhangs = ['ATGT', 'AAGC', 'GGGA', 'TAGT']
  back_overhangs = ['GTTT']
  for f in front_overhangs:
    for b in back_overhangs:
      oligo = spacer + f + target.target + b + recaps
      if bsai in oligo or iasb in oligo:
        return True
  return False


def main():
  args = parse_args()
  logging.info('Reading in targets from {0}'.format(args.input_tsv_file.name))
  targets = (sgrna_target.from_tsv(x)
                  for x in args.input_tsv_file
                  if not x.startswith('#'))
  filtered = itertools.filterfalse(contains_site, targets)
  logging.info('Writing uncut targets to {0}'.format(args.output_tsv_file))
  with open(args.output_tsv_file, 'w') as output_file:
    for x in filtered:
      output_file.write(str(x) + '\n')


##############################################
if __name__ == "__main__":
  sys.exit(main())
