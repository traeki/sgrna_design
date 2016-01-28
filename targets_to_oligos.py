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


def parse_args():
  """Read in the arguments for the sgrna library construction code."""
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--input_tsv_file', type=file, help=':')
  parser.add_argument('--output_oligo_file', type=str, help=':', default=None)
  parser.add_argument('--front_overhang', type=str, required=True,
                      help='pVeg:ATGT; pLepA:AAGC; pMelanie:GGGA; pJason:TAGT')
  args = parser.parse_args()
  if args.output_oligo_file is None:
    base, ext = os.path.splitext(args.input_tsv_file.name)
    args.output_oligo_file =  base + '.oligos'
  return args


def main():
  args = parse_args()
  logging.info('Reading in targets from {0}'.format(args.input_tsv_file.name))
  targets = (sgrna_target.from_tsv(x)
                  for x in args.input_tsv_file
                  if not x.startswith('#'))
  logging.info('Writing oligos to {0}'.format(args.output_oligo_file))
  with open(args.output_oligo_file, 'w') as output_file:
    for x in targets:
      bsai = 'GGTCTC'
      spacer = 'T'
      iasb = revcomp(bsai)
      recaps = revcomp(spacer)
      back_overhang = 'GTTT'
      core = args.front_overhang + x.target + back_overhang
      if bsai in core or iasb in core:
        logging.warn('BsaI cut site in target {x}'.format(**vars()))
      output_file.write(bsai + spacer + core + recaps + iasb + '\n')


##############################################
if __name__ == "__main__":
  sys.exit(main())
