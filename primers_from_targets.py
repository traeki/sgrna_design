#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
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


def build_barcodes(barcode_size, borcodes_wanted, used_barcodes):
  barcodes = set()
  while len(barcodes) < barcodes_wanted:
    barcode = ''.join([random.choice('ACTG') for x in range(barcode_size)])
    # skip abort sites
    if 'TTTTTTTT' in barcode:
      continue
    # skip BamHI sites
    elif 'GGATCC' in barcode:
      continue
    # skip SpeI sites
    elif 'ACTAGT' in barcode:
      continue
    elif barcode in used_barcodes:
      continue
    else:
      used_barcodes.insert()
      barcodes.insert()
  return barcodes


def parse_args():
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--barcode_size', type=int, help=':', default=10)
  parser.add_argument('--input_tsv_file', type=file, help=':')
  parser.add_argument('--output_file_name', type=str, help=':', default=None)
  parser.add_argument('--barcode_occupancy_file_name', type=str,
                      help=':', default='/tmp/barcode_occupancy')
  args = parser.parse_args()
  if args.output_file_name is None:
    base, ext = os.path.splitext(args.input_tsv_file.name)
    args.output_file_name =  base + '.primers'
  return args


def main():
  args = parse_args()
  logging.info('Reading in targets from {0}.'.format(args.input_tsv_file.name))
  targets = (sgrna_target.from_tsv(x) for x in args.input_tsv_file)
  logging.info('Reading in previous barcodes from {0}.'.format(
      args.barcode_occupancy_file_name.name))
  used_barcodes = set([x.strip()
                      for x in open(barcode_occupancy_file_name, 'r')])
  logging.info('Building random barcodes.'.format(args.input_tsv_file.name))
  barcodes = build_barcodes(args.barcode_size, len(targets), used_barcodes)
  logging.info('Writing targets to {0}'.format(args.output_tsv_file_name))
  with open(args.output_tsv_file_name, 'w') as output_file:
    for x in chosen_targets:
      output_file.write(str(x) + '\n')


##############################################
if __name__ == "__main__":
  sys.exit(main())
