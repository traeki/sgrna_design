#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import collections
import logging
import os.path
import sys

from Bio import SeqIO

from sgrna_target import sgrna_target


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

class Error(Exception):
  pass

class SampleError(Error):
  pass


def parse_args():
  """Read in the arguments for the sgrna library construction code."""
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--input_tsv_file', type=file, help=':', required=True)
  parser.add_argument(
          '--comparison_tsv_file', type=file, help=':', required=True)
  parser.add_argument(
      '--output_tsv_file_name', type=str, help=':', default=None)
  args = parser.parse_args()
  if args.output_tsv_file_name is None:
    base, ext = os.path.splitext(args.input_tsv_file.name)
    args.output_tsv_file_name =  base + '.matched' + ext
  return args


def main():
  args = parse_args()
  logging.info('Reading targets from {0}'.format(args.input_tsv_file.name))
  targets = (sgrna_target.from_tsv(x) for x in args.input_tsv_file)
  logging.info('Reading targets from {0}'.format(args.comparison_tsv_file.name))
  comps = (sgrna_target.from_tsv(x) for x in args.comparison_tsv_file)
  comp_seqs = collections.defaultdict(int)
  for x in comps:
      comp_seqs[x.target] += 1
  chosen = [x for x in targets if
               1 == comp_seqs[x.target] and
               39 == x.specificity]
  logging.info('Writing targets to {0}'.format(args.output_tsv_file_name))
  with open(args.output_tsv_file_name, 'w') as output_file:
    for x in chosen:
      output_file.write(str(x) + '\n')


##############################################
if __name__ == "__main__":
  sys.exit(main())
