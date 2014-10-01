#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import bisect
import collections
import copy
import itertools
import logging
import os.path
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


def parse_jason_genes(jason_file):
  """Extract jason-format file genes into a more usable form.

  Args:
    jason_file: Name of input file.
  Returns:
    gff_genes: list of (gene, chrom, start, end, strand) entries.
  """
  gff_genes = list()
  for x in open(jason_file):
    if x.startswith('#'):
      continue
    parts = x.strip().split('\t')
    try:
      (name,synonyms,names,start,end,strand,uniqid) = parts
    except ValueError:
      logging.error('Could not parse: {x}'.format(**vars()))
      sys.exit(1)
    gff_genes.append((name, 'NC_000913.2', start, end, strand))
  return gff_genes


def parse_gff_genes(gff_file):
  """Extract gff file genes into a more usable form.

  Args:
    gff_file: Name of gff input file.
  Returns:
    gff_genes: list of (gene, chrom, start, end, strand) entries.
  """
  gff_genes = list()
  for x in open(gff_file):
    if x.startswith('#'):
      continue
    parts = x.strip().split('\t')
    try:
      (chrom,source,feature,start,end,score,strand,frame,attributes) = parts
    except ValueError:
      logging.error('Could not parse: {x}'.format(**vars()))
      sys.exit(1)
    if feature.lower() == 'gene':
      for kv in attributes.split(';'):
        gene = None
        k, v = kv.split('=')
        if k.lower() == 'name':
          gene = v
          break
      if gene: # Skip over genes that didn't have a gene field for now
        gff_genes.append((gene, chrom, int(start), int(end), strand))
  return gff_genes


def parse_args():
  """Read in the arguments for the sgrna library construction code."""
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--input_gff_name', type=str, help=':')
  group.add_argument('--input_jason_name', type=str, help=':')
  parser.add_argument('--gene_list', type=file, help=':', default=None)
  args = parser.parse_args()
  return args


def main():
  args = parse_args()
  gene_list = set([x.strip() for x in args.gene_list])
  if args.input_gff_name:
    genes = parse_gff_genes(args.input_gff_name)
  else:
    genes = parse_jason_genes(args.input_jason_name)
  gff_set = set([x[0].strip() for x in genes])
  difference = gene_list - gff_set
  for x in difference:
    logging.warn('{x} from gff file was not in gene list'.format(**vars()))


##############################################
if __name__ == "__main__":
  sys.exit(main())
