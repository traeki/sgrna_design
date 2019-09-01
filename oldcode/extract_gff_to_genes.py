#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import logging
import os.path
import sys

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
  parser.add_argument('gff', type=str, help=':')
  args = parser.parse_args()
  return args


def parse_gff_genes(gff_file):
  """Extract gff file genes into a more usable form.

  Args:
    gff_file: Name of gff input file.
  Returns:
    genes: list of (gene, chrom, start, end, strand) entries.
  """
  logging.info('Parsing gff file.')
  genes = list()
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
      gene = None
      for kv in attributes.split(';'):
        k, v = kv.split('=')
        if k.lower() == 'name':
          gene = v
          break
      if gene: # Skip over genes that didn't have a gene field for now
        genes.append((gene, chrom, int(start), int(end), strand))
  return genes


def main():
  args = parse_args()
  genes = parse_gff_genes(args.gff)
  for gene in genes:
    sys.stdout.write('\t'.join([str(x) for x in gene]))
    sys.stdout.write('\n')

##############################################
if __name__ == "__main__":
  sys.exit(main())
