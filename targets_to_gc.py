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
  args = parser.parse_args()
  return args


def main():
  args = parse_args()
  for line in sys.stdin:
      field = line.strip().split('\t')
      target = field[2]
      gc = (target.count('g') + target.count('G') +
            target.count('c') + target.count('C')) * 5
      sys.stdout.write('\t'.join([target, str(gc)]) + '\n')


##############################################
if __name__ == "__main__":
  sys.exit(main())
