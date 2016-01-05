#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

import logging
import random
import string
import sys

from sgrna_target import sgrna_target

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

class Error(Exception):
  pass

class SampleError(Error):
  pass


DNA_PAIRINGS = string.maketrans('atcgATCG', 'tagcTAGC')

BASES = 'ATCG'

COST_VECTOR=[10,10,10,10,10,10,10,10,19,19,19,19,19,28,28,28,28,28,28,28]

def revcomp(x):
  return x.translate(DNA_PAIRINGS)[::-1]

def all_single_variants(target):
  if len(target) != len(COST_VECTOR):
    logging.fatal('target and COST_VECTOR length not equal.'.format(**vars()))
    sys.exit(2)
  target = target.upper()
  for i in xrange(len(target)):
    for x in BASES:
      if x != target[i]:
        new_target = list(target)
        new_target[i] = x
        yield (''.join(new_target), COST_VECTOR[i])

def n_double_variants(target, n):
  generator = random_double_variants(target)
  for i in xrange(n):
    yield generator.next()

def random_double_variants(target):
  if len(target) != len(COST_VECTOR):
    logging.fatal('target and COST_VECTOR length not equal.'.format(**vars()))
    sys.exit(2)
  target = target.upper()
  while True:
    i, j = random.sample(xrange(len(target)), 2)
    new_target = list(target)
    new_target[i] = random.choice(BASES)
    new_target[j] = random.choice(BASES)
    if i == j or new_target[i] == target[i] or new_target[j] == target[j]:
      continue
    yield (''.join(new_target), COST_VECTOR[i] + COST_VECTOR[j])
