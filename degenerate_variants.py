#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

import logging
import math
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
TAIL = range(0,8)
MIDDLE = range(8, 13)
SEED = range(13,20)

def revcomp(x):
  return x.translate(DNA_PAIRINGS)[::-1]

def fractioned_variants(target, count):
  orig_count = int(count * 0.2 * 0.4)
  tail_count = int(count * 0.2 * 0.6)
  middle_count = int(count * 0.2)
  seed_count = int(count * 0.2)
  double_count = count - sum([orig_count, tail_count, middle_count, seed_count])
  variants = set()
  single_iterator = random_variants(target, 1)
  double_iterator = random_variants(target, 2)
  current_goal = double_count
  while len(variants) < current_goal:
    variants.add(double_iterator.next())
  current_goal += tail_count
  while len(variants) < current_goal:
    variant = single_iterator.next()
    (candidate, indices, weakness) = variant
    if indices[0] in TAIL:
      variants.add(variant)
  current_goal += middle_count
  while len(variants) < current_goal:
    variant = single_iterator.next()
    (candidate, indices, weakness) = variant
    if indices[0] in MIDDLE:
      variants.add(variant)
  current_goal += seed_count
  while len(variants) < current_goal:
    variant = single_iterator.next()
    (candidate, indices, weakness) = variant
    if indices[0] in SEED:
      variants.add(variant)
  variants = list(variants)
  while len(variants) < count:
    variants.append((target, [], 0))
  return variants

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

def random_variants(target, number):
  if len(target) != len(COST_VECTOR):
    logging.fatal('target and COST_VECTOR length not equal.'.format(**vars()))
    sys.exit(2)
  target = target.upper()
  while True:
    keep = True
    indices = tuple(random.sample(xrange(len(target)), number))
    new_target = list(target)
    for index in indices:
      new_target[index] = random.choice(BASES)
      if new_target[index] == target[index]:
        keep = False
    if keep:
      yield (''.join(new_target), indices, sum([COST_VECTOR[index] for index in indices]))
