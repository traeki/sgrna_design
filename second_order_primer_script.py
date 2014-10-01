#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

import os.path

from sgrna_target import sgrna_target

tsv_file = '/Users/jsh/gd/proj/P0001/second_order_planning/bsu_168.targets.second_order.tsv'
primer_template = '{target}GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC'

files = [tsv_file]
templates = [primer_template]

for f, t in zip(files, templates):
  base, ext = os.path.splitext(f)
  organism = os.path.basename(f).split('.')[0]
  outfile = open(base + '.oligos', 'w')
  for line in open(f):
    sgt = sgrna_target.from_tsv(line)
    target = sgt.target
    gene = sgt.gene
    offset = sgt.offset
    sdir = 'antisense'
    if sgt.sense_strand:
      sdir = 'sense'
    name = '{organism}.{gene}.{offset}.{sdir}'.format(**vars())
    outfile.write(name + '\t' + t.format(**vars()) + '\n')
