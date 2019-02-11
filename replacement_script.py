#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

import os.path

from .sgrna_target import sgrna_target

bsu_essential = '/Users/jsh/gd/qilab/jsh/refseq/bsu_168/bsu_168.targets.essential.sub.tsv'
bsu_essential_template = 'GCATCGAGCTGAAGGGCAACTAGT{target}GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGAGACCATCCGCCACATGCAGCTC'
bsu_non_essential = '/Users/jsh/gd/qilab/jsh/refseq/bsu_168/bsu_168.targets.non-essential.sub.tsv'
bsu_non_essential_template = 'CGGCAAGCTGACCCTGAAACTAGT{target}GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGAGACCGGTGCCCATCCTGGTCGA'
eco_essential = '/Users/jsh/gd/qilab/jsh/refseq/eco_k12_mg1655/eco_k12_mg1655.targets.essential.sub.tsv'
eco_essential_template = 'CGGCAACATCCTGGGGCAACTAGT{target}GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGAGACCACATCGAGGACGGCAGCG'
eco_non_essential = '/Users/jsh/gd/qilab/jsh/refseq/eco_k12_mg1655/eco_k12_mg1655.targets.non-essential.sub.tsv'
eco_non_essential_template = 'CGCCATGCCCGAAGGCTAACTAGT{target}GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGAGACCTACAAGACCCGCGCCGAG'

files = [bsu_essential,
         bsu_non_essential,
         eco_essential,
         eco_non_essential]
templates = [bsu_essential_template,
             bsu_non_essential_template,
             eco_essential_template,
             eco_non_essential_template]

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
