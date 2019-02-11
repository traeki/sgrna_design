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


def eval_specificity(controls, fasta):
  specific_controls = list()
  # Is there a bowtie index yet?
  if not os.path.exists(fasta + '.1.ebwt'):
    command = ['bowtie-build', fasta, fasta]
    build_job = subprocess.Popen(command)
    if build_job.wait() != 0:
      logging.fatal('Failed to build bowtie index')
      sys.exit(build_job.returncode)
  # Generate faked FASTQ file
  phredString = 'I4!=======44444++++++++'  # 33333333222221111111NGG
  threshold = 39
  fastq_tempfile, fastq_name = tempfile.mkstemp()
  with contextlib.closing(os.fdopen(fastq_tempfile, 'w')) as fastq_file:
    for name, t in controls.items():
      fullseq = t.sequence_with_pam()
      fastq_file.write(
          '@{name}\n{fullseq}\n+\n{phredString}\n'.format(**vars()))
  (al_tempfile, al_name) = tempfile.mkstemp()
  (max_tempfile, max_name) = tempfile.mkstemp()
  (un_tempfile, un_name) = tempfile.mkstemp()
  (specific_tempfile, specific_name) = tempfile.mkstemp()
  command = ['bowtie']
  command.extend(['-S'])  # output SAM
  command.extend(['--nomaqround'])  # don't do rounding
  command.extend(['-q'])  # input is fastq
  command.extend(['-a'])  # report each non-specific hit
  command.extend(['--best'])  # judge the *closest* non-specific match
  command.extend(['--tryhard'])  # judge the *closest* non-specific match
  command.extend(['--chunkmbs', 256])  # memory setting for --best flag
  command.extend(['-p', 6])  # how many processors to use
  command.extend(['-n', 3])  # allowable mismatches in seed
  command.extend(['-l', 15])  # size of seed
  command.extend(['-e', threshold])  # dissimilarity sum before not non-specific hit
  command.extend(['-m', 1])  # discard reads with >1 alignment
  command.extend(['--al', al_name])
  command.extend(['--un', un_name])
  command.extend(['--max', max_name])
  command.append(fasta)  # index base, built above
  command.append(fastq_name)  # faked fastq temp file
  command.append(specific_name)  # unique hits
  command = [str(x) for x in command]
  logging.info(' '.join(command))
  bowtie_job = subprocess.Popen(command)
  if bowtie_job.wait() != 0:
    sys.exit(bowtie_job.returncode)
  unaligned_reads = SeqIO.parse(un_name, 'fastq-sanger')
  for x in unaligned_reads:
    c = controls[x.name]
    c.specificity = 39
    specific_controls.append(c)
  return specific_controls


def shuffle_targets(targets, needed):
  controls = dict()
  for i in range(needed):
    t = copy.deepcopy(random.choice(targets))
    bases = list(t.target)
    random.shuffle(bases)
    t.target = ''.join(bases)
    t.specificity = 0
    controls[t.id_str()] = t
  return controls


def parse_args():
  """Read in the arguments for the sgrna library construction code."""
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--genomic_background', type=str, required=True,
                      help='FASTA genome to check for specificity eval.')
  parser.add_argument('--input_tsv_file', type=file, help=':')
  parser.add_argument('--output_tsv_file', type=str, help=':', default=None)
  parser.add_argument('--needed', type=int, help=':', default=100)
  args = parser.parse_args()
  if args.output_tsv_file is None:
    base, ext = os.path.splitext(args.input_tsv_file.name)
    args.output_tsv_file =  base + '.controls'
  return args


def main():
  args = parse_args()
  logging.info('Reading in targets from {0}'.format(args.input_tsv_file.name))
  targets = [sgrna_target.from_tsv(x)
                  for x in args.input_tsv_file
                  if not x.startswith('#')]
  controls = shuffle_targets(targets, args.needed)
  specific_controls = eval_specificity(controls, args.genomic_background)
  logging.info('Writing controls to {0}'.format(args.output_tsv_file))
  with open(args.output_tsv_file, 'w') as output_file:
    for x in specific_controls:
      output_file.write(str(x) + '\n')


##############################################
if __name__ == "__main__":
  sys.exit(main())
