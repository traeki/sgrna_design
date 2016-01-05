#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import bisect
import collections
import copy
import itertools
import logging
import os.path
import random
import re
import shutil
import subprocess
import string
import sys
import tempfile

from Bio import SeqIO
import pysam

from sgrna_target import sgrna_target
import degenerate_variants as degvar


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

class Error(Exception):
  pass

class SampleError(Error):
  pass


DNA_PAIRINGS = string.maketrans('atcgATCG', 'tagcTAGC')

def revcomp(x):
  return x.translate(DNA_PAIRINGS)[::-1]


def extract_targets(infile_name, pam, target_len):
  """Generate the complete list of pam-adjacent potential targets in a genome.

  Args:
    infile_name [str]:  Name of the file containing the source genome.
    pam [str]:          Regexp DNA pattern for the PAM sequence.
    target_len [int]:   How many bases to pull from the adjacent region.
  Returns:
    Iterable sequence of sgrna targets.
  Notes:
    Discards targets containing 'N' bases.
  """
  # TODO(jsh): Do something with "bases" other than N, ATCG.
  logging.info('Extracting target set from {infile_name}.'.format(**vars()))
  fasta_sequences = SeqIO.parse(infile_name, 'fasta')
  raw_targets = dict()
  for seq_record in fasta_sequences:
    genome = seq_record.seq.upper()
    chrom = seq_record.name
    pam = pam.upper()
    reversed_pam = revcomp(pam)
    block = r'(.{' + str(target_len) + r'})'
    pam_pattern = r'(?=(' + block + pam + r'))'
    rev_pattern = r'(?=(' + reversed_pam + block + r'))'
    for hit in re.finditer(pam_pattern, str(genome)):
      if 'N' in hit.group(1):
        continue  # ...Don't target unknown genetic material.
      t = sgrna_target(
                hit.group(2),
                hit.group(1)[-len(pam):],
                chrom,
                hit.start() + 1,
                hit.start() + 1 + target_len,
                False)
      name = t.id_str()
      raw_targets[name] = t
    for hit in re.finditer(rev_pattern, str(genome)):
      if 'N' in hit.group(1):
        continue
      t = sgrna_target(
                revcomp(hit.group(2)),
                revcomp(hit.group(1))[-len(pam):],
                chrom,
                hit.start() + 1 + len(pam),
                hit.start() + 1 + len(pam) + target_len,
                True)
      name = t.id_str()
      raw_targets[name] = t
  logging.info('{0} raw targets.'.format(len(raw_targets)))
  return raw_targets


def parse_target_regions(target_regions_file):
  """Extract target regions into a more usable form.

  Args:
    target_regions_file: Name of input file.
  Returns:
    target_regions: list of (gene, chrom, start, end, strand) entries.
  """
  logging.info('Parsing target region file.')
  target_regions = list()
  for x in open(target_regions_file):
    if x.startswith('#'):
      continue
    parts = x.strip().split('\t')
    try:
      (name,chrom,start,end,strand) = parts
    except ValueError:
      trf = target_regions_file
      logging.error('Could not parse from {trf}: {x}'.format(**vars()))
      sys.exit(1)
    try:
      target_regions.append((name, chrom, int(start), int(end), strand))
    except ValueError:
      x = x.strip()
      logging.warning('Could not fully parse: {x}'.format(**vars()))
      continue
  logging.info(
      'Found {0} target regions in region file.'.format(len(target_regions)))
  return target_regions


def chrom_lengths(fasta_file_name):
  """Get lengths of chromosomes (entries) for fasta file.

  Args:
    fasta_file_name [str]:  Name of the file containing the source genome.
  Returns:
    chrom_lens: dict mapping fasta entry name (chrom) to sequence length.
  """
  chrom_lens = dict()
  fasta_sequences = SeqIO.parse(fasta_file_name, 'fasta')
  for seq_record in fasta_sequences:
    chrom_lens[seq_record.name] = len(seq_record.seq)
  return chrom_lens


def ascribe_specificity(targets, genome_fasta_name, sam_copy):
  """Set up bowtie stuff and repeatedly call mark_specificity_tier."""
  # Is there a bowtie index yet?
  if not os.path.exists(genome_fasta_name + '.1.ebwt'):
    command = ['bowtie-build', genome_fasta_name, genome_fasta_name]
    build_job = subprocess.Popen(command)
    if build_job.wait() != 0:
      logging.fatal('Failed to build bowtie index')
      sys.exit(build_job.returncode)
  # Generate faked FASTQ file
  phredString = '++++++++44444=======!4I'  # 33333333222221111111NGG
  _, fastq_name = tempfile.mkstemp()
  with open(fastq_name, 'w') as fastq_file:
    for name, t in targets.iteritems():
      fullseq = t.sequence_with_pam()
      fastq_file.write(
          '@{name}\n{fullseq}\n+\n{phredString}\n'.format(**vars()))
  for threshold in (95,90,80,70,60,50,40,30,20,11,1):
    mark_unadjusted_specificity_threshold(
        targets, fastq_name, genome_fasta_name, threshold, sam_copy)
  for _, t in targets.iteritems():
    if t.specificity != 0:
      t.specificity -= t.weakness

def mark_unadjusted_specificity_threshold(
        targets, fastq_name, genome_name, threshold, sam_copy):
  """Marks the indicated threshold for any targets that are newly lapsed.

  'Unadjusted' because it does not account for a non-zero relative weakness
  of the guide to the intended target.
  """
  # prep output files
  (_, specific_name) = tempfile.mkstemp()
  # Filter based on specificity
  command = ['bowtie']
  command.extend(['-S'])  # output SAM
  command.extend(['--nomaqround'])  # don't do rounding
  command.extend(['-q'])  # input is fastq
  command.extend(['-a'])  # report each non-specific hit
  command.extend(['--best'])  # judge the *closest* non-specific match
  command.extend(['--chunkmbs', 128])  # memory setting for --best flag
  command.extend(['-p', 7])  # how many processors to use
  command.extend(['-n', 3])  # allowable mismatches in seed
  command.extend(['-l', 15])  # size of seed
  command.extend(['-e', threshold])  # dissimilarity sum before not non-specific hit
  command.extend(['-m', 1])  # discard reads with >1 alignment
  command.append(genome_name)  # index base, built above
  command.append(fastq_name)  # faked fastq temp file
  command.append(specific_name)  # unique hits
  command = [str(x) for x in command]
  logging.info(' '.join(command))
  bowtie_job = subprocess.Popen(command)
  # Check for problems
  if bowtie_job.wait() != 0:
    sys.exit(bowtie_job.returncode)
  if sam_copy:
    shutil.copyfile(specific_name, sam_copy)
  aligned_reads = pysam.Samfile(specific_name)
  for x in aligned_reads:
    # flag 4 means unaligned, so skip those
    if not x.flag & 4:
      t = targets[x.qname]
      if t.specificity < threshold:
        t.specificity = threshold


def label_targets(targets,
                  target_regions,
                  chrom_lens,
                  include_unlabeled,
                  allow_partial_overlap):
  """Annotate targets according to overlaps with gff entries.
  Args:
    targets: the targets to annotate.
    target_regions: the target regions for which to produce annotations
    chrom_lens: mapping from chrom name to sequence length.
    include_unlabeled: if true, add unlabeled versions of full target set.
    allow_partial_overlap: Include targets which only partially overlap region.
  Returns:
    anno_targets: list of targets with added region annotations
  """
  logging.info(
      'Labeling targets based on region file.'.format(**vars()))
  anno_targets = list()
  counter = 0
  # Organize targets by chromosome and then start location.
  per_chrom_sorted_targets = collections.defaultdict(list)
  for name, x in targets.iteritems():
    per_chrom_sorted_targets[x.chrom].append(x)
    if include_unlabeled:
      anno_targets.append(copy.deepcopy(x))
  for x in per_chrom_sorted_targets:
    per_chrom_sorted_targets[x].sort(key=lambda x:(x.start, x.end))
  front, back = 0, 0
  for i, x in enumerate(target_regions):
    (gene, chrom, gene_start, gene_end, gene_strand) = x
    if i % 100 is 0:
      logging.info('Examining gene {i} [{gene}].'.format(**vars()))
    reverse_strand_gene = gene_strand == '-'
    if gene_start >= chrom_lens[chrom]:
      continue
    chrom_targets = per_chrom_sorted_targets[chrom]
    if allow_partial_overlap:
      # Shift back index until target.start >= gene_end
      while (back < len(chrom_targets) and
             chrom_targets[back].start < gene_end):
        back += 1
      # Shift front index until target.end > gene_start
      while (front < len(chrom_targets) and
             chrom_targets[front].end <= gene_start):
        front += 1
    else:
      # Shift back index until target.end > gene_end
      while (back < len(chrom_targets) and
             chrom_targets[back].end <= gene_end):
        back += 1
      # Shift front index until target.start >= gene_start
      while (front < len(chrom_targets) and
             chrom_targets[front].start < gene_start):
        front += 1
    overlap = chrom_targets[front:back]
    if len(overlap) == 0:
      logging.warn('No overlapping targets for gene {gene}.'.format(**vars()))
    for target in overlap:
      if reverse_strand_gene:
        offset = gene_end - target.end
      else:
        offset = target.start - gene_start
      returnable = copy.deepcopy(target)
      returnable.gene = gene
      returnable.offset = offset
      returnable.sense_strand = (reverse_strand_gene == target.reverse)
      anno_targets.append(returnable)
  return anno_targets


def parse_args():
  """Read in the arguments for the sgrna library construction code."""
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--input_fasta_genome_name', type=str, required=True,
                      help='Location of genome file in FASTA format.')
  parser.add_argument('--sam_copy', type=str,
                      help='Copy of sam file from (final) bowtie run.',
                      default=None)
  parser.add_argument('--tsv_file_name', type=str,
                      help='Output file to create.', default=None)
  parser.add_argument('--target_regions_file', type=str, required=True,
                      help='Location of target regions file in tsv format.')
  parser.add_argument('--include_unlabeled', action='store_true',
      default=False,
      help='Output targets even if they overlapped no target region.')
  parser.add_argument('--only_include_fully_overlapping', action='store_false',
      dest='allow_partial_overlap', default=True,
      help='Only label targets which are fully contained in the region.')
  # TODO(jsh) Before using different PAMs, need phred faking flag for PAM.
  parser.add_argument('--pam', default='.gg', type=str,
                      help='NOT YET IMPLEMENTED DO NOT USE!')
  parser.add_argument('--double_variants', default=20, type=int,
                      help='How many random samples of double-degenerate variants to use.')
  # TODO(jsh) Need to take correct fraction of phred-score string to unbreak.
  parser.add_argument('--target_len', default=20, type=int,
                      help='NOT YET IMPLEMENTED DO NOT USE!')
  args = parser.parse_args()
  if args.tsv_file_name is None:
    base = os.path.splitext(args.input_fasta_genome_name)[0]
    args.tsv_file_name =  base + '.targets.all.tsv'
  return args


def main():
  args = parse_args()
  # Build initial list
  clean_targets = extract_targets(args.input_fasta_genome_name,
                                  args.pam,
                                  args.target_len)
  # Spawn degenerate variants
  all_targets=dict()
  counter = 0
  for name, t in clean_targets.iteritems():
    counter += 1
    if random.random() < 0.001:
      logging.info('Varying target {0}: {1}'.format(counter, name))
    all_targets[name] = t
    for variant, weakness in degvar.all_single_variants(t.target):
      new_target = copy.deepcopy(t)
      new_target.target = variant
      new_target.weakness = weakness
      name = new_target.id_str()
      all_targets[name] = new_target
    for variant, weakness in degvar.n_double_variants(t.target, args.double_variants):
      new_target = copy.deepcopy(t)
      new_target.target = variant
      new_target.weakness = weakness
      name = new_target.id_str()
      all_targets[name] = new_target
  logging.info('{0} all targets.'.format(len(all_targets)))
  # Score list
  ascribe_specificity(all_targets, args.input_fasta_genome_name, args.sam_copy)
  # Annotate list
  chrom_lens = chrom_lengths(args.input_fasta_genome_name)
  target_regions = parse_target_regions(args.target_regions_file)
  all_targets = label_targets(all_targets,
                              target_regions,
                              chrom_lens,
                              args.include_unlabeled,
                              args.allow_partial_overlap)
  # Generate output
  total_count = len(all_targets)
  logging.info(
      'Writing {total_count} annotated targets to {args.tsv_file_name}'.format(
          **vars()))
  with open(args.tsv_file_name, 'w') as tsv_file:
    tsv_file.write('#' + sgrna_target.header() + '\n')
    for target in all_targets:
      tsv_file.write(str(target) + '\n')

##############################################
if __name__ == "__main__":
  sys.exit(main())
