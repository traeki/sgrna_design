#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

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


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

class Error(Exception):
  pass

class SampleError(Error):
  pass


DNA_PAIRINGS = string.maketrans('atcgATCG', 'tagcTAGC')

def revcomp(x):
  return x.translate(DNA_PAIRINGS)[::-1]


def none_or_bool(x):
  if x == 'None':
    return None
  if type(x) is bool:
    return x
  elif x == 'True':
    return True
  elif x == 'False':
    return False
  else:
    logging.fatal('Unrecognized bool signifier {x}'.format(**vars()))


def none_or_int(x):
  if x == 'None':
    return None
  if type(x) is int:
    return x
  else:
    return int(x)


def none_or_str(x):
  if x == 'None':
    return None
  if type(x) is str:
    return x
  else:
    return str(x)


class sgrna_target(object):
  def __init__(self, target, pam, chrom, start, end, reverse):
    self.gene = None
    self.offset = None
    self.target = str(target)
    self.pam = str(pam)
    self.chrom = str(chrom)
    self.start = int(start)
    self.end = int(end)
    self.reverse = none_or_bool(reverse)
    self.sense_strand = None
    self.specificity = 0

  @classmethod
  def from_tsv(cls, tsv, sep='\t'):
    """Alternate factory constructor from serialized sgrna_target string.

    Intended for reading back in the result of sgrna_target.__str__()
    """
    (gene,
     offset,
     target,
     pam,
     chrom,
     start,
     end,
     reverse,
     sense_strand,
     specificity) = tsv.strip().split(sep)
    t = sgrna_target(target, pam, chrom, start, end, reverse)
    t.gene = none_or_str(gene)
    t.offset = none_or_int(offset)
    t.sense_strand = none_or_bool(sense_strand)
    t.specificity= none_or_int(specificity)
    return t

  def __str__(self, sep='\t'):
    return sep.join([str(x) for x in [
        self.gene,
        self.offset,
        self.target,
        self.pam,
        self.chrom,
        self.start,
        self.end,
        self.reverse,
        self.sense_strand,
        self.specificity]])

  def id_str(self, sep=';'):
    return sep.join([str(x) for x in [
        self.target,
        self.pam,
        self.chrom,
        self.start,
        self.reverse]])

  def primer_str(self, homology, sep='\t'):
    """Generate the primer for this target.

    Args:
      homology: homology component of the primer.
      sep: optional separator for return value fields.
    Returns:
      <sep>-delimited string describing primer info.
    """
    primer_seq = self.target.lower() + homology.upper()
    return sep.join([
        self.gene,
        self.offset,
        primer_seq])

  def sequence_with_pam(self):
    """DNA sequence with trailing PAM in place."""
    return self.target + self.pam

  def build_plasmid(self, homology, rev_primer, template, output_dir):
    """Write out the ape plasmid for a target.

    Args:
      rev_primer: reverse primer to be used.
      homology: homology component of the forward primer.
      template: SeqIO.parsed genbank/ape record of template plasmid.
      output_dir: path to directory in which to place new plasmid.
    """
    template = copy.deepcopy(template)
    s = template.seq
    start = s.find(revcomp(rev_primer)) + len(rev_primer)
    end = s.find(homology)
    template.seq = s[:start] + target + s[end:]
    destination_file_name = os.path.join(
        output_dir,
        '{gene}.{offset}.sgrna.ape'.format(vars(self)))
    SeqIO.write([template], destination_file_name, 'genbank')
