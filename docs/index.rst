.. sgrna_design documentation master file, created by
   sphinx-quickstart on Fri Dec  1 15:26:47 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to sgrna_design's documentation!
========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

# sgRNA Design Scripts

Author: John S. Hawkins [really@gmail.com]

## Things you will need before you can run this code

* OS X (might run on Linux/Windows, but caveat emptor

* bowtie (not bowtie2)

* bowtie-build (should come with bowtie)

* Biopython

* pysam

## How to use this code

Primarily you will use

    build_sgrna_targets.py

which you can call with -h to learn about its flags.

The normal usage case is to call the script with a genbank file like so:

    ./build_sgrna_library.py --input_genbank_genome_name testdata/U00096.3_full_sequence.gb

Which will generate a file called

    testdata/U00096.3_full_sequence.targets.all.tsv

Specifying all of the targets for the provided genome (the test genome, in this
case), annotated with the locus tag, and scored for specificity (the final
column)

For bacteria we suggest using guides that
*   have a small, positive offset
*   are on the antisense strand ('anti' in the STRAND column)
*   have a SPECIFICITY score of 39

If a guide meeting those criteria is not available, lower specificity can be
used, but you should check for near-matches elsewhere in the genome.

Alternately you can use a guide with a higher offset, or a negative offset, or
on the template strand, but you should expect lower-fold or even negligible knockdown.

(the WEAKNESS column is not relevant in the current output, ignore it)

This code will currently only work for the S. pyogenes NGG PAM sequence.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
