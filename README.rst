sgRNA Design Scripts
====================
Author: John S. Hawkins [really@gmail.com]
Prerequisites
-------------------------------------------------

*   bowtie

**IMPORTANT NOTE**: **bowtie2** is not just "the new version", it's actually
functionally different.  This code **specifically requires** that you have
installed **bowtie**.  bowtie2's differences preclude its use for this purpose.
It's fine if you have *both* installed, provided that 'bowtie' resolves to the
non-bowtie2 version in your enviroment.

* bowtie-build (should come with bowtie)

* the Biopython library suite for python (can be installed with pip)

* the pysam library for python (can be installed with pip)

How to use this code
--------------------

Primarily you will use the script buid_sgrna_library.py.  See

::

    build_sgrna_library.py -h

for usage information.

The normal usage case is to call the script with a genbank file like so:

::

    ./build_sgrna_library.py --input_genbank_genome_name testdata/U00096.3_full_sequence.gb

Which will generate an adjacent file called

::

    testdata/U00096.3_full_sequence.targets.all.tsv

which specified all of the targets for the provided genome (the test genome, in
this case), annotated with the locus_tag, and scored for specificity (the final
column)

For bacteria we suggest using guides that

*   have a small, positive offset

*   are on the antisense strand ('anti' in the 'transdir' column)

*   have a SPECIFICITY score of 39

If a guide meeting these criteria is not available, lower specificity can be
used, but you should check for near-matches elsewhere in the genome to see if
they are likely to cause issues.  Guides on the 'sense' strand are not
recommended.  They generally have a greatly reduced, and hard to predict, level
of effect.  If reduced effect is desired, we suggest the use of
http://www.github.com/traeki/mismatch_crispri to achieve more reliable
outcomes.
