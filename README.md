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

which you can call with -h to learn about its flags.  The two required arguments are

    input_fasta_genome_name [testdata/test.fna]

and

    target_regions_file [testdata/test_regions.tsv]

so to run on the sample files, you would call:

    ./build_sgrna_library.py --target_regions_file testdata/test_regions.tsv --input_fasta_genome_name testdata/doubletest.fna

If you have a gff file and want to use all the genes, you can try using (possibly with some modification)

    extract_gff_to_genes.py

Once you have created your target list with the build function, you can either search the file manually for your gene of interest, or run

    subselect_sgrna_library.py

which groups the sgrna targets by gene and calls a chosen function to choose a subset for a final library.
