helper-code-for-uparse
======================

This is a compilation of Python scripts to help with the UPARSE pipeline (http://drive5.com/uparse/), which is used to process high-throughput sequence data (i.e. from Illumina sequencing instruments) from environmental microbial communities. The end goal of this pipeline is an 'OTU table', a table that lists OTUs (aka, phylotypes, species, taxa, etc.) as rows, samples as columns, and sequence counts as values.

The two main scripts in this compilation are:

1. 'prep_fastq_for_uparse_paired.py' -- This script is useful for demultiplexing raw data from the sequencing instrument using the index reads and a mapping file. Use the '-h' option with the script for more info.

2. 'create_otu_table_from_uc_file.py' -- This script converts a .uc readmap file (output of usearch_global) to an OTU table (tab delimited). This script is much faster than the Drive5 version ('uc2otutab.py').