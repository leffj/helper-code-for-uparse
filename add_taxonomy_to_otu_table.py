#!/usr/bin/env python

__author__ = "Jonathan Leff"
__email__ = "jonathan.leff@colorado.edu"
__version__ = "0.0.1"

"""
Add taxonomy strings to the last column of a text formatted OTU table based
on corresponding entries in a database (such as Greengenes)
"""

import argparse
import gzip
import subprocess
import os
from collections import defaultdict
from collections import Counter


def main():
    parser = argparse.ArgumentParser(description=\
        'Add taxonomy strings to the last column of a text formatted OTU table. \
        will overwrite original file.')
    req = parser.add_argument_group('required arguments')
    req.add_argument('-i', '--input_fp', required=True,
        type=str, help='The filepath of a tab-delimited OTU table')
    req.add_argument('-d', '--taxonomy_db_fp', required=True,
        help='The taxonomy database filepath. So far, only tested with \
        greengenes clustered at 97%% similarity. Format is tab-delimited \
        with the OTU IDs in the first column and the taxonomy strings in the \
        second.')
    # req.add_argument('-o', '--output_fp', required=True,
    #     help='The output filepath.')

    args = parser.parse_args()

    add_taxonomy_to_otu_table(args.input_fp, args.taxonomy_db_fp)


def add_taxonomy_to_otu_table(otu_table_fp, taxonomy_fp):
    otu_table = open(otu_table_fp, 'U')
    taxonomy = open(taxonomy_fp, 'U')
    taxonomy_dict = {}
    table_out = []
    for line in taxonomy:
        line_split = line.strip().split('\t')
        taxonomy_dict[line_split[0]] = line_split[1]
    for line in otu_table:
        line_split = line.strip().split('\t')
        if line_split[0].startswith('#'):
            table_out.append(line.strip() + "\ttaxonomy\n")
        else:
            otu_tax = taxonomy_dict[line_split[0]]
            table_out.append("%s\t%s\n" %(line.strip(), otu_tax))
    otu_table.close()
    output = open(otu_table_fp, 'w')
    output.write("".join(table_out))
    output.close()

if __name__ == "__main__":
    main()
