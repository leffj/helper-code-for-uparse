#!/usr/bin/env python

__author__ = "Jamie Morton, Jonathan Leff"
__email__ = "jonathan.leff@colorado.edu"
__version__ = "0.0.1"

"""
Converts an uc file to an otu table
"""

import argparse
import re
from collections import defaultdict
from collections import Counter
bar_re = re.compile("barcodelabel=(\S+);")


def main():
    parser = argparse.ArgumentParser(description=\
        'Create otu table from uc file')
    parser.add_argument('-i', '--uc_file', required=True,
        type=str, help='Input uc file')
    parser.add_argument('-o', '--output_fp', required=True,
        help='The output filepath')

    args = parser.parse_args()

    create_otu_table(args.uc_file, args.output_fp)


"""Creates OTU table using barcode and otu id"""
def create_otu_table(uc_file, out_fp):
    table = defaultdict(Counter)
    otus,barcodes = set(), set()
    
    with open(uc_file,'r') as handle:
        for ln in handle:
            ln = ln.rstrip()
            toks = ln.split('\t')
            barcode = bar_re.findall(toks[8])[0]
            otu = toks[9]
            if otu=="*": continue
            table[otu][barcode]+=1
            otus.add(otu)
            barcodes.add(barcode)
    
    # write table
    output = open(out_fp, "w")
    line = "OTUId\t%s"%('\t'.join(barcodes))
    output.write(line)
    for otu in otus:
        counts = [0]*len(barcodes)
        for idx in xrange(barcodes):
            barcode = barcodes[idx]
            counts[idx]=table[otu][barcode]
        
        line = "%s\t%s"%(otu,'\t'.join(map(str,counts)))
        output.write(line)


if __name__ == "__main__":
    main()
