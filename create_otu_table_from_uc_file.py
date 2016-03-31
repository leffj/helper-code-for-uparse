#!/usr/bin/env python

__author__ = "Jamie Morton, Jonathan Leff"
__email__ = "jonathan.leff@colorado.edu"
__version__ = "0.0.1"

"""
Converts an uc file to an otu table. This is an alternative to the drive5 python script.
"""

import argparse
import sys
import re
from collections import defaultdict
from collections import Counter
bar_re = re.compile("barcodelabel=(\S+);")


def main():
    parser = argparse.ArgumentParser(description=\
        'Create otu table from uc file')
    parser.add_argument('-i', '--input_fp', required=True,
        type=str, help='Input uc file')
    parser.add_argument('-o', '--output_fp', required=True,
        help='The output filepath')

    args = parser.parse_args()

    input_fp = args.input_fp
    output_fp = args.output_fp

    create_otu_table(input_fp, output_fp)


"""Creates OTU table using barcode and otu id"""
def create_otu_table(uc_file, out_fp):
    table = defaultdict(Counter)
    otus,barcodes = set(), set()

    handle = open(uc_file,'r')
    fileSize = get_file_size(uc_file)
    printcounter = 0
    for ln in handle:
        ln = ln.rstrip()
        toks = ln.split('\t')
        try:
            barcode = bar_re.findall(toks[8])[0]
        except IndexError:
            print "Error in uc file formating. Check for spaces in sample IDs and to make sure there is a semicolon after sample IDs."
            print "First line with issue:\n%s" %ln
            break
        otu = toks[9]
        if otu=="*": continue
        table[otu][barcode]+=1
        otus.add(otu)
        barcodes.add(barcode)
        # progress
        if printcounter == 1000:
            pos = handle.tell()
            display_progress(pos, fileSize)
            printcounter = 0
        printcounter += 1
    display_progress(handle.tell(), fileSize)

    # write table
    print "\nWriting table..."
    output = open(out_fp, "w")
    line = "#OTUId\t%s\n"%('\t'.join(barcodes))
    output.write(line)
    for otu in otus:
        counts = [0]*len(barcodes)
        for idx, barcode in enumerate(barcodes):
            counts[idx]=table[otu][barcode]
        line = "%s\t%s\n"%(otu,'\t'.join(map(str,counts)))
        output.write(line)


def get_file_size(FileName):
    File = open(FileName,'U')
    Pos = File.tell()
    File.seek(0, 2)
    FileSize = File.tell()
    File.seek(Pos)
    File.close()
    return FileSize


def display_progress(position,fileSize):
    pct = (100.0*position)/fileSize
    sys.stdout.write('\r%5.1f%%' % (pct))
    sys.stdout.flush()


if __name__ == "__main__":
    main()
