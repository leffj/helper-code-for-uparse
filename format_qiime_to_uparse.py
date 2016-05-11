#!/usr/bin/env python

import os
import argparse

def main():
    parser = argparse.ArgumentParser(description=\
        'Take QIIME / QIITA demultiplexed sequences in fastq format and \
        reformat their headers for use in UPARSE.')
    req = parser.add_argument_group('required arguments')
    req.add_argument('-i', '--input_fp', required=True,
        type=str, help='A fastq file with headers formatted by QIIME \
        (i.e. Sample IDs first in header lines).')
    req.add_argument('-o', '--output_fp', required=True,
        help='The output file path.')

    args = parser.parse_args()

    seqs_in = open(args.input_fp, 'U')
    seqs_out = open(args.output_fp, 'w')

    for head, seq, qual in basic_fastq_parser(seqs_in):
    	sampleID = head.split('_')[0].split('@')[1]
    	new_head = head + ";barcodelabel=" + sampleID + ";"
    	write_fastq(new_head, seq, qual, seqs_out)


def basic_fastq_parser(in_f):
    lineno, head, seq, qual = 0, "", "", ""
    for l in in_f:
        lineno += 1
        if lineno % 4 == 1:
            head = l.strip()
        elif lineno % 4 == 2:
            seq = l.strip()
        elif lineno % 4 == 0:
            qual = l.strip()
            yield head, seq, qual


def write_fastq(header, seq, qual, out_f):
    out_f.write('%s\n%s\n+\n%s\n' % (header, seq, qual))


if __name__ == "__main__":
    main()