#!/usr/bin/env python

import os
import argparse
import gzip

def main():
    parser = argparse.ArgumentParser(description=\
        'Take a directory containing demultiplexed sequences in fastq files \
        (one per sample) and reformat their headers for use in UPARSE. ')
    req = parser.add_argument_group('required arguments')
    req.add_argument('-i', '--input_fp', required=True,
        type=str, help='A tab-delimited file with sample IDs in the first \
        column, corresponding filepaths for the forward reads and reverse \
        reads in the second and third columns, respectively. No header or \
        header must start with "#".')
    req.add_argument('-o', '--output_dir', required=True,
        help='The output directory.')

    args = parser.parse_args()

    # create output directory if doesn't exist
    out_dir = args.output_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    else:
        overwrite = raw_input('Output directory already exists. Do you want to \
          overwrite? (y/n)')
        if overwrite == 'y':
            pass
        else:
            raise IOError('Not overwritting. Please write to new \
            directory.')

    fwd_out = open(out_dir + '/demultiplexed_seqs_1.fq', 'w')
    rev_out = open(out_dir + '/demultiplexed_seqs_2.fq', 'w')
    for line in open(args.input_fp, 'U'):
        if line.startswith('#'):
            continue
        sampleID = line.split('\t')[0]
        fwd_fp = line.strip().split('\t')[1]
        rev_fp = line.strip().split('\t')[2]
        # for each file relabel headers and append to appropriate out file
        print 'Reformatting ' + sampleID + '...'
        if os.path.exists(fwd_fp):
            reformat(fwd_fp, fwd_out, sampleID)
        else:
            print sampleID + ' forward reads not found'
        if os.path.exists(rev_fp):
            reformat(rev_fp, rev_out, sampleID)
        else:
            print sampleID + ' reverse reads not found'

def reformat(in_fp, out_f, sampleID):
    if in_fp.endswith('.gz'):
        seqs = gzip.open(in_fp, 'rb')
    else:
        seqs = open(in_fp, 'U')
    for head, seq, qual in basic_fastq_parser(seqs):
        new_label = head + ";barcodelabel=" + sampleID + ";"
        write_fastq(new_label, seq, qual, out_f)

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
