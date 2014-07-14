#!/usr/bin/env python
from __future__ import division

__author__ = "Jonathan Leff"
__version__ = "0.0.1"

from os.path import split, splitext
import argparse
import re
import sys
from Bio import SeqIO
from random import random



def main():
    parser = argparse.ArgumentParser(description="Print sequences per sample from a uparse formated "\
        "demultiplexed fastq or fasta file")
    parser.add_argument('-i','--input_seqs_fp',required=True,\
        help='The input sequences in fastq or fasta format')
    parser.add_argument('-n','--number_subsample',required=True,\
        help='Specify the number of sequences per sample to subsample')
    parser.add_argument('-o', '--output_fp',
        help='The output fp')
    parser.add_argument('-t','--file_type',\
        help='Filetype detected from file extension, but can be overridden here.'\
        'The input sequences file type (either "fastq" or "fasta")')

    args = parser.parse_args()

    in_fp = args.input_seqs_fp

    nSub = args.number_subsample
    if nSub <= 0:
        raise ValueError,('number_subsample must be >0')

    # determine filetype
    if args.file_type:
        fileType = args.file_type
    else:
        fileType = splitext(in_fp)[1].split('.')[1]

    if fileType == 'fasta' or fileType == 'fa':
        ftype = 'fasta'
    elif fileType == 'fastq' or fileType == 'fq':
        ftype = 'fastq'
    
    fileSize = get_file_size(in_fp)

    # if the output fp isn't specified, create one
    output_fp = args.output_fp
    if not output_fp:
        input_file_basename, input_file_ext = \
         splitext(split(in_fp)[1])
        output_fp = '%s_sub_%s.%s' % (input_file_basename,nSub,ftype)

    # count the number of seqs for each unique sample
    in_f = open(in_fp, "U")
    number_seqs_bySample = {}
    print "Counting seqs per sample ..."
    printcounter = 0
    for seq_rec in SeqIO.parse(in_f,ftype):
        matchID = re.match('^.*barcodelabel=(.*);$',seq_rec.description)
        sampleID = matchID.group(1)
        if sampleID in number_seqs_bySample:
            number_seqs_bySample[sampleID] += 1
        else:
            number_seqs_bySample[sampleID] = 1
        if printcounter == 1000:
            pos = in_f.tell()
            display_progress(pos, fileSize)
            printcounter = 0
        printcounter += 1
    sys.stdout.write('\r%5.1f%%' % (100))
    sys.stdout.write('\n')
    in_f.close()


    # write a seq the proportion of the time (for each sample) that will give 
    # approx the right amount of seqs per sample
    print "Subsampling ..."
    in_f = open(in_fp, "U")
    out = open(output_fp, "w")
    out_seqs = []
    i = 0
    for seq_rec in SeqIO.parse(in_f,ftype):
        try:
            matchID = re.match('.*barcodelabel=(.*);',seq_rec.description)
        except AttributeError:
            matchID = re.match('(^.*?)_(.*)$',seq_rec.description)
        sampleID = matchID.group(1)
        prop = int(nSub) / number_seqs_bySample[sampleID]
        if random() <= prop:
            out_seqs.append(seq_rec)
        # for progress status
        i += 1
        if printcounter == 1000:
            pos = in_f.tell()
            display_progress(pos, fileSize)
            printcounter = 0
        printcounter += 1
    sys.stdout.write('\r%5.1f%%' % (100))
    sys.stdout.write('\n')
    in_f.close()

    SeqIO.write(out_seqs,out,'fastq')
    out.close()


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
