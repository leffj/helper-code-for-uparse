#!/usr/bin/env python
from __future__ import division

__author__ = "Jonathan Leff"
__version__ = "0.0.1"

from os.path import split, splitext

from qiime.util import parse_command_line_parameters, get_options_lookup,\
 make_option
from cogent.parse.fasta import MinimalFastaParser
from cogent.parse.fastq import MinimalFastqParser
import re
import sys

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Print sequences per sample from a uparse formated demultiplexed fastq or fasta file"""
script_info['script_description']="""seqs per sample"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""Print sequences per sample""","""%prog -i $PWD/seqs.fq -o $PWD/seqs_counts.txt"""))
script_info['output_description']=""""""
script_info['required_options']=[\
   make_option('-i','--input_seqs_fp',type="existing_filepath",\
        help='The input sequences in fastq or fasta format')
]
script_info['optional_options']=[\
   options_lookup['output_fp'],
   make_option('-t','--file_type',\
        help='Filetype detected from file extension, but can be overridden here.'\
        'The input sequences file type (either "fastq" or "fasta")')
] 
script_info['version'] = __version__


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

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
      
    verbose = opts.verbose
    
    input_seqs_fp = opts.input_seqs_fp
    fileSize = get_file_size(input_seqs_fp)
    if opts.file_type:
        fileType = opts.file_type
    else:
        fileType = splitext(input_seqs_fp)[1].split('.')[1]
    output_fp = opts.output_fp
    
    # if the output fp isn't specified, create one
    if not output_fp:
        input_file_basename, input_file_ext = \
         splitext(split(input_seqs_fp)[1])
        output_fp = '%s_counts.txt' % (input_file_basename)

    input_seqs = open(input_seqs_fp, "U")

    output = open(output_fp, "w")

    # count the number of seqs for each unique sample
    number_seqs_bySample = {}
    printcounter = 0
    if fileType == 'fasta' or fileType == 'fa':
        for label, seq in MinimalFastaParser(input_seqs):
            matchID = re.match('^.*barcodelabel=(.*);$',label)
            sampleID = matchID.group(1)
            if sampleID in number_seqs_bySample:
                number_seqs_bySample[sampleID] += 1
            else:
                number_seqs_bySample[sampleID] = 1
            if printcounter == 1000:
                pos = input_seqs.tell()
                display_progress(pos, fileSize)
                printcounter = 0
            printcounter += 1
    elif fileType == 'fastq' or fileType == 'fq':
        for label, seq, qual in MinimalFastqParser(input_seqs,strict=False):
            matchID = re.match('^.*barcodelabel=(.*);$',label)
            sampleID = matchID.group(1)
            if sampleID in number_seqs_bySample:
                number_seqs_bySample[sampleID] += 1
            else:
                number_seqs_bySample[sampleID] = 1
            if printcounter == 1000:
                pos = input_seqs.tell()
                display_progress(pos, fileSize)
                printcounter = 0
            printcounter += 1
    else:
        print "Invalid file type"
    for key in number_seqs_bySample:
        output.write('%s\t%s\n' %(key,number_seqs_bySample[key]))

    sys.stdout.write('\n')
    input_seqs.close()
    output.close()


if __name__ == "__main__":
    main()
