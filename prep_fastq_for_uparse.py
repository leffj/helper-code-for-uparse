#!/usr/bin/env python

from __future__ import division

__author__ = "Jonathan Leff"
__version__ = "0.0.1"

import sys
from time import gmtime, strftime
from cogent.parse.fastq import MinimalFastqParser
from cogent import DNA
from itertools import izip
import gzip
from qiime.util import parse_command_line_parameters, get_options_lookup,\
 make_option


options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Prepare raw Illumina sequences for processing with UPARSE by demultiplexing."""
script_info['script_description']="""prep fastq for uparse"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""Prep fastq sequences for UPARSE""","""%prog -i $PWD/seqs.fastq -b $PWD/barcode_seqs.fastq -m $PWD/mapping_file.txt -o $PWD/seqs_for_uparse.fastq"""))
script_info['output_description']=""""""
script_info['required_options']=[\
   make_option('-i','--sequence_reads_fp',type="existing_filepath",\
        help='The sequence reads in fastq format'),
   make_option('-b','--barcode_reads_fp',type="existing_filepath",\
        help='The index (barcode) reads in fastq format'),
   make_option('-m','--mapping_file_fp',type="existing_filepath",\
        help='The the mapping file with barcodes as second column in txt format'),
   options_lookup['output_fp']
]
script_info['optional_options']=[
   make_option('-c','--rev_comp_mapping_barcodes',action='store_true',\
   	    help='should the mapping file barcodes be reverse complemented?',default=False)
] 
script_info['version'] = __version__


def get_file_size(sequence_reads_fp):
	File = open(sequence_reads_fp,'U')
	Pos = File.tell()
	File.seek(0, 2)
	FileSize = File.tell()
	File.seek(Pos)
	return FileSize

def main():
	option_parser, opts, args =\
		parse_command_line_parameters(**script_info)

	sequence_reads_fp = opts.sequence_reads_fp
	output_fp = opts.output_fp
	barcode_reads_fp = opts.barcode_reads_fp
	mapping_fp = opts.mapping_file_fp
	rc = opts.rev_comp_mapping_barcodes

	time = str(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
	sys.stdout.write('Start time: '+time+'\n')

	if sequence_reads_fp.endswith('.gz'):
		seqs = gzip.open(sequence_reads_fp,'rb')
	else:
		seqs = open(sequence_reads_fp,'U')
	if barcode_reads_fp.endswith('.gz'):
		barcodes = gzip.open(barcode_reads_fp,'rb')
	else:
		barcodes = open(barcode_reads_fp,'U')
	mapping = open(mapping_fp,'U')
	outSeqs = open(output_fp,'w')

	# create dictionary with sample IDs and barcodes
	barcodeDict = {}
	for line in mapping:
		if list(line)[0] == '#':
			continue
		sampleID = line.strip().split('\t')[0]
		barcode = line.strip().split('\t')[1]
		if rc:
			barcode = str(DNA.makeSequence(barcode).reversecomplement())
		barcodeDict[barcode] = sampleID

	# relabel sequences with sample ID as barcode
	number_seqs = 0
	number_matched = 0
	printcounter = 0
	for i,((label,seq,qual),(label2,seq2,qual2)) in enumerate(izip(MinimalFastqParser(seqs,strict=False),MinimalFastqParser(barcodes,strict=False))):
		if len(seq2) == 13:
			seq2 = seq2[:12]
		if seq2 in barcodeDict:
			number_matched += 1
			sampleID = barcodeDict[seq2]
			NewLabel = label + ";barcode=" + seq2 + ";barcodelabel=" + sampleID + ";"
			outSeqs.write('@%s\n%s\n+\n%s\n' % (NewLabel, seq, qual))
			if(printcounter == 1000):
				pct_kept = number_matched / number_seqs * 100
				sys.stdout.write('\r')
				sys.stdout.write('Seqs processed: %d Percent kept: %5.1f%%' % (number_seqs,pct_kept))
				sys.stdout.flush()
				printcounter = 0
			printcounter += 1
		else:
			if(printcounter == 1000):
				pct_kept = number_matched / number_seqs * 100
				sys.stdout.write('\r')
				sys.stdout.write('Seqs processed: %d Percent kept: %5.1f%%' % (number_seqs,pct_kept))
				sys.stdout.flush()
				printcounter = 0
			printcounter += 1
		number_seqs += 1

	time = str(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
	sys.stdout.write('\n'+'End time: '+time+'\n')

if __name__ == "__main__":
    main()
