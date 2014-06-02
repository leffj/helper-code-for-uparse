#!/usr/bin/env python

__author__ = "Jonathan Leff"
__version__ = "0.0.1"
__email__ = "jleff@symbiota.ag"

"""Filter fasta/q file to keep sequences from a list of samples"""

from Bio import SeqIO
import random
import argparse
import re
import sys
from os.path import split, splitext


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_fp', required=True,
		help='The input filepath')
	parser.add_argument('-s', '--samples_fp', required=True,
		help='Samples to keep')
	parser.add_argument('-o', '--output_fp',
		help='The output filepath')
	parser.add_argument('-t','--file_type',\
		help='Filetype detected from file extension, but can be overridden here.'\
		'The input sequences file type (either "fastq" or "fasta")')

	args = parser.parse_args()


	input_fp = args.input_fp
	fileSize = get_file_size(input_fp)

	# determine filetype
	if args.file_type:
		fileType = args.file_type
	else:
		fileType = splitext(input_fp)[1].split('.')[1]

	# if the output fp isn't specified, create one
	output_fp = args.output_fp
	if not output_fp:
		input_file_basename, input_file_ext = \
			splitext(split(input_fp)[1])
		output_fp = '%s_filtered.%s' % (input_file_basename,fileType)

	# open input and output files
	input_seqs = open(input_fp, "U")
	filterIDs = open(args.samples_fp, "U")
	output = open(output_fp, "w")

	# get list of samples
	filterIDsList = []
	for line in filterIDs:
		if line[0] == '#':
			continue
		filterIDsList.append(line.strip().split('\t')[0])

	# filter input based on list
	good_seqs = []
	printcounter = 0
	if fileType == 'fasta' or fileType == 'fa':
		ftype = 'fasta'
	elif fileType == 'fastq' or fileType == 'fq':
		ftype = 'fastq'

	seq_iterator = SeqIO.parse(input_seqs,ftype)
	for seq_record in seq_iterator:
		matchID = re.match('^.*barcodelabel=(.*);$',seq_record.description)
		sampleID = matchID.group(1)
		if sampleID in filterIDsList:
			good_seqs.append(seq_record)
		if printcounter == 1000:
			pos = input_seqs.tell()
			display_progress(pos, fileSize)
			printcounter = 0
		printcounter += 1
	SeqIO.write(good_seqs, output, ftype)



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


