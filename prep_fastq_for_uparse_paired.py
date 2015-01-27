#!/usr/bin/env python

from __future__ import division

__author__ = "Jonathan Leff"
__email__ = "jonathan.leff@colorado.edu"
__version__ = "0.0.1"

import sys
import os
from time import gmtime, strftime
from itertools import izip
import gzip
import argparse


parser = argparse.ArgumentParser(description=
	"Prepare raw Illumina sequences for processing with UPARSE by demultiplexing. This \
	version can be used for paired end reads that will later be merged but doesn't have to \
	be. NOTE: demultiplexing is based on order of sequences in the index and read files and \
	DOES NOT check headers.")
parser.add_argument('-i','--sequence_reads_fp',required=True,\
	help='The sequence reads in fastq format'),
parser.add_argument('-r','--reverse_reads_fp',required=True,\
	help='The reverse sequence reads in fastq format'),
parser.add_argument('-b','--barcode_reads_fp',required=True,\
	help='The index (barcode) reads in fastq format'),
parser.add_argument('-m','--mapping_file_fp',required=True,\
	help='The the mapping file with barcodes as second column in txt format'),
parser.add_argument('-o','--output_dir',required=True,\
	help='The output directory')
parser.add_argument('-c','--rev_comp_mapping_barcodes',action='store_true',\
	help='should the mapping file barcodes be reverse complemented?',default=False)


def main():
	args = parser.parse_args()
	
	sequence_reads_fp = args.sequence_reads_fp
	reverse_reads_fp = args.reverse_reads_fp
	output_dir = args.output_dir
	barcode_reads_fp = args.barcode_reads_fp
	mapping_fp = args.mapping_file_fp
	rc = args.rev_comp_mapping_barcodes

	time = str(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
	sys.stdout.write('Start time: '+time+'\n')

	seqs_fwd, seqs_rev, barcodes, mapping, output_seqs_fwd, output_seqs_rev = \
		setup(sequence_reads_fp, reverse_reads_fp, output_dir, barcode_reads_fp, mapping_fp)

	barcode_dictionary = create_barcode_dictionary(mapping, rc)

	relabel_to_demultiplex(barcode_dictionary, seqs_fwd, barcodes, seqs_rev, output_seqs_fwd, 
		output_seqs_rev)

	time = str(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
	sys.stdout.write('\n'+'End time: '+time+'\n')


def setup(sequence_reads_fp, reverse_reads_fp, output_dir, barcode_reads_fp, mapping_fp):
	if sequence_reads_fp.endswith('.gz'):
		seqs = gzip.open(sequence_reads_fp,'rb')
	else:
		seqs = open(sequence_reads_fp,'U')
	if reverse_reads_fp.endswith('.gz'):
		revSeqs = gzip.open(reverse_reads_fp,'rb')
	else:
		revSeqs = open(reverse_reads_fp,'U')
	if barcode_reads_fp.endswith('.gz'):
		barcodes = gzip.open(barcode_reads_fp,'rb')
	else:
		barcodes = open(barcode_reads_fp,'U')
	mapping = open(mapping_fp,'U')
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	outSeqs = open(output_dir + '/demultiplexed_seqs_1.fq','w')
	outSeqsRev = open(output_dir + '/demultiplexed_seqs_2.fq','w')
	return seqs, revSeqs, barcodes, mapping, outSeqs, outSeqsRev


def create_barcode_dictionary(mapping, rc):
	# create dictionary with sample IDs and barcodes
	barcode_dict = {}
	for line in mapping:
		if list(line)[0] == '#':
			continue
		sampleID = line.strip().split('\t')[0]
		barcode = line.strip().split('\t')[1]
		if rc:
			barcode = reverse_complement(barcode)
		barcode_dict[barcode] = sampleID
	return barcode_dict


def relabel_to_demultiplex(bc_dict, seqs, barcodes, revSeqs, seqs_out_fwd, seqs_out_rev, 
							show_progress=True):
	# relabel sequences with sample ID as barcode
	number_seqs = 0
	number_matched = 0
	printcounter = 0
	for i, ((label,seq,qual), (label2,seq2,qual2), (label3,seq3,qual3)) in \
		enumerate(izip(basic_fastq_parser(seqs), 
			basic_fastq_parser(barcodes),
			basic_fastq_parser(revSeqs))):
		if len(seq2) == 13:
			seq2 = seq2[:12]
		if seq2 in bc_dict:
			number_matched += 1
			sampleID = bc_dict[seq2]
			NewLabel = label + ";barcode=" + seq2 + ";barcodelabel=" + sampleID + ";"
			write_fastq(NewLabel, seq, qual, seqs_out_fwd)
			NewLabelRev = label3 + ";barcode=" + seq2 + ";barcodelabel=" + sampleID + ";"
			write_fastq(NewLabelRev, seq3, qual3, seqs_out_rev)
			if(printcounter == 1000):
				if show_progress:
					pct_kept = number_matched / number_seqs * 100
					sys.stdout.write('\r')
					sys.stdout.write('Seqs processed: %d Percent kept: %5.1f%%' % (number_seqs,pct_kept))
					sys.stdout.flush()
					printcounter = 0
			printcounter += 1
		else:
			if(printcounter == 1000):
				if show_progress:
					pct_kept = number_matched / number_seqs * 100
					sys.stdout.write('\r')
					sys.stdout.write('Seqs processed: %d Percent kept: %5.1f%%' % (number_seqs,pct_kept))
					sys.stdout.flush()
					printcounter = 0
			printcounter += 1
		number_seqs += 1
	if show_progress:
		pct_kept = number_matched / number_seqs * 100
		sys.stdout.write('\r')
		sys.stdout.write('Seqs processed: %d Percent kept: %5.1f%%' % (number_seqs,pct_kept))
		sys.stdout.flush()


def basic_fastq_parser(in_f):
	lineno, head, seq, qual = 0, "", "", ""
	for l in in_f:
		lineno += 1
		if lineno%4 == 1: head = l.strip()
		elif lineno%4 == 2: seq = l.strip()
		elif lineno%4 == 0:
			qual = l.strip()
			yield head, seq, qual


def write_fastq(header, seq, qual, out_f):
	out_f.write('%s\n%s\n+\n%s\n' % (header, seq, qual))

reverse_complement = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])


if __name__ == "__main__":
    main()
