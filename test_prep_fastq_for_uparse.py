#!/usr/bin/env python

__author__ = "Jonathan Leff"
__email__ = "jonathan.leff@colorado.edu"
__version__ = "0.0.1"

import unittest
import prep_fastq_for_uparse_paired as testing_code
import os
import shutil

def correct_number_of_index_reads_recognized(directory, correct_no):
	seq_ct = 0
	for a_file in os.listdir(directory):
		if a_file.endswith("_1.fq"):
			for header, sequence, qualscore in basic_fastq_parser(
					open(directory + "/" + a_file, 'U')):
				seq_ct += 1
	if seq_ct == correct_no:
		return True
	else:
		return False


def correct_number_of_reads_per_sample(directory, table_w_amounts):
	# need to add
	pass


def correct_seqs_in_each_output_file(directory, dict_to_check):
	test_case = True
	for a_file in os.listdir(directory):
		# sample = os.path.splitext(a_file)[0].split('_')[0]
		for header, sequence, qualscore in basic_fastq_parser(
			open(directory + "/" + a_file, 'U')):
			sample = header.split(';')[2].split('=')[1]
			if sequence.startswith(dict_to_check[sample]):
				continue
			else:
				test_case = False
	return test_case


def basic_fastq_parser(in_f):
	lineno, head, seq, qual = 0, "", "", ""
	for l in in_f:
		lineno += 1
		if lineno%4 == 1: head = l.strip()
		elif lineno%4 == 2: seq = l.strip()
		elif lineno%4 == 0:
			qual = l.strip()
			yield head, seq, qual

class test_case(unittest.TestCase):
	# test cases
	def setUp(self):
		seqs, revSeqs, barcodes, mapping, outSeqsFwd, outSeqsRev = testing_code.setup(
			'files_for_unit_tests/test_R1.fq.gz', 
			'files_for_unit_tests/test_R2.fq.gz', 
			'files_for_unit_tests/test_out', 
			'files_for_unit_tests/test_I1.fq.gz', 
			'files_for_unit_tests/test_map.txt')
		bc_dict = testing_code.create_barcode_dictionary(mapping, rc=False)
		testing_code.relabel_to_demultiplex(bc_dict, seqs, barcodes, revSeqs, 
			outSeqsFwd, outSeqsRev, show_progress=False)

	def tearDown(self):
		shutil.rmtree('files_for_unit_tests/test_out')

	def test_correct_number_recognized_index_reads(self):
		self.assertTrue(correct_number_of_index_reads_recognized(
			'files_for_unit_tests/test_out', 8))

	def test_correct_seqs_in_each_output_file(self):
		to_check = {}
		for line in open('files_for_unit_tests/test_map.txt', 'U'):
			if line.startswith('#') == False:
				sample = line.split('\t')[0]
				sample_code = line.split('\t')[3]
				to_check[sample] = sample_code
		self.assertTrue(correct_seqs_in_each_output_file(
			'files_for_unit_tests/test_out',
			to_check))



if __name__ == '__main__':
	unittest.main()