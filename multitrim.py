#!/usr/bin/env python3

import sys
import os
import subprocess
import tempfile
import argparse
import multiprocessing
import re
import shutil
from datetime import datetime

#Reads a file with adapters and uses them as the starting set for adapter identification. By default, uses the current MiGA adapter list as of Feb. 23, 2021
def read_adapters(adapters_fasta):

	cleanup = False
	if adapters_fasta == "internal":
		adapters, adapters_fasta = generate_adapters_temporary_file()

		cleanup = True
	else:
		adapters = {}
		current_seq = ""
		current_id = ""

		adapt = open(adapters_fasta, "r")

		for line in adapt:
			if line.startswith(">"):
				if len(current_seq) > 0:
					adapters[current_id] = current_seq
				current_id = line.strip()[1:]
				current_seq = ""
			else:
				current_seq += line.strip()

		adapters[current_id] = current_seq

		adapt.close()

	return adapters, adapters_fasta, cleanup

#Only contains adapters we already recognize as part of a kit. It will need updated as new ones may be added.
def family_detection(adapter_seqs):
	#Currently acceptable fams:
	'''
	singleend
	pairedend
	dpnII
	smallrna
	multiplex
	pcr
	dpnIIgex
	otherrna
	trueseq
	rnapcr
	trueseq2
	nextera
	cre-loxp
	truseq1
	pcr_primer
	nextera_junction
	'''

	#There are some repeats in adapters. All are added - this meant to make the program as conservative as possible.
	fam_to_id_to_seq = {}
	#MiGA adapters
	fam_to_id_to_seq['singleend'] = {'Illumina_Single_End_Apapter_1': 'ACACTCTTTCCCTACACGACGCTGTTCCATCT', 'Illumina_Single_End_Apapter_2': 'CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT', 'Illumina_Single_End_PCR_Primer_1': 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'Illumina_Single_End_PCR_Primer_2': 'CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT', 'Illumina_Single_End_Sequencing_Primer': 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'}
	fam_to_id_to_seq['pairedend'] = {'Illumina_Paired_End_Adapter_1': 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'Illumina_Paired_End_Adapter_2': 'CTCGGCATTCCTGCTGAACCGCTCTTCCGATCT', 'Illumina_Paried_End_PCR_Primer_1': 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'Illumina_Paired_End_PCR_Primer_2': 'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT', 'Illumina_Paried_End_Sequencing_Primer_1': 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'Illumina_Paired_End_Sequencing_Primer_2': 'CGGTCTCGGCATTCCTACTGAACCGCTCTTCCGATCT'}
	fam_to_id_to_seq['dpnII'] = {'Illumina_DpnII_expression_Adapter_1': 'ACAGGTTCAGAGTTCTACAGTCCGAC', 'Illumina_DpnII_expression_Adapter_2': 'CAAGCAGAAGACGGCATACGA', 'Illumina_DpnII_expression_PCR_Primer_1': 'CAAGCAGAAGACGGCATACGA', 'Illumina_DpnII_expression_PCR_Primer_2': 'AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA', 'Illumina_DpnII_expression_Sequencing_Primer': 'CGACAGGTTCAGAGTTCTACAGTCCGACGATC', 'Illumina_NlaIII_expression_Adapter_1': 'ACAGGTTCAGAGTTCTACAGTCCGACATG', 'Illumina_NlaIII_expression_Adapter_2': 'CAAGCAGAAGACGGCATACGA', 'Illumina_NlaIII_expression_PCR_Primer_1': 'CAAGCAGAAGACGGCATACGA', 'Illumina_NlaIII_expression_PCR_Primer_2': 'AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA', 'Illumina_NlaIII_expression_Sequencing_Primer': 'CCGACAGGTTCAGAGTTCTACAGTCCGACATG'}
	fam_to_id_to_seq['smallrna'] = {'Illumina_Small_RNA_Adapter_1': 'GTTCAGAGTTCTACAGTCCGACGATC', 'Illumina_Small_RNA_Adapter_2': 'TCGTATGCCGTCTTCTGCTTGT', 'Illumina_Small_RNA_RT_Primer': 'CAAGCAGAAGACGGCATACGA', 'Illumina_Small_RNA_PCR_Primer_1': 'CAAGCAGAAGACGGCATACGA', 'Illumina_Small_RNA_PCR_Primer_2': 'AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA', 'Illumina_Small_RNA_Sequencing_Primer': 'CGACAGGTTCAGAGTTCTACAGTCCGACGATC'}
	fam_to_id_to_seq['multiplex'] = {'Illumina_Multiplexing_Adapter_1': 'GATCGGAAGAGCACACGTCT', 'Illumina_Multiplexing_Adapter_2': 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'Illumina_Multiplexing_PCR_Primer_1.01': 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'Illumina_Multiplexing_PCR_Primer_2.01': 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT', 'Illumina_Multiplexing_Read1_Sequencing_Primer': 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'Illumina_Multiplexing_Index_Sequencing_Primer': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC', 'Illumina_Multiplexing_Read2_Sequencing_Primer': 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'}
	fam_to_id_to_seq['pcr'] = {'Illumina_PCR_Primer_Index_1': 'CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC', 'Illumina_PCR_Primer_Index_2': 'CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC', 'Illumina_PCR_Primer_Index_3': 'CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC', 'Illumina_PCR_Primer_Index_4': 'CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC', 'Illumina_PCR_Primer_Index_5': 'CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC', 'Illumina_PCR_Primer_Index_6': 'CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC', 'Illumina_PCR_Primer_Index_7': 'CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC', 'Illumina_PCR_Primer_Index_8': 'CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC', 'Illumina_PCR_Primer_Index_9': 'CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC', 'Illumina_PCR_Primer_Index_10': 'CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC', 'Illumina_PCR_Primer_Index_11': 'CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC', 'Illumina_PCR_Primer_Index_12': 'CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC'}
	fam_to_id_to_seq['dpnIIgex'] = {'Illumina_DpnII_Gex_Adapter_1': 'GATCGTCGGACTGTAGAACTCTGAAC', 'Illumina_DpnII_Gex_Adapter_1.01': 'ACAGGTTCAGAGTTCTACAGTCCGAC', 'Illumina_DpnII_Gex_Adapter_2': 'CAAGCAGAAGACGGCATACGA', 'Illumina_DpnII_Gex_Adapter_2.01': 'TCGTATGCCGTCTTCTGCTTG', 'Illumina_DpnII_Gex_PCR_Primer_1': 'CAAGCAGAAGACGGCATACGA', 'Illumina_DpnII_Gex_PCR_Primer_2': 'AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA', 'Illumina_DpnII_Gex_Sequencing_Primer': 'CGACAGGTTCAGAGTTCTACAGTCCGACGATC', 'Illumina_NlaIII_Gex_Adapter_1.01': 'TCGGACTGTAGAACTCTGAAC', 'Illumina_NlaIII_Gex_Adapter_1.02': 'ACAGGTTCAGAGTTCTACAGTCCGACATG', 'Illumina_NlaIII_Gex_Adapter_2.01': 'CAAGCAGAAGACGGCATACGA', 'Illumina_NlaIII_Gex_Adapter_2.02': 'TCGTATGCCGTCTTCTGCTTG', 'Illumina_NlaIII_Gex_PCR_Primer_1': 'CAAGCAGAAGACGGCATACGA', 'Illumina_NlaIII_Gex_PCR_Primer_2': 'AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA', 'Illumina_NlaIII_Gex_Sequencing_Primer': 'CCGACAGGTTCAGAGTTCTACAGTCCGACATG'}
	fam_to_id_to_seq['otherrna'] = {'Illumina_5p_RNA_Adapter': 'GTTCAGAGTTCTACAGTCCGACGATC', 'Illumina_RNA_Adapter1': 'TCGTATGCCGTCTTCTGCTTGT', 'Illumina_Small_RNA_3p_Adapter_1': 'ATCTCGTATGCCGTCTTCTGCTTG'}
	fam_to_id_to_seq['trueseq'] = {'TruSeq_Universal_Adapter': 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'TruSeq_Adapter_Index_1': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG', 'TruSeq_Adapter_Index_2': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG', 'TruSeq_Adapter_Index_3': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG', 'TruSeq_Adapter_Index_4': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG', 'TruSeq_Adapter_Index_5': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG', 'TruSeq_Adapter_Index_6': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG', 'TruSeq_Adapter_Index_7': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG', 'TruSeq_Adapter_Index_8': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG', 'TruSeq_Adapter_Index_9': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG', 'TruSeq_Adapter_Index_10': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG', 'TruSeq_Adapter_Index_11': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG', 'TruSeq_Adapter_Index_12': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG'}
	fam_to_id_to_seq['rnapcr'] = {'Illumina_RNA_RT_Primer': 'GCCTTGGCACCCGAGAATTCCA', 'Illumina_RNA_PCR_Primer': 'AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA', 'RNA_PCR_Primer_Index_1': 'CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_2': 'CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_3': 'CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_4': 'CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_5': 'CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_6': 'CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_7': 'CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_8': 'CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_9': 'CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_10': 'CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_11': 'CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_12': 'CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_13': 'CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_14': 'CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_15': 'CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_16': 'CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_17': 'CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_18': 'CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_19': 'CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_20': 'CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_21': 'CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_22': 'CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_23': 'CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_24': 'CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_25': 'CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_26': 'CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_27': 'CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_28': 'CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_29': 'CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_30': 'CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_31': 'CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_32': 'CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_33': 'CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_34': 'CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_35': 'CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_36': 'CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_37': 'CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_38': 'CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_39': 'CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_40': 'CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_41': 'CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_42': 'CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_43': 'CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_44': 'CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_45': 'CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_46': 'CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_47': 'CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA', 'RNA_PCR_Primer_Index_48': 'CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'}
	fam_to_id_to_seq['abi'] = {'ABI_Dynabead_EcoP_Oligo': 'CTGATCTAGAGGTACCGGATCCCAGCAGT', 'ABI_Solid3_Adapter_A': 'CTGCCCCGGGTTCCTCATTCTCTCAGCAGCATG', 'ABI_Solid3_Adapter_B': 'CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT', 'ABI_Solid3_5_AMP_Primer': 'CCACTACGCCTCCGCTTTCCTCTCTATG', 'ABI_Solid3_3_AMP_Primer': 'CTGCCCCGGGTTCCTCATTCT', 'ABI_Solid3_EF1_alpha_Sense_Primer': 'CATGTGTGTTGAGAGCTTC', 'ABI_Solid3_EF1_alpha_Antisense_Primer': 'GAAAACCAAAGTGGTCCAC', 'ABI_Solid3_GAPDH_Forward_Primer': 'TTAGCACCCCTGGCCAAGG', 'ABI_Solid3_GAPDH_Reverse_Primer': 'CTTACTCCTTGGAGGCCATG'}
	fam_to_id_to_seq['trueseq2'] = {'TruSeq2_SE': 'AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG', 'TruSeq2_PE_f': 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT', 'TruSeq2_PE_r': 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG', 'TruSeq3_IndexedAdapter': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', 'TruSeq3_UniversalAdapter': 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'}
	fam_to_id_to_seq['nextera'] = {'Nextera_PE_PrefixNX/1': 'AGATGTGTATAAGAGACAG', 'Nextera_PE_PrefixNX/2': 'AGATGTGTATAAGAGACAG', 'Nextera_PE_Trans1': 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG', 'Nextera_PE_Trans1_rc': 'CTGTCTCTTATACACATCTGACGCTGCCGACGA', 'Nextera_PE_Trans2': 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG', 'Nextera_PE_Trans2_rc': 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'}
	#FaQCs adapters for safety
	fam_to_id_to_seq['cre-loxp'] = {'cre-loxp-forward' : 'TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG', 'cre-loxp-reverse' : 'AGCATATTGAAGCATATTACATACGATATGCTTCAATAATGC'}
	fam_to_id_to_seq['truseq1'] = {'TruSeq-adapter-1' : 'GGGGTAGTGTGGATCCTCCTCTAGGCAGTTGGGTTATTCTAGAAGCAGATGTGTTGGCTGTTTCTGAAACTCTGGAAAA', 'TruSeq-adapter-3' : 'CAACAGCCGGTCAAAACATCTGGAGGGTAAGCCATAAACACCTCAACAGAAAA'}
	fam_to_id_to_seq['pcr_primer'] = {'PCR-primer-1' : 'CGATAACTTCGTATAATGTATGCTATACGAAGTTATTACG', 'PCR-primer-2' : 'GCATAACTTCGTATAGCATACATTATACGAAGTTATACGA'}
	fam_to_id_to_seq['nextera_junction'] = {'Nextera-junction-adapter-1' : 'CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG'}
	fam_to_id_to_seq['Nextera-primer-adapter'] = {'Nextera-primer-adapter-1' : 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC', 'Nextera-primer-adapter-2' : 'GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'}

	detected_fams = []
	#for each user sequence
	for seq in adapter_seqs.values():
		#for each family
		for family in fam_to_id_to_seq:
			#Check if the sequence appears as a value in the current family; add it to the detected family list if it's not already there.
			if seq in fam_to_id_to_seq[family].values() and family not in detected_fams:
				detected_fams.append(family)

	#If a family was detected, add ALL sequences from that family to the final list, except the ones that the user already supplied.
	for fam in detected_fams:
		for id in fam_to_id_to_seq[fam]:
			sequence = fam_to_id_to_seq[fam][id]
			#User seqs came in with adapter_seqs, so we want to skip adding the preset one; it would be redundant, change names. Just add the others.
			if sequence not in adapter_seqs.values():
				adapter_seqs[id] = sequence

	#add '>' to the start of each seq.
	easy_print = {}

	for id in adapter_seqs:
		easy_print[">"+id] = adapter_seqs[id]

	return easy_print

#This contains code which generates a complete list of illumina adapters from scratch
def generate_adapters_temporary_file():

	#print("Preparing adapter file for you.")
	adapters_dict = {}

	'''
	I identify the adapter families here with comments. Any adapter recognized in one of these during preprocessing will include
	all of the members of its family in final, e.g. seeing Illumina_Single_End_Apapter_1 will include the following:
	Illumina_Single_End_Apapter_1, Illumina_Single_End_Apapter_2, Illumina_Single_End_PCR_Primer_1, Illumina_Single_End_PCR_Primer_2, and Illumina_Single_End_Sequencing_Primer
	in the final filtering fasta
	'''

	#Single end family
	adapters_dict["Illumina_Single_End_Apapter_1"] = "ACACTCTTTCCCTACACGACGCTGTTCCATCT"
	adapters_dict["Illumina_Single_End_Apapter_2"] = "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"
	adapters_dict["Illumina_Single_End_PCR_Primer_1"] = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapters_dict["Illumina_Single_End_PCR_Primer_2"] = "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"
	adapters_dict["Illumina_Single_End_Sequencing_Primer"] = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"

	#Paired end family
	adapters_dict["Illumina_Paired_End_Adapter_1"] = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapters_dict["Illumina_Paired_End_Adapter_2"] = "CTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"
	adapters_dict["Illumina_Paried_End_PCR_Primer_1"] = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapters_dict["Illumina_Paired_End_PCR_Primer_2"] = "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"
	adapters_dict["Illumina_Paried_End_Sequencing_Primer_1"] = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapters_dict["Illumina_Paired_End_Sequencing_Primer_2"] = "CGGTCTCGGCATTCCTACTGAACCGCTCTTCCGATCT"

	#DpnII family
	adapters_dict["Illumina_DpnII_expression_Adapter_1"] = "ACAGGTTCAGAGTTCTACAGTCCGAC"
	adapters_dict["Illumina_DpnII_expression_Adapter_2"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict["Illumina_DpnII_expression_PCR_Primer_1"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict["Illumina_DpnII_expression_PCR_Primer_2"] = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"
	adapters_dict["Illumina_DpnII_expression_Sequencing_Primer"] = "CGACAGGTTCAGAGTTCTACAGTCCGACGATC"
	adapters_dict["Illumina_NlaIII_expression_Adapter_1"] = "ACAGGTTCAGAGTTCTACAGTCCGACATG"
	adapters_dict["Illumina_NlaIII_expression_Adapter_2"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict["Illumina_NlaIII_expression_PCR_Primer_1"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict["Illumina_NlaIII_expression_PCR_Primer_2"] = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"
	adapters_dict["Illumina_NlaIII_expression_Sequencing_Primer"] = "CCGACAGGTTCAGAGTTCTACAGTCCGACATG"

	#Small RNA family
	adapters_dict["Illumina_Small_RNA_Adapter_1"] = "GTTCAGAGTTCTACAGTCCGACGATC"
	adapters_dict["Illumina_Small_RNA_Adapter_2"] = "TCGTATGCCGTCTTCTGCTTGT"
	adapters_dict["Illumina_Small_RNA_RT_Primer"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict["Illumina_Small_RNA_PCR_Primer_1"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict["Illumina_Small_RNA_PCR_Primer_2"] = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"
	adapters_dict["Illumina_Small_RNA_Sequencing_Primer"] = "CGACAGGTTCAGAGTTCTACAGTCCGACGATC"


	#Multiplexing Family
	adapters_dict["Illumina_Multiplexing_Adapter_1"] = "GATCGGAAGAGCACACGTCT"
	adapters_dict["Illumina_Multiplexing_Adapter_2"] = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapters_dict["Illumina_Multiplexing_PCR_Primer_1.01"] = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapters_dict["Illumina_Multiplexing_PCR_Primer_2.01"] = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
	adapters_dict["Illumina_Multiplexing_Read1_Sequencing_Primer"] = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapters_dict["Illumina_Multiplexing_Index_Sequencing_Primer"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
	adapters_dict["Illumina_Multiplexing_Read2_Sequencing_Primer"] = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"


	#PCR primer family
	adapters_dict["Illumina_PCR_Primer_Index_1"] = "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC"
	adapters_dict["Illumina_PCR_Primer_Index_2"] = "CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC"
	adapters_dict["Illumina_PCR_Primer_Index_3"] = "CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC"
	adapters_dict["Illumina_PCR_Primer_Index_4"] = "CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC"
	adapters_dict["Illumina_PCR_Primer_Index_5"] = "CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC"
	adapters_dict["Illumina_PCR_Primer_Index_6"] = "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC"
	adapters_dict["Illumina_PCR_Primer_Index_7"] = "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC"
	adapters_dict["Illumina_PCR_Primer_Index_8"] = "CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC"
	adapters_dict["Illumina_PCR_Primer_Index_9"] = "CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC"
	adapters_dict["Illumina_PCR_Primer_Index_10"] = "CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC"
	adapters_dict["Illumina_PCR_Primer_Index_11"] = "CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC"
	adapters_dict["Illumina_PCR_Primer_Index_12"] = "CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC"


	#DpnII Gex family
	adapters_dict["Illumina_DpnII_Gex_Adapter_1"] = "GATCGTCGGACTGTAGAACTCTGAAC"
	adapters_dict["Illumina_DpnII_Gex_Adapter_1.01"] = "ACAGGTTCAGAGTTCTACAGTCCGAC"
	adapters_dict["Illumina_DpnII_Gex_Adapter_2"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict["Illumina_DpnII_Gex_Adapter_2.01"] = "TCGTATGCCGTCTTCTGCTTG"
	adapters_dict["Illumina_DpnII_Gex_PCR_Primer_1"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict["Illumina_DpnII_Gex_PCR_Primer_2"] = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"
	adapters_dict["Illumina_DpnII_Gex_Sequencing_Primer"] = "CGACAGGTTCAGAGTTCTACAGTCCGACGATC"
	adapters_dict["Illumina_NlaIII_Gex_Adapter_1.01"] = "TCGGACTGTAGAACTCTGAAC"
	adapters_dict["Illumina_NlaIII_Gex_Adapter_1.02"] = "ACAGGTTCAGAGTTCTACAGTCCGACATG"
	adapters_dict["Illumina_NlaIII_Gex_Adapter_2.01"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict["Illumina_NlaIII_Gex_Adapter_2.02"] = "TCGTATGCCGTCTTCTGCTTG"
	adapters_dict["Illumina_NlaIII_Gex_PCR_Primer_1"] = "CAAGCAGAAGACGGCATACGA"
	adapters_dict["Illumina_NlaIII_Gex_PCR_Primer_2"] = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"
	adapters_dict["Illumina_NlaIII_Gex_Sequencing_Primer"] = "CCGACAGGTTCAGAGTTCTACAGTCCGACATG"

	#Other RNA family
	adapters_dict["Illumina_5p_RNA_Adapter"] = "GTTCAGAGTTCTACAGTCCGACGATC"
	adapters_dict["Illumina_RNA_Adapter1"] = "TCGTATGCCGTCTTCTGCTTGT"
	adapters_dict["Illumina_Small_RNA_3p_Adapter_1"] = "ATCTCGTATGCCGTCTTCTGCTTG"

	#TrueSeq family
	adapters_dict["TruSeq_Universal_Adapter"] = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapters_dict["TruSeq_Adapter_Index_1"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict["TruSeq_Adapter_Index_2"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict["TruSeq_Adapter_Index_3"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict["TruSeq_Adapter_Index_4"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict["TruSeq_Adapter_Index_5"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict["TruSeq_Adapter_Index_6"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict["TruSeq_Adapter_Index_7"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict["TruSeq_Adapter_Index_8"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict["TruSeq_Adapter_Index_9"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict["TruSeq_Adapter_Index_10"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict["TruSeq_Adapter_Index_11"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict["TruSeq_Adapter_Index_12"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG"

	#RNA PCR family
	adapters_dict["Illumina_RNA_RT_Primer"] = "GCCTTGGCACCCGAGAATTCCA"
	adapters_dict["Illumina_RNA_PCR_Primer"] = "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA"
	adapters_dict["RNA_PCR_Primer_Index_1"] = "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_2"] = "CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_3"] = "CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_4"] = "CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_5"] = "CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_6"] = "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_7"] = "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_8"] = "CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_9"] = "CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_10"] = "CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_11"] = "CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_12"] = "CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_13"] = "CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_14"] = "CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_15"] = "CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_16"] = "CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_17"] = "CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_18"] = "CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_19"] = "CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_20"] = "CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_21"] = "CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_22"] = "CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_23"] = "CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_24"] = "CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_25"] = "CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_26"] = "CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_27"] = "CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_28"] = "CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_29"] = "CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_30"] = "CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_31"] = "CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_32"] = "CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_33"] = "CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_34"] = "CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_35"] = "CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_36"] = "CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_37"] = "CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_38"] = "CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_39"] = "CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_40"] = "CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_41"] = "CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_42"] = "CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_43"] = "CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_44"] = "CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_45"] = "CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_46"] = "CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_47"] = "CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
	adapters_dict["RNA_PCR_Primer_Index_48"] = "CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"

	#ABI family
	adapters_dict["ABI_Dynabead_EcoP_Oligo"] = "CTGATCTAGAGGTACCGGATCCCAGCAGT"
	adapters_dict["ABI_Solid3_Adapter_A"] = "CTGCCCCGGGTTCCTCATTCTCTCAGCAGCATG"
	adapters_dict["ABI_Solid3_Adapter_B"] = "CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT"
	adapters_dict["ABI_Solid3_5_AMP_Primer"] = "CCACTACGCCTCCGCTTTCCTCTCTATG"
	adapters_dict["ABI_Solid3_3_AMP_Primer"] = "CTGCCCCGGGTTCCTCATTCT"
	adapters_dict["ABI_Solid3_EF1_alpha_Sense_Primer"] = "CATGTGTGTTGAGAGCTTC"
	adapters_dict["ABI_Solid3_EF1_alpha_Antisense_Primer"] = "GAAAACCAAAGTGGTCCAC"
	adapters_dict["ABI_Solid3_GAPDH_Forward_Primer"] = "TTAGCACCCCTGGCCAAGG"
	adapters_dict["ABI_Solid3_GAPDH_Reverse_Primer"] = "CTTACTCCTTGGAGGCCATG"

	#TrueSeq2 family
	adapters_dict["TruSeq2_SE"] = "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"
	adapters_dict["TruSeq2_PE_f"] = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
	adapters_dict["TruSeq2_PE_r"] = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"
	adapters_dict["TruSeq3_IndexedAdapter"] = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
	adapters_dict["TruSeq3_UniversalAdapter"] = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"

	#Nextera Family
	adapters_dict["Nextera_PE_PrefixNX/1"] = "AGATGTGTATAAGAGACAG"
	adapters_dict["Nextera_PE_PrefixNX/2"] = "AGATGTGTATAAGAGACAG"
	adapters_dict["Nextera_PE_Trans1"] = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
	adapters_dict["Nextera_PE_Trans1_rc"] = "CTGTCTCTTATACACATCTGACGCTGCCGACGA"
	adapters_dict["Nextera_PE_Trans2"] = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
	adapters_dict["Nextera_PE_Trans2_rc"] = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"

	all_adapters = tempfile.NamedTemporaryFile(mode = "w", delete = False)

	for adapt in adapters_dict:
		print(">"+adapt, file = all_adapters)
		print(adapters_dict[adapt], file = all_adapters)

	name = all_adapters.name
	all_adapters.close()

	return adapters_dict, name

#FaCQs supports external adapter sequences, but has no option to EXCLUDE its own internal adapters while doing so.
#This function returns a dict of ID:adapter for the FaQCs internal sequences so that Multitrim doesn't break should a FaQCs adapter
#appear in parse_adapters
def faqcs_internal_adapters():
	adapters_dict = {}
	#Below are the adapters present in FaQCs by default.

	#Cre-loxp family
	adapters_dict["cre-loxp-forward"] = "TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG"
	adapters_dict["cre-loxp-reverse"] = "AGCATATTGAAGCATATTACATACGATATGCTTCAATAATGC"

	#TruSeq 1 family
	adapters_dict["TruSeq-adapter-1"] = "GGGGTAGTGTGGATCCTCCTCTAGGCAGTTGGGTTATTCTAGAAGCAGATGTGTTGGCTGTTTCTGAAACTCTGGAAAA"
	adapters_dict["TruSeq-adapter-3"] = "CAACAGCCGGTCAAAACATCTGGAGGGTAAGCCATAAACACCTCAACAGAAAA"

	#PCR primers
	adapters_dict["PCR-primer-1"] = "CGATAACTTCGTATAATGTATGCTATACGAAGTTATTACG"
	adapters_dict["PCR-primer-2"] = "GCATAACTTCGTATAGCATACATTATACGAAGTTATACGA"

	#Nextera Junction family
	adapters_dict["Nextera-junction-adapter-1"] = "CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG"

	#Nextera-primer-adapter family; these are copies of earlier adapters in this list, but I want to make sure they're detectable since they're internal to FaQCs
	adapters_dict["Nextera-primer-adapter-1"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
	adapters_dict["Nextera-primer-adapter-2"] = "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

	return(adapters_dict)

#Get file names right up front for ease of use
def names_pe(forward, reverse, outdir = ".", prefix = ""):
	forward_basename = os.path.basename(os.path.normpath(forward))
	if forward_basename.endswith(".gz"):
		forward_basename = forward_basename[:-3]
	forward_basename = os.path.splitext(forward_basename)[0]

	reverse_basename = os.path.basename(os.path.normpath(reverse))
	if reverse_basename.endswith(".gz"):
		reverse_basename = reverse_basename[:-3]
	reverse_basename = os.path.splitext(reverse_basename)[0]

	pre_qc_f = outdir + "/" + prefix + "1.pre_trim_QC_" + forward_basename
	pre_qc_r = outdir + "/" + prefix + "2.pre_trim_QC_" + reverse_basename

	post_qc_f = outdir + "/" + prefix + "1.post_trim_QC_" + forward_basename
	post_qc_r = outdir + "/" + prefix + "2.post_trim_QC_" + reverse_basename

	post_trim_reads_f = outdir + "/" + prefix + "1.post_trim_" + forward_basename + ".fq"
	post_trim_reads_r = outdir + "/" + prefix + "2.post_trim_" + reverse_basename + ".fq"

	return pre_qc_f, pre_qc_r, post_qc_f, post_qc_r, post_trim_reads_f, post_trim_reads_r

#Get file names right up front for ease of use
def names_se(reads, outdir = ".", prefix = ""):
	base_name = os.path.basename(os.path.normpath(reads))
	if base_name.endswith(".gz"):
		base_name = base_name[:-3]
	base_name = os.path.splitext(base_name)[0]

	pre_qc = outdir + "/" + prefix + "unpaired.pre_trim_QC_" + base_name
	post_qc = outdir + "/" + prefix + "unpaired.post_trim_QC_" + base_name
	post_trim_reads = outdir + "/" + prefix + "unpaired.post_trim_" + base_name + ".fq"

	return pre_qc, post_qc, post_trim_reads

#DSRC needs its own. Whoops.
def do_falco(read_name_tool):
	'''
	Falco does not support naming files, but does support selecting output directory.

	As we are possibly generating multiple falco reports simultaneously,
	we get around this issue by generating the generically named files in a temp dir
	and then move the results to the final location with an appropriate rename.
	'''

	#temp directory
	loc = tempfile.mkdtemp()

	reads = read_name_tool[0]
	output_name = read_name_tool[1]
	falco_path = read_name_tool[2]

	#falco command
	command = [falco_path, "--quiet", "-o", loc, reads]
	#command = [falco_path, "-o", loc, reads]


	#run the command
	#Working perfectly, the falco call should not produce any output. Until falco has bugs patched, it's not working perfectly
	#subprocess.call(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	subprocess.call(command)


	#move the results and rename
	#I'm just gonna move the html.
	#shutil.move(loc+"/fastqc_data.txt",    output_name + ".data.txt")
	shutil.move(loc+"/fastqc_report.html", output_name + ".html")

	#Cleanup
	shutil.rmtree(loc)

	return None

#do all QC at once - now old
def falco_qc_pe(pre_trim_reads_f, pre_trim_reads_r, post_trim_reads_f, post_trim_reads_r, pre_name_f, post_name_f, pre_name_r, post_name_r, threads, falco_binary):
	pre_forward = [pre_trim_reads_f, pre_name_f, falco_binary]
	pre_reverse = [pre_trim_reads_r, pre_name_r, falco_binary]
	post_forward = [post_trim_reads_f+".gz", post_name_f, falco_binary]
	post_reverse = [post_trim_reads_r+".gz", post_name_r, falco_binary]

	commands = [pre_forward, pre_reverse, post_forward, post_reverse]

	print("Generating QC reports.")

	pool = multiprocessing.Pool(min(4, threads))

	pool.map(do_falco, commands)

	pool.close()

#do all QC at once - now old
def falco_qc_se(pre_trim_reads, post_trim_reads, pre_name, post_name, threads, falco_binary):
	pre = [pre_trim_reads, pre_name, falco_binary]
	post = [post_trim_reads+".gz", post_name, falco_binary]

	commands = [pre, post]

	print("Generating QC reports.")

	pool = multiprocessing.Pool(min(2, threads))

	pool.map(do_falco, commands)

	pool.close()

def do_seqtk(read_tool):
	sample = read_tool[0]
	seqtk_path = read_tool[1]

	print("Subsampling:", sample)

	#-s 100 specifies seed as 100. The number chosen is arbitrary, and I only spcify it so that results are deterministic and reproducible.
	command = [seqtk_path, "sample", "-s", "100", sample, "100000"]

	temp = tempfile.NamedTemporaryFile("w", delete=False)

	ps = subprocess.run(command, stdout=subprocess.PIPE, universal_newlines = True)
	temp.write(ps.stdout)

	name = temp.name

	temp.close()

	return name

#Subsample reads; identify adapters with FaQCs
def adapter_identification_pe(artificial_artifacts, seqtk_binary, faqcs_binary, forward = "", reverse = "", threads = 1, output = ".", minimum_presence = 0.1, prefix = "", phred_fmt = "33"):
	#seqtk forward and reverse
	subsample_f = [forward, seqtk_binary]
	subsample_r = [reverse, seqtk_binary]

	seqtk_commands = [subsample_f, subsample_r]

	pool = multiprocessing.Pool(min(2, threads))

	seqtk_samples = pool.map(do_seqtk, seqtk_commands)

	pool.close()

	#FaQCs PE with adapter file
	faqcs_subset_command = [faqcs_binary, "-t", str(threads), "--qc_only", "-d", output, "--artifactFile", artificial_artifacts, "--ascii", phred_fmt]

	#proper naming
	if prefix != "":
		faqcs_subset_command.append("--prefix")
		faqcs_subset_command.append(prefix + "Subsample_Adapter_Detection")
		pdf_name = prefix + "Subsample_Adapter_Detection_qc_report.pdf"
	else :
		faqcs_subset_command.append("--prefix")
		faqcs_subset_command.append("Subsample_Adapter_Detection")
		pdf_name = "Subsample_Adapter_Detection_qc_report.pdf"

	#forward strand
	faqcs_subset_command.append("-1")
	faqcs_subset_command.append(seqtk_samples[0])

	#reverse strand
	faqcs_subset_command.append("-2")
	faqcs_subset_command.append(seqtk_samples[1])

	print("Detecting adapters now... ", end = "")
	ps = subprocess.Popen(faqcs_subset_command)
	ps.wait()

	os.remove(output + "/" + pdf_name)

	#Adapter detection from output of FaQCs
	detection_report = open(output + "/" + prefix + "Subsample_Adapter_Detection.stats.txt")

	detected_adapters = {}
	begin_assessment = False
	for line in detection_report:
		if not begin_assessment:
			if line.strip().startswith("Reads with Adapters/Primers:"):
				begin_assessment = True
		else:
			segment = line.strip().split()
			detected_adapters[segment[0]] = float(re.findall("\d+\.\d+", segment[3])[0])

	detection_report.close()

	clean_detection = []

	for adapter in detected_adapters:
		if detected_adapters[adapter] >= minimum_presence:
			clean_detection.append(adapter)

	#Cleans up after itself.
	for item in seqtk_samples:
		os.remove(item)

	print("Detection done!")

	#Return adapter file
	return clean_detection

#Subsample reads; identify adapters with FaQCs
def adapter_identification_se(artificial_artifacts, seqtk_binary, faqcs_binary, unpaired = "", threads = 1, output = ".", minimum_presence = 0.1, prefix = "", phred_fmt = "33"):
	#seqtk forward and reverse
	subsample = [unpaired, seqtk_binary]

	seqtk_samples = do_seqtk(subsample)

	#FaQCs SE with adapter file
	faqcs_subset_command = [faqcs_binary, "-t", str(threads), "--qc_only", "-d", output, "--artifactFile", artificial_artifacts, "--ascii", phred_fmt]

	#proper naming
	if prefix != "":
		faqcs_subset_command.append("--prefix")
		faqcs_subset_command.append(prefix + "Subsample_Adapter_Detection")
		pdf_name = prefix + "Subsample_Adapter_Detection_qc_report.pdf"
	else :
		faqcs_subset_command.append("--prefix")
		faqcs_subset_command.append("Subsample_Adapter_Detection")
		pdf_name = "Subsample_Adapter_Detection_qc_report.pdf"

	#forward strand
	faqcs_subset_command.append("-u")
	faqcs_subset_command.append(seqtk_samples)


	print("Detecting adapters now... ", end = "")
	ps = subprocess.Popen(faqcs_subset_command)
	ps.wait()

	os.remove(output + "/" + pdf_name)

	#Adapter detection from output of FaQCs
	detection_report = open(output + "/" + prefix + "Subsample_Adapter_Detection.stats.txt")

	detected_adapters = {}
	begin_assessment = False
	for line in detection_report:
		if not begin_assessment:
			if line.strip().startswith("Reads with Adapters/Primers:"):
				begin_assessment = True
		else:
			segment = line.strip().split()
			detected_adapters[segment[0]] = float(re.findall("\d+\.\d+", segment[3])[0])

	detection_report.close()

	clean_detection = []

	for adapter in detected_adapters:
		if detected_adapters[adapter] >= minimum_presence:
			clean_detection.append(adapter)

	#Cleans up after itself.
	os.remove(seqtk_samples)

	print("Detection done!")

	#Return adapter file
	return clean_detection

#gets adapter families for later use
def parse_adapters(full_list, detected_adapters, output, prefix = ""):
	print("Creating specific adapters file for you.")

	#detected adapters is just a list of the user's detected adapters by ID.

	faqcs_internal_adapter_list = faqcs_internal_adapters()

	found = False
	detected_seqs = {}
	for id in detected_adapters:
		found = False
		if id in full_list:
			found = True
			print("Adapter sequence:", id, "detected.")
			detected_seqs[id] = full_list[id]
		if id in faqcs_internal_adapter_list:
			found = True
			print("Adapter sequence:", id, "detected. This adapter is part of a non-optional internal list used by FaQCs and will be included.")
			detected_seqs[id] = faqcs_internal_adapter_list[id]
		#Skip adapter if it cannot be found. Should never happen, now that FaQCs' adapters will always be found and other seqs must be from internal or supplied sequences file
		if not found:
			print("Adapter sequence:", id, "not found in Multitrim's adapter list! It will NOT be included in trimming.")

	adapters_by_family = family_detection(detected_seqs)

	#This is a file I don't want to be temporary. It both helps identify the adapters present in a dataset and provides a fasta for a user to reuse
	subset = open(output + "/" + prefix + "detected_adapters.fasta", "w")

	for adapter in adapters_by_family:
		print(adapter, file = subset)
		print(adapters_by_family[adapter], file = subset)

	subset.close()

	return(output+"/"+ prefix + "detected_adapters.fasta")

#paired end version of the full trim; trims using detected adapters with FaQCs -q 27, then fastp --cut_right window 3 qual 20
def full_trim_pe(forward_in, reverse_in, forward_out, reverse_out, directory, adapters, threads, faqcs, fastp, score, minlen, window, window_qual, prefix, compressor, compress_level, phred_fmt = "33", advanced = False, skip_fastp = False, skip_faqcs = False):
	'''
	Command structure:

	The primary purpose is to issue a FaQCs call on the untrimmed reads, then a subsequent fastp call on the outputs from the FaQCs call.
	Additionally, supports using only one of the two tools. Commands will be built even if the tool is to be skipped, but the call will never be issued.
	'''

	faqcs_command = [faqcs, "-t", str(threads), "-1", forward_in, "-2", reverse_in, "--artifactFile", adapters, "-q", str(score), "--min_L", str(minlen), "--prefix", "reads", "--trim_only", "-d", directory, "--ascii", phred_fmt]
	fastp_command = [fastp, "--thread", str(threads), "--adapter_fasta", adapters, "-l", str(minlen), "--json", directory + "/" + prefix + "post_trim_fastp.json", "--html", directory + "/" + prefix + "post_trim_fastp.html"]

	#Args can be added to fastp command with no consequences if fastp is skipped; command simply won't issue so they will be silent
	if skip_faqcs:
		#This handles taking the input reads directly
		fastp_command.append("-i")
		fastp_command.append(forward_in)
		fastp_command.append("-I")
		fastp_command.append(reverse_in)
	else:
		#FaQCs goes first; this is how I coerce FaQCs reads to look afterwards
		fastp_command.append("-i")
		fastp_command.append(directory+"/reads.1.trimmed.fastq")
		fastp_command.append("-I")
		fastp_command.append(directory+"/reads.2.trimmed.fastq")

	#Outputs are the same regardless of inputs
	fastp_command.append("-o")
	fastp_command.append(forward_out)
	fastp_command.append("-O")
	fastp_command.append(reverse_out)

	if int(window) > 0:
		fastp_command.append("--cut_right")
		fastp_command.append("--cut_right_window_size")
		fastp_command.append(str(window))
		fastp_command.append("--cut_right_mean_quality")
		fastp_command.append(str(window_qual))

	if phred_fmt != "33":
		fastp_command.append("--phred64")

	if advanced:
		fastp_command.append("--trim_poly_g")
		fastp_command.append("--low_complexity_filter")

	time_format = "%d/%m/%Y %H:%M:%S"

	#Manage issuing of commands
	if not skip_faqcs:
		timer = datetime.now()
		printable_time = timer.strftime(time_format)
		print("Trimming with FaQCs. Started at:", printable_time)
		subprocess.run(faqcs_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		os.remove(directory + "/" + "reads.stats.txt")

	if not skip_fastp:
		timer = datetime.now()
		printable_time = timer.strftime(time_format)
		print("Trimming with Fastp. Started at:", printable_time)
		subprocess.run(fastp_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		os.remove(directory + "/" + prefix + "post_trim_fastp.json")
		os.remove(directory + "/" + prefix + "post_trim_fastp.html")



	#We want to rename the non-fastp files, then pass all files and threads to compress with pigz under the nice, neat names

	if skip_fastp:
		#rename FaQCs files to correct names; compress

		#remove this one in any event. We don't want any unpaireds with paired end
		os.remove(directory+"/reads.unpaired.trimmed.fastq")
		shutil.move(directory+"/reads.1.trimmed.fastq", forward_out)
		shutil.move(directory+"/reads.2.trimmed.fastq", reverse_out)
		#compress_commands = [[directory+"/reads.1.trimmed.fastq", forward_out], [directory+"/reads.2.trimmed.fastq", reverse_out]]
		#might as well be parallel
		#pool = multiprocessing.Pool(min(2, threads))
		#pool.map(compress_faqcs, compress_commands)
		#pool.close()

	elif not skip_faqcs:
		#remove FaQCs files if fastp has results or skip if FaQCs not done.
		os.remove(directory+"/reads.1.trimmed.fastq")
		os.remove(directory+"/reads.2.trimmed.fastq")
		#remove this one in any event. We don't want any unpaireds with paired end - the call has to be duplicated, unfortunately.
		os.remove(directory+"/reads.unpaired.trimmed.fastq")

	compress_results([forward_out, reverse_out], threads, compressor, compress_level)

	return None

#single end version of the full trim; trims using detected adapters with FaQCs -q 27, then fastp --cut_right window 3 qual 20
def full_trim_se(reads_in, reads_out, directory, adapters, threads, faqcs, fastp, score, minlen, window, window_qual, prefix, compressor, compress_level, phred_fmt = "33", advanced = False, skip_fastp = False, skip_faqcs = False):
	'''
	Command structure:

	The primary purpose is to issue a FaQCs call on the untrimmed reads, then a subsequent fastp call on the outputs from the FaQCs call.
	Additionally, supports using only one of the two tools. Commands will be built even if the tool is to be skipped, but the call will never be issued.
	'''
	faqcs_command = [faqcs, "-t", str(threads), "-u", reads_in, "--artifactFile", adapters, "-q", str(score), "--min_L", str(minlen), "--prefix", "reads", "--trim_only", "-d", directory, "--ascii", phred_fmt]
	fastp_command = [fastp, "--thread", str(threads), "--adapter_fasta", adapters, "-l", str(minlen), "--json", directory + "/" + prefix + "post_trim_fastp.json", "--html", directory + "/" + prefix + "post_trim_fastp.html"]

	#Args can be added to fastp command with no consequences if fastp is skipped; command simply won't issue so they will be silent
	if skip_faqcs:
		#This handles taking the input reads directly
		fastp_command.append("-i")
		fastp_command.append(reads_in)
	else:
		#FaQCs goes first; this is how I coerce FaQCs reads to look afterwards
		fastp_command.append("-i")
		fastp_command.append(directory+"/reads.unpaired.trimmed.fastq")

	#Outputs are the same regardless of inputs
	fastp_command.append("-o")
	fastp_command.append(reads_out)

	if int(window) > 0:
		fastp_command.append("--cut_right")
		fastp_command.append("--cut_right_window_size")
		fastp_command.append(str(window))
		fastp_command.append("--cut_right_mean_quality")
		fastp_command.append(str(window_qual))

	if phred_fmt != "33":
		fastp_command.append("--phred64")

	if advanced:
		fastp_command.append("--trim_poly_g")
		fastp_command.append("--low_complexity_filter")

	time_format = "%d/%m/%Y %H:%M:%S"

	#Manage issuing of commands
	if not skip_faqcs:
		timer = datetime.now()
		printable_time = timer.strftime(time_format)
		print("Trimming with FaQCs. Started at:", printable_time)
		subprocess.run(faqcs_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		os.remove(directory + "/" + "reads.stats.txt")

	if not skip_fastp:
		timer = datetime.now()
		printable_time = timer.strftime(time_format)
		print("Trimming with Fastp. Started at:", printable_time)
		subprocess.run(fastp_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		os.remove(directory + "/" + prefix + "post_trim_fastp.json")
		os.remove(directory + "/" + prefix + "post_trim_fastp.html")


	if skip_fastp:
		#compress the result
		#remove this one in any event. We don't want any unpaireds with paired end
		shutil.move(directory+"/reads.unpaired.trimmed.fastq", reads_out)
	elif not skip_faqcs:
		#remove FaQCs files if fastp has results or skip if FaQCs not run.
		os.remove(directory+"/reads.unpaired.trimmed.fastq")

	compress_results([reads_out], threads, compressor, compress_level)

	return None

#compress results using selected compressor.
def compress_results(output_files, threads, compressor, level):

	#Just for printing feedback.
	pretty_compressor = ["GZIP", "PIGZ", "DSRC-2"][["gzip", "pigz", "dsrc"].index(compressor)]

	time_format = "%d/%m/%Y %H:%M:%S"
	timer = datetime.now()
	printable_time = timer.strftime(time_format)
	print("Beginning compression of trimmed reads using", pretty_compressor, "at:", printable_time)

	#These get the file sizes, runtimes, compression ratio
	if compressor == "gzip":
		gzip_compress_module(output_files, threads, level)

	if compressor == "pigz":
		pigz_compress_module(output_files, threads, level)

	'''
	if compressor == "dsrc":
		dsrc_compress_module(output_files, threads, level)

	'''

	print("Outputs compressed!")

	return None

#gzip is NOT threaded, so we open up to 4 threads and compress each input simultaneously, feeding the results to falco as we go.
def gzip_compress_module(outputs, threads, level):
	#Get the gzip set up for each file
	gzip_arguments = []
	for file in outputs:
		gzip_arguments.append(["gzip", "-f"+str(level), file])

	#Don't open more threads than you have to.
	num_files = len(outputs)
	#Run args
	pool = multiprocessing.Pool(min(threads, num_files))
	pool.map(do_gzip_pretty, gzip_arguments)
	pool.close()

	return None

#The particular parallelization for this is a bother.
def do_gzip_pretty(compress_argument):
	file = compress_argument[2]
	start_time = datetime.now()
	initial_size = os.path.getsize(file)

	subprocess.call(compress_argument)

	final_size = os.path.getsize(file+".gz")
	end_time = datetime.now()

	pretty_print_file_size(file, initial_size, final_size, start_time, end_time)

	return None

#pigz is threaded, so compression happens 1 file at a time using all threads, then falco QC 4 using the gzip approach above since the result is in gzip format
def pigz_compress_module(outputs, threads, level):
	for file in outputs:
		start_time = datetime.now()
		initial_size = os.path.getsize(file)

		pigz_argument = ["pigz", "-f", "-"+str(level), "-p", str(threads), file]
		subprocess.call(pigz_argument)

		final_size = os.path.getsize(file+".gz")
		end_time = datetime.now()

		pretty_print_file_size(file, initial_size, final_size, start_time, end_time)

	return None

#Unfinished, Has more moving parts to take care of.
#DSRC-2 is threaded, but the compressed format is not supported by falco. Thus, we run QC, THEN compress each file 1 at a time using all threads.
def dsrc_compress_module(inputs, outputs, threads, level):
	print("DSRC-2 will also produce QC reports at this time!")

	#DSRC only accepts up to 64 threads
	if threads > 64:
		threads = 64

	#DSRC-formatted args
	threads = "-t"+str(threads)
	level   = "-m"+str(level)

	#falco goes here for DSRC-2, must be uncompressed files.
	num_files = min(threads, len(inputs)+len(outputs))


	for file in files:
		output_file_name = file+".dsrc"

		start_time = datetime.now()
		initial_size = os.path.getsize(files[i])

		compress_command = ["dsrc", "c", threads, level, file, output_file_name]
		subprocess.run(compress_command)

		ending_size = os.path.getsize(output_file_name)
		end_time = datetime.now()

		pretty_print_file_size(files[i], initial_size, ending_size, start_time, end_time)

	print("Compression and QC complete!")

#Unfinished.
#Function for checking if an input file is a DSRC archive - these have to be decompressed for trimming, since the tools don't directly support such archives.
def check_is_dsrc(file):
	#We're going to make a file in a temporary directory and use it
	base_name = os.path.basename(file)
	loc = tempfile.mkdtemp()
	tempout = loc + "/" + base_name
	is_dsrc = False

	#Attempt to DSRC decompress into the temp file
	try:
		#Multiple reasons this could fail, including tool absence. All should be handled by this except.
		dsrc_decomp = ["dsrc", "d", file, tempout]
		subprocess.run(dsrc_decomp, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
		#DSRC only creates the file if it's successful in opening the file and DSRC can be called in the first place.
		is_dsrc = os.path.exists(tempout)
	#If the file cannot be decompressed, delete the temp file and return self.
	except:
		shutil.rmtree(loc)

	if is_dsrc:
		dsrc_file = tempout
	else:
		dsrc_file = file

	return is_dsrc, dsrc_file

#Convert a file's size in bytes to human-readable format.
def humansize(nbytes):
	suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
	i = 0
	while nbytes >= 1024 and i < len(suffixes)-1:
		nbytes /= 1024.
		i += 1
	f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
	return '%s %s' % (f, suffixes[i])

#Print well-formatted compression time info.
def pretty_print_file_size(name, start, end, start_time, end_time):
	runtime = end_time - start_time

	try:
		hours = runtime.hours
	except:
		hours = 0

	try:
		minutes = runtime.minutes
	except:
		minutes = 0

	try:
		seconds = runtime.seconds
	except:
		seconds = 0

	runtime = '%02d:%02d:%02d' % (hours, minutes, seconds)

	print(name, "compressed! Compression took:", runtime, "and the file was compressed to",  str(round((end/start)*100, 2)), "percent of original size from", humansize(start), "to", humansize(end))

	return None

#Stolen from a SO thread on how to issue usage information on an error.
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

#Option parsing
def gather_opts():
	parser = MyParser(description=''' This program is designed to facilitate effective trimming of your reads.
	It will help to identify the presence of adapters in your reads, trim those adapters and the reads efficiently,
	and produce several bfore and after quality reports in addition to the trimmed reads. This is a pipeline incorporating
	FaQCs, falco, and seqtk commands, in addition to several python operations which exist to facilitate adapter finding and
	subsetting. --user and --UNLIMITED_POWER are jokes, but you should usually use --UNLIMITED_POWER.''')
	#Use all available cores.
	parser.add_argument("--max", dest = "Sheev", action = 'store_true', help = "Attempts to detect and use all available processors for threading.")
	#Or this many threads. Laaaaame
	parser.add_argument("--threads", "-t", dest = "threads", default = 1, help = "Number of threads to use for parallel processes. Default 1")

	#file inputs
	parser.add_argument("--forward", "-1", dest = "f", default = "", help = "Forward Strand Reads (use -u for unpaired reads)")
	parser.add_argument("--reverse", "-2", dest = "r", default = "", help = "Reverse Strand Reads (use -u for unpaired reads)")
	parser.add_argument("--unpaired", "-u", dest = "u", default = "", help = "Unpaired Reads")

	#final out directory
	parser.add_argument("--output", "-o", dest = "outdir", default = ".", help = "Directory to send final outputs.")
	#naming convention
	parser.add_argument("--prefix", "-p", dest = "pref", default = "", help = "Prefix to place on outputs.")

	#Adapter detection opts
	parser.add_argument("--min_adapt_pres", "-m", dest = "minpres", default = 0.1, help = "Minimum presence of an adapter for it to be considered present in a set of reads. Default 0.1, so an adapter is considered present if detected in 0.1 percent of reads.")

	parser.add_argument("--adapters", "-a", dest = "adapter_fasta", default = "internal", help = "Supply a custom set of adapters for adapter detection. Detected adapters can come only from this set. Multitrim uses the MiGA adapter set by defualt.")
	#parser.add_argument("--kits", dest = "adapter_families", default = "internal", help = "Supply a 2-column, comma separated list of adapter IDs and kits of origin. When an adapter is detected, all adapters in the same seq. prep kit are also considered detected when using the default MiGA adapters.")

	#Shared options
	parser.add_argument("--min_L", "-l", dest = "length", default = "50", help = "Minimum read length. Default 50 base pairs.")
	parser.add_argument("--phred_fmt", dest = "phred", default = "33", help = "Phred q score format (default 33)")
	parser.add_argument("--advanced", dest = "advanced", action = 'store_true', help = "Apply advanced trimming options (poly-G tail, low-complexity). Only useful for reads sequenced with 2-dye chemistry.")

	#FaQCs opts
	parser.add_argument("--score", "-s", dest= "score", default = "27", help = "FaQCs quality target. Default 27")
	parser.add_argument("--skip_faqcs", dest = "skip_fq", action = 'store_true', help = "Do not trim with FaQCs (use fastp only). Cannot skip both.")

	#fastp opts
	parser.add_argument("--window", "-w", dest = "mid", default = "3", help = "Trimmomatic-like sliding window. Default 3.")
	parser.add_argument("--window_qual", "-q", dest = "mid_q", default = "20", help = "Trim quality cutoff for trimmomatic window. Default 20.")
	parser.add_argument("--skip_fastp", dest = "skip_fp", action = 'store_true', help = "Do not trim with fastp (use FaQCs only). Cannot skip both.")

	#parser.add_argument("--falco", dest = "falco_path", default = "falco", help = "Location of Falco QC binary.")
	#parser.add_argument("--seqtk", dest = "seqtk_path", default = "seqtk", help = "Location of SeqTK binary.")
	#parser.add_argument("--faqcs", dest = "faqcs_path", default = "FaQCs", help = "Location of FaQCs binary.")
	#parser.add_argument("--fastp", dest = "fastp_path", default = "fastp", help = "Location of fastp binary.")

	#Update with DSRC later
	#parser.add_argument("--zip", dest = "compressor", default = "gzip", help = "Select a compressor for outputs. Supported options are: 'gzip' (default), 'pigz', 'dsrc'")
	#-1 default value used for automating pigz/gzip vs dsrc selection
	#parser.add_argument("--level", dest = "zip_level", default = -1, help = "Choose a compression level for outputs. gzip, pigz take values 1-9 with default 6. DSRC takes 0-2 with default 0. Higher values are slower but compress better.")

	parser.add_argument("--zip", dest = "compressor", default = "gzip", help = "Select a compressor for outputs. Supported options are: 'gzip' (default), 'pigz'")
	#-1 default value used for automating pigz/gzip vs dsrc selection
	parser.add_argument("--level", dest = "zip_level", default = -1, help = "Choose a compression level for outputs. gzip and pigz take values 1-9 with default 6")

	parser.add_argument("--resources", dest = "resource_list", action = 'store_true', help = "Print a list of resources used by Multitrim and quit.")


	return(parser, parser.parse_args())

def print_resources():
	print("Multitrim github: https://github.com/KGerhardt/multitrim")
	print("MiGA adapters available at: https://github.com/bio-miga/miga/blob/main/utils/adapters.fa")
	internal_adapters = faqcs_internal_adapters()
	print("FaQCs mandatory adapters are:")
	for id in internal_adapters:
		print(id, internal_adapters[id])
	print("FaQCs github: https://github.com/LANL-Bioinformatics/FaQCs")
	print("fastp github: https://github.com/OpenGene/fastp")
	print("Falco github: https://github.com/smithlabcode/falco")

#Program Control
def main():
	#Keep the parser on hand so I can prent usage as needed.s
	help_message, options = gather_opts()

	resources = options.resource_list
	if resources:
		print_resources()
		quit()


	#Allows for the script to take no inputs and print help/usage
	if len(sys.argv)==1:
		help_message.print_help(sys.stderr)
		quit()


	skip_fq = options.skip_fq
	skip_fp = options.skip_fp

	if skip_fp and skip_fq:
		print("Cannot skip both trimming tools. This would result in no trim at all. Exiting program.")
		sys.exit(1)

	#file name prefix
	prefix = str(options.pref)

	#Make it more convenient for me later
	if prefix != "":
		if not prefix.endswith("_"):
			prefix = prefix + "_"

	#Tool names
	fp = "fastp"
	fq = "FaQCs"
	stk = "seqtk"
	that_aint_falco = "falco"

	#Get the reads
	f = options.f
	r = options.r
	u = options.u

	#phred format
	phred = str(options.phred)

	#num threads
	threads = int(options.threads)
	#Check for --max flag
	if options.Sheev:
		#Detects and uses all the threads a system has available.
		try:
			threads = len(os.sched_getaffinity(0))
		except:
			print("Cannot detect how many cores are available! Defaulting to 1. Use --threads to specify more cores if you see this message.")
			threads = 1
	else:
		#Check to ensure a user doesn't request more procs than available.
		try:
			threads = min(threads, len(os.sched_getaffinity(0)))
		except:
			#Handle case where the sched getaffinity function is unavailable.
			threads = int(options.threads)



	#No reads shorter than minlen
	minlen = str(options.length)

	#advanced trimming opts
	#FaQCs:
		# currently no advanced opts
	#Fastp:
		# --trim_poly_g
		# --low_complexity_filter
	advanced = options.advanced

	#These options control the trimming behavior for fastp, correspond to sliding window width and avg. quality min.
	mid = str(options.mid)
	mid_q = str(options.mid_q)

	#faqcs target score. Lower = less aggressive, higher = more aggressive
	score = str(options.score)

	#directory to place results. Creates if needed, but won't create multiple dirs.
	final_output = options.outdir

	#Autocomplete may include the slash, but I don't want it
	if final_output.endswith("/"):
		final_output = final_output[:-1]

	#Check to make sure it actually has data, or exits
	if f == "" and r == "" and u == "":
		print("I need to be given reads! Use -1 and -2 for paired-end reads, or -u for unpaired reads. Exiting program.")
		quit()

	#Check to make sure that a forward read is paired with a reverse read if either is supplied
	if f == "" and r != "" or f != "" and r == "":
		print("If you have paired reads, I need both the forward and reverse files. If you just want to process one, use -u to specify it. Exiting program.")
		quit()

	#fastp cannot take both unpaired and paired simultaneously
	if u != "" and (r != "" or f != ""):
		print("If you have paired reads, I need both the forward and reverse files. If you just want to process one, use -u to specify it. Exiting program.")
		quit()

	#Determine single or paired end mode
	if u == "":
		paired_end = True
	else:
		paired_end = False

	if paired_end:
		quit_out = False
		if not os.path.exists(f):
			print("Forward Reads: " + f + " not found. Multitrim will exit.")
			quit_out = True
		if not os.path.exists(r):
			print("Reverse Reads: " + r + " not found. Multitrim will exit.")
			quit_out = True

		if quit_out:
			quit()

	else:
		if not os.path.exists(u):
			print("Unpaired Reads: " + u + " not found. Multitrim will exit.")
			quit()


	#Check if a directory is specified and which doesn't exist; try to create if needed or exit gracefully.
	#This has to happen last, or it risks making the directory when the program is otherwise going to quit.
	if final_output != ".":
		if not os.path.exists(final_output):
			try:
				os.mkdir(final_output)
			except:
				print("Multitrim wasn't able to find or create the specified output directory. Exiting program.")
				quit()


	compressor = options.compressor
	compression_level = int(options.zip_level)

	#if compressor not in ['gzip', 'pigz', 'dsrc']:
	if compressor not in ['gzip', 'pigz']:
		#print("Chosen compressor '"+compressor+"' not supported! Supported options are: 'gzip', 'pigz', 'dsrc'")
		print("Chosen compressor '"+compressor+"' not supported! Supported options are: 'gzip', 'pigz'")
		quit()

	if compressor in ['gzip', 'pigz']:
		if compression_level == -1:
			compression_level = 6

		if not 1 <= compression_level <= 9:
			print("Compression level", compression_level, "not acceptable! For GZIP and PIGZ, supported compression levels are 1-9!")
			quit()

	#For DSRC development later
	'''
	if compressor in ['dsrc']:
		if compression_level == -1:
			compression_level = 0

		if not 0 <= compression_level <= 3:
			print("Compression level", compression_level, "not acceptable! For DSRC-2, supported compression levels are 0-2!")
			quit()
	'''

	#Check for input adapter file or default to internal list.
	input_adapters = options.adapter_fasta
	if not os.path.exists(input_adapters) and input_adapters != "internal":
		print("Adapters file", input_adapters, "not found! Exiting.")
		quit()
	#Adapters detected if minpres% of reads have that specific adapter present according to FaQCs stats.
	minpres = float(options.minpres)

	#Reads user adapters or creates all adapters file from scratch
	adapter_set, complete_adapter_file_name, needs_cleanup = read_adapters(input_adapters)

	#User feedback
	if final_output == ".":
		print("Placing results in:", os.getcwd())
	else:
		print("Placing results in:", final_output)


	if options.Sheev:
		print("Using all available cores. Number of cores:", threads)
	else:
		print("Working with", threads, "threads.")

	if paired_end:
		#two inputs; paired end behavior
		print("Primary Strand Reads:", f, "\nReverse Strand Reads:", r)
		#User feedback
		print("Adapters considered detected if present in "+ str(minpres) + " % of reads.")

		pre_qc_f, pre_qc_r, post_qc_f, post_qc_r, post_trim_f, post_trim_r = names_pe(f, r, final_output, prefix)

		adapters_detected = adapter_identification_pe(complete_adapter_file_name, stk, fq, f, r, threads, final_output, minpres, prefix, phred)
		cleaned_adapters = parse_adapters(adapter_set, adapters_detected, final_output, prefix)

		if needs_cleanup:
			print("Removing automatically generated adapters...")
			os.remove(complete_adapter_file_name)

		full_trim_pe(f, r, post_trim_f, post_trim_r, final_output, cleaned_adapters, threads, fq, fp, score, minlen, mid, mid_q, prefix, compressor, compression_level, phred, advanced, skip_fp, skip_fq)

		if compressor not in ["dsrc"]:
			falco_qc_pe(f, r, post_trim_f, post_trim_r, pre_qc_f, post_qc_f, pre_qc_r, post_qc_r, threads, that_aint_falco)

	else:
		#one input; SE behavior
		print("Unpaired Reads:", u)
		#User feedback
		print("Adapters considered detected if present in "+ str(minpres) + " % of reads.")

		pre_qc, post_qc, post_trim = names_se(u, final_output, prefix)

		adapters_detected = adapter_identification_se(complete_adapter_file_name, stk, fq, u, threads, final_output, minpres, prefix, phred)
		cleaned_adapters = parse_adapters(adapter_set, adapters_detected, final_output, prefix)

		if needs_cleanup:
			print("Removing automatically generated adapters...")
			os.remove(complete_adapter_file_name)

		full_trim_se(u, post_trim, final_output, cleaned_adapters, threads, fq, fp, score, minlen, mid, mid_q, prefix, compressor, compression_level, phred, advanced, skip_fp, skip_fq)

		if compressor not in ["dsrc"]:
			falco_qc_se(u, post_trim, pre_qc, post_qc, threads, that_aint_falco)


	print("Trimming complete.")

#just runs main
if __name__ == "__main__":
	main()

#End of functional components of Multitrim.

#Below are leftover creation functions that could be used to update the list of adapters. Not used in the program proper.

#Regenerate adatapers file output from a fasta. This is a utility function I do not expect to see used in the final product.
def fasta_to_permanent_python(original_adapters_fasta):
	fasta = open(original_adapters_fasta, "r")

	fasta_seq_dict = {}

	current_line = fasta.readline().strip()

	current_id = current_line
	current_seq = ""

	current_line = fasta.readline().strip()

	while current_line:
		if current_line.startswith(">"):
			fasta_seq_dict[current_id] = current_seq
			current_id = current_line
			current_seq = ""
		else:
			current_seq += current_line

		current_line = fasta.readline().strip()

	#Finally, python needs this logic.
	fasta_seq_dict[current_id] = current_seq

	fasta.close()

	for contig in fasta_seq_dict:
		print("adapters_dict[\""+contig+"\"] = \""+ fasta_seq_dict[contig]+"\"")

#These spit out spoofed python code for the in-built creation of an adapters file to supply tools, without the external need for this file.
#I just copy-paste the results to tbe generate_adapters_temporary_file function's body to get the results.

#As above, this prints out a python-correct set of commands for me to copy-paste.
#This one produces a set of "families" for adapters, where each is the set of adapters in a kit.
def fasta_to_families(original_adapters_fasta):
	fasta = open(original_adapters_fasta, "r")

	whichfam = [5, 6, 10, 6, 7, 12, 14, 3, 13, 50, 9, 5, 6]

	families = []
	families.append("singleend")
	families.append("pairedend")
	families.append("dpnII")
	families.append("smallrna")
	families.append("multiplex")
	families.append("pcr")
	families.append("dpnIIgex")
	families.append("otherrna")
	families.append("trueseq")
	families.append("rnapcr")
	families.append("abi")
	families.append("trueseq2")
	families.append("nextera")

	current_fam = 0

	famlist = []

	for family_size in whichfam:
		for i in range(0, family_size):
			famlist.append(families[current_fam])
		current_fam += 1

	fasta_fam_dict = {}

	current_fam = 0

	for line in fasta:
		if line.strip().startswith(">"):
			fasta_fam_dict[line.strip()[1:]] = famlist[current_fam]
			current_fam += 1

	fasta.close()

	for contig in fasta_fam_dict:
		print("adapters_fam_dict[\""+contig+"\"] = \""+ fasta_fam_dict[contig]+"\"")

#identifies the same adapters as in the full file with a family of origin, so that all adapters in a family can be selected.
def create_adapter_families():
	adapters_fam_dict = {}

	adapters_fam_dict["Illumina_Single_End_Apapter_1"] = "singleend"
	adapters_fam_dict["Illumina_Single_End_Apapter_2"] = "singleend"
	adapters_fam_dict["Illumina_Single_End_PCR_Primer_1"] = "singleend"
	adapters_fam_dict["Illumina_Single_End_PCR_Primer_2"] = "singleend"
	adapters_fam_dict["Illumina_Single_End_Sequencing_Primer"] = "singleend"
	adapters_fam_dict["Illumina_Paired_End_Adapter_1"] = "pairedend"
	adapters_fam_dict["Illumina_Paired_End_Adapter_2"] = "pairedend"
	adapters_fam_dict["Illumina_Paried_End_PCR_Primer_1"] = "pairedend"
	adapters_fam_dict["Illumina_Paired_End_PCR_Primer_2"] = "pairedend"
	adapters_fam_dict["Illumina_Paried_End_Sequencing_Primer_1"] = "pairedend"
	adapters_fam_dict["Illumina_Paired_End_Sequencing_Primer_2"] = "pairedend"
	adapters_fam_dict["Illumina_DpnII_expression_Adapter_1"] = "dpnII"
	adapters_fam_dict["Illumina_DpnII_expression_Adapter_2"] = "dpnII"
	adapters_fam_dict["Illumina_DpnII_expression_PCR_Primer_1"] = "dpnII"
	adapters_fam_dict["Illumina_DpnII_expression_PCR_Primer_2"] = "dpnII"
	adapters_fam_dict["Illumina_DpnII_expression_Sequencing_Primer"] = "dpnII"
	adapters_fam_dict["Illumina_NlaIII_expression_Adapter_1"] = "dpnII"
	adapters_fam_dict["Illumina_NlaIII_expression_Adapter_2"] = "dpnII"
	adapters_fam_dict["Illumina_NlaIII_expression_PCR_Primer_1"] = "dpnII"
	adapters_fam_dict["Illumina_NlaIII_expression_PCR_Primer_2"] = "dpnII"
	adapters_fam_dict["Illumina_NlaIII_expression_Sequencing_Primer"] = "dpnII"
	adapters_fam_dict["Illumina_Small_RNA_Adapter_1"] = "smallrna"
	adapters_fam_dict["Illumina_Small_RNA_Adapter_2"] = "smallrna"
	adapters_fam_dict["Illumina_Small_RNA_RT_Primer"] = "smallrna"
	adapters_fam_dict["Illumina_Small_RNA_PCR_Primer_1"] = "smallrna"
	adapters_fam_dict["Illumina_Small_RNA_PCR_Primer_2"] = "smallrna"
	adapters_fam_dict["Illumina_Small_RNA_Sequencing_Primer"] = "smallrna"
	adapters_fam_dict["Illumina_Multiplexing_Adapter_1"] = "multiplex"
	adapters_fam_dict["Illumina_Multiplexing_Adapter_2"] = "multiplex"
	adapters_fam_dict["Illumina_Multiplexing_PCR_Primer_1.01"] = "multiplex"
	adapters_fam_dict["Illumina_Multiplexing_PCR_Primer_2.01"] = "multiplex"
	adapters_fam_dict["Illumina_Multiplexing_Read1_Sequencing_Primer"] = "multiplex"
	adapters_fam_dict["Illumina_Multiplexing_Index_Sequencing_Primer"] = "multiplex"
	adapters_fam_dict["Illumina_Multiplexing_Read2_Sequencing_Primer"] = "multiplex"
	adapters_fam_dict["Illumina_PCR_Primer_Index_1"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_2"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_3"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_4"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_5"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_6"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_7"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_8"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_9"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_10"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_11"] = "pcr"
	adapters_fam_dict["Illumina_PCR_Primer_Index_12"] = "pcr"
	adapters_fam_dict["Illumina_DpnII_Gex_Adapter_1"] = "dpnIIgex"
	adapters_fam_dict["Illumina_DpnII_Gex_Adapter_1.01"] = "dpnIIgex"
	adapters_fam_dict["Illumina_DpnII_Gex_Adapter_2"] = "dpnIIgex"
	adapters_fam_dict["Illumina_DpnII_Gex_Adapter_2.01"] = "dpnIIgex"
	adapters_fam_dict["Illumina_DpnII_Gex_PCR_Primer_1"] = "dpnIIgex"
	adapters_fam_dict["Illumina_DpnII_Gex_PCR_Primer_2"] = "dpnIIgex"
	adapters_fam_dict["Illumina_DpnII_Gex_Sequencing_Primer"] = "dpnIIgex"
	adapters_fam_dict["Illumina_NlaIII_Gex_Adapter_1.01"] = "dpnIIgex"
	adapters_fam_dict["Illumina_NlaIII_Gex_Adapter_1.02"] = "dpnIIgex"
	adapters_fam_dict["Illumina_NlaIII_Gex_Adapter_2.01"] = "dpnIIgex"
	adapters_fam_dict["Illumina_NlaIII_Gex_Adapter_2.02"] = "dpnIIgex"
	adapters_fam_dict["Illumina_NlaIII_Gex_PCR_Primer_1"] = "dpnIIgex"
	adapters_fam_dict["Illumina_NlaIII_Gex_PCR_Primer_2"] = "dpnIIgex"
	adapters_fam_dict["Illumina_NlaIII_Gex_Sequencing_Primer"] = "dpnIIgex"
	adapters_fam_dict["Illumina_5p_RNA_Adapter"] = "otherrna"
	adapters_fam_dict["Illumina_RNA_Adapter1"] = "otherrna"
	adapters_fam_dict["Illumina_Small_RNA_3p_Adapter_1"] = "otherrna"
	adapters_fam_dict["TruSeq_Universal_Adapter"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_1"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_2"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_3"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_4"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_5"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_6"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_7"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_8"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_9"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_10"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_11"] = "trueseq"
	adapters_fam_dict["TruSeq_Adapter_Index_12"] = "trueseq"
	adapters_fam_dict["Illumina_RNA_RT_Primer"] = "rnapcr"
	adapters_fam_dict["Illumina_RNA_PCR_Primer"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_1"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_2"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_3"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_4"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_5"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_6"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_7"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_8"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_9"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_10"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_11"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_12"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_13"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_14"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_15"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_16"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_17"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_18"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_19"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_20"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_21"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_22"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_23"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_24"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_25"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_26"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_27"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_28"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_29"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_30"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_31"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_32"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_33"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_34"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_35"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_36"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_37"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_38"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_39"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_40"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_41"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_42"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_43"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_44"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_45"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_46"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_47"] = "rnapcr"
	adapters_fam_dict["RNA_PCR_Primer_Index_48"] = "rnapcr"
	adapters_fam_dict["ABI_Dynabead_EcoP_Oligo"] = "abi"
	adapters_fam_dict["ABI_Solid3_Adapter_A"] = "abi"
	adapters_fam_dict["ABI_Solid3_Adapter_B"] = "abi"
	adapters_fam_dict["ABI_Solid3_5_AMP_Primer"] = "abi"
	adapters_fam_dict["ABI_Solid3_3_AMP_Primer"] = "abi"
	adapters_fam_dict["ABI_Solid3_EF1_alpha_Sense_Primer"] = "abi"
	adapters_fam_dict["ABI_Solid3_EF1_alpha_Antisense_Primer"] = "abi"
	adapters_fam_dict["ABI_Solid3_GAPDH_Forward_Primer"] = "abi"
	adapters_fam_dict["ABI_Solid3_GAPDH_Reverse_Primer"] = "abi"
	adapters_fam_dict["TruSeq2_SE"] = "trueseq2"
	adapters_fam_dict["TruSeq2_PE_f"] = "trueseq2"
	adapters_fam_dict["TruSeq2_PE_r"] = "trueseq2"
	adapters_fam_dict["TruSeq3_IndexedAdapter"] = "trueseq2"
	adapters_fam_dict["TruSeq3_UniversalAdapter"] = "trueseq2"
	adapters_fam_dict["Nextera_PE_PrefixNX/1"] = "nextera"
	adapters_fam_dict["Nextera_PE_PrefixNX/2"] = "nextera"
	adapters_fam_dict["Nextera_PE_Trans1"] = "nextera"
	adapters_fam_dict["Nextera_PE_Trans1_rc"] = "nextera"
	adapters_fam_dict["Nextera_PE_Trans2"] = "nextera"
	adapters_fam_dict["Nextera_PE_Trans2_rc"] = "nextera"

	return(adapters_fam_dict)

#creates python code for family:id:sequence dictionary for kit detection code. Relies on the current state of the generate_adapters_temporary_file function.
def create_seq_to_fam():
	dict, name = generate_adapters_temporary_file()

	fams = create_adapter_families()

	seq_to_fam = {}

	for id in dict:
		family = fams[id[1:]]
		seq_to_fam[dict[id]] = family
		#print("seq_to_fam['"+dict[id]+"'] = '" + family+ "'")

	fam_to_id_to_seq = {}

	for id in dict:
		family = fams[id[1:]]
		seq = dict[id]

		if family in fam_to_id_to_seq:
			fam_to_id_to_seq[family][id[1:]] = seq
		else:
			fam_to_id_to_seq[family] = {}
			fam_to_id_to_seq[family][id[1:]] = seq

	for fam in fam_to_id_to_seq:
		print("fam_to_id_to_seq['"+fam+"'] =", fam_to_id_to_seq[fam])








	os.remove(name)
