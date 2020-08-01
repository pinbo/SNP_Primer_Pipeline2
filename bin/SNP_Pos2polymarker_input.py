#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  SNP_Pos2polymarker_input.py
#  
#  Copyright 2017 Junli Zhang <zhjl86@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
# example: input file should have NO header
# SNP_Pos2polymarker_input.py <input> <output> <reference number>
# SNP_Pos2polymarker_input.py snp_infor.txt for_polymarker.csv reference
##############################################
# Modules needed
from itertools import groupby
import sys
from subprocess import call

##############################################
# SNP class to store all the related information for a SNP
class SNP:
	''' Object representing a SNP record. '''
	def __init__(self, contig, ref_pos, ref_allele, alt_allele):
		self.chr = contig
		self.name = contig + "-" + str(ref_pos)
		xstream = 50 # get 50 bps on each side
		self.contig = contig
		self.ref_pos = ref_pos
		self.ref_allele = ref_allele
		self.alt_allele = alt_allele
		self.seq = ""
		self.leftpos = max(1, ref_pos - xstream)
		self.rightpos = ref_pos + xstream

# function to parse the snp information file from the exon capture
# example (NO header please): Contig  Ref_Pos      Ref     Alternative
#IWGSC_CSS_7AL_scaff_4491815     4083       T       A
#IWGSC_CSS_7AL_scaff_4552312     3240       G       T
#chr7A	1111400000	G	C
def parse_exon_snp(snpinfo):
	snpdict = {} # dictionary of snp information
	seq_name_list = [] # for changing the sequence names from the blastdbcmd output
	with open(snpinfo) as infile:
		#next(infile) # skip header
		for line in infile:
			line = line.strip()
			if not line:
				continue
			contig, ref_pos, ref_allele, alt_allele = line.rstrip().split() # split with white space, tab or space are all okay
			snpdict[contig + "-" + ref_pos] = SNP(contig, int(ref_pos), ref_allele, alt_allele)
			seq_name_list.append(contig + "-" + ref_pos)
	return snpdict, seq_name_list

# function to prepare file for blastdbcmd to get the flanking sequences of SNPs
def prepare_seq_range(snpdict, seq_name_list, outfile): # seq_name_list to make sure the output order is the same as input order
	out = open(outfile, "w")
	for k in seq_name_list:
		v = snpdict[k]
		contig = "-".join(k.split("-")[:-1]) # k = contig + "-" + ref_pos, in case there are "-" in the chromosome/contig names
		out.write(contig + "\t" + str(v.leftpos) + "-" + str(v.rightpos) + "\n")
	out.close()
	return 0

# function to get the flanking sequences
def get_flanking(range_file, flanking_file, reference):
	cmd = "blastdbcmd -entry_batch " + range_file + " -db " + reference + " > " + flanking_file
	print("Command to get the flanking sequences for each SNP\n", cmd)
	call(cmd, shell=True)

# function to extract sequences from a fasta file 
def get_fasta(infile, seq_name_list):
	fasta = {} # dictionary for alignment
	n = 0 # for sequence name list index
	with open(infile) as file_one:
		for line in file_one:
			line = line.strip()
			if not line:
				continue
			if line.startswith(">"):
				sequence_name = line.split()[0].split(":")[0].lstrip(">") # so both blast+2.9 or earlier work
				if sequence_name in seq_name_list[n]: # in case seq name mismatch
					sequence_name = seq_name_list[n]
					n += 1
			else:
				fasta.setdefault(sequence_name, "")
				fasta[sequence_name] += line.rstrip()
	return fasta

# main function
def main(argv):
	#print("\n\n!!!Input file should have NO header!!!\nEach line is a SNP!\n\n")
	snpinfo = argv[1] # input file
	outfile = argv[2] # output file
	reference = argv[3]
	snpdict, seq_name_list = parse_exon_snp(snpinfo)
	print("length of snpdict ", len(snpdict))
	range_file = "temp_range.txt"
	prepare_seq_range(snpdict, seq_name_list, range_file)
	flanking_file = "flanking_seq.fa"
	get_flanking(range_file, flanking_file, reference)
	seq_fasta = get_fasta(flanking_file, seq_name_list)
	out = open(outfile, "w")
	#for header, seq in fasta:
	for i in seq_name_list:
		snp = snpdict[i]
		seq = seq_fasta[i]
		if (snp.leftpos == 1):
			snp.seq = seq[0:-52] + "[" + snp.ref_allele + "/" + snp.alt_allele + "]" + seq[-50:]
		else:
			snp.seq = seq[0:50] + "[" + snp.ref_allele + "/" + snp.alt_allele + "]" + seq[51:]
		out.write(",".join([snp.name,snp.chr,snp.seq]) + "\n")
	out.close()
	return 0

#############################################
# run the script
if __name__ == "__main__":
	main(sys.argv)
