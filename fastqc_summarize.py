#!/usr/bin/env python3

# Script written by Nick Bailey

import argparse
import os
import glob
import re

def obtain_value(Statement):
        return Statement[0].split('\t')[1]

def main():

# Define arguments from argparse module, the purpose of which are specified in the below lines

	parser = argparse.ArgumentParser()
	parser.add_argument('-P','--path', help='Path to directory with FastQC directories. If no argument is given this script will search current directory.', required=False)
	parser.add_argument('-O','--output', help='Name of output file. If no argument is given this script will produce/overwrite FastQC_Summary.txt.', required=False)
	args = parser.parse_args()

# If a path is given then change to that directory

	if args.path:
		os.chdir(args.path)

# Put all fastqc_data.txt files into an array.
# This assumes uncompressed fastqc folders with standard naming.

	Files = sorted(glob.glob('*_fastqc/fastqc_data.txt'))

# If there are no files then give relevant error and suggestion.

	if not Files:
		print('There appear to be no FastQC directories here. Either move this script to a directory with these files or specify path using -P (or --path) option.')
	else:

# If output file is specified then use that, otherwise use default

		if args.output:
			Output = open(args.output, 'w')
		else:
			Output = open('FastQC_Summary.txt', 'w')

# Write out the header for each category to the file

		Output.write('Read	SeqNumber	NumPoorQualSeq	SeqLength	PerBaseQual	OverRepSeq	Adapter' + '\n')

# Loop through the fastqc_data.txt files, make a name for sequence reads and then read in file
# Then define the relevant fastqc information and output it in file
# the "try" block handles fastqc version 0.11.8 data while the "except" block excludes data not available in version 0.10.1

		for File in Files:

			ReadName = File.replace('_fastqc/fastqc_data.txt','')
			FastQC = open(File).read().rstrip()
			try:
				SeqNumber = obtain_value(re.findall('Total Sequences\t[0-9]*', FastQC))
				PoorQual = obtain_value(re.findall('Sequences flagged as poor quality\t[0-9]*', FastQC))
				SeqLength = obtain_value(re.findall('Sequence length\t[0-9]*', FastQC))
				OverRep = obtain_value(re.findall('>>Overrepresented sequences\t[a-z]*', FastQC))
				Adapter = obtain_value(re.findall('>>Adapter Content\t[a-z]*', FastQC))
				BaseQual = obtain_value(re.findall('>>Per base sequence quality\t[a-z]*', FastQC))
				Output.write(ReadName + '\t' + SeqNumber + '\t' + PoorQual + '\t' + SeqLength + '\t' + BaseQual + '\t' + OverRep + '\t' + Adapter + '\n')
			except IndexError:
				SeqNumber = obtain_value(re.findall('Total Sequences\t[0-9]*', FastQC))
				SeqLength = obtain_value(re.findall('Sequence length\t[0-9]*', FastQC))
				OverRep = obtain_value(re.findall('>>Overrepresented sequences\t[a-z]*', FastQC))
				BaseQual = obtain_value(re.findall('>>Per base sequence quality\t[a-z]*', FastQC))
				Output.write(ReadName + '\t' + SeqNumber + '\t' + 'NA' + '\t' + SeqLength + '\t' + BaseQual + '\t' + OverRep + '\t' + 'NA' + '\n')

		Output.close()

if __name__ == "__main__":
	main()
