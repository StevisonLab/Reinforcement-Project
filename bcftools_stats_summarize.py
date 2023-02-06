#!/usr/bin/env python3

# Script written by Nick Bailey

import argparse
import sys
import os
import glob
import re
import math

#def ObtainValue(Statement):
#        return Statement[0].split('\t')

def main():

# Define arguments from argparse module, the purpose of which are specified in the below lines

	parser = argparse.ArgumentParser()
	parser.add_argument('-P','--path', help='Path to directory with bcftools stats files. If no argument is given this script will search current directory.', required=False)
	parser.add_argument('-S','--samples', help='Path to file with sample names listed. This is required or else the script can\'t sum values per sample.', required=True)
	parser.add_argument('-O','--output', help='Name of output file. If no argument is given this script will produce/overwrite FastQC_Summary.txt.', required=False)
	args = parser.parse_args()

	Samples = open(args.samples, 'r')

# If a path is given then change to that directory

	if args.path:
		os.chdir(args.path)

# Put all stats.txt files into an array.

	Files = sorted(glob.glob('*.stats.txt'))

# If there are no files then give relevant error and suggestion.

	if not Files:
		print('There appear to be no bcftools stats files here. Either move this script to a directory with these files or specify path using -P (or --path) option.')
	else:

# If output file is specified then use that, otherwise use default

		if args.output:
			Output = open(args.output, 'w')
		else:
			Output = open('BCFtools_Summary.txt', 'w')

# Write out the header for each category to the file

		Output.write('sample	RefHom	NonRefHom	Hets	Transitions	Transversions	Singletons	Missing' + '\n')

# Loop through the stats.txt files, make a name for sequence reads and then read in file
# Then define the relevant fastqc information and output it in file
# the "try" block handles fastqc version 0.11.8 data while the "except" block excludes data not available in version 0.10.1
		for Sample in Samples:

			Sample = Sample.strip('\n')
			SampleData = []
			RefHom = []
			NonRefHom = []
			Hets = []
			Transitions = []
			Transversions = []
			Singletons = []
			Missing = []
			for File in Files:

				Stats = open(File).read().rstrip()
				PSC = re.findall('PSC\t[0-9]*\t' + Sample + '\t[0-9]*\t[0-9]*\t[0-9]*\t[0-9]*\t[0-9]*\t[0-9]*\t[0-9]*.[0-9]\t[0-9]*\t[0-9]*\t[0-9]*\t[0-9]*', Stats)
				try:
					SampleData.append(PSC[0].split('\t'))
				except IndexError:
					print('Sample', Sample, 'does not exist in', File)

			for Data in SampleData:
				RefHom.append(int(Data[3]))
				NonRefHom.append(int(Data[4]))
				Hets.append(int(Data[5]))
				Transitions.append(int(Data[6]))
				Transversions.append(int(Data[7]))
				Singletons.append(int(Data[10]))
				Missing.append(int(Data[13]))
			

			Output.write(Sample + '\t' + str(sum(RefHom)) + '\t' + str(sum(NonRefHom)) + '\t' + str(sum(Hets)) + '\t' + str(sum(Transitions)) + '\t' + str(sum(Transversions)) + '\t' + str(sum(Singletons)) + '\t' + str(sum(Missing)) + '\n')

		Output.close()

if __name__ == "__main__":
	main()
