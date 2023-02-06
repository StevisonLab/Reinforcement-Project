#!/usr/bin/env python3

# Script written by Nick Bailey

import numpy
import argparse
import os
import glob
import re

def obtain_value(Statement):
        return Statement[0].split(' ')[0]

def main():

# Define arguments from argparse module, the purpose of which are specified in the below lines

	parser = argparse.ArgumentParser()
	parser.add_argument('-P','--path', help='Path to directory with flagstat files. If no argument is given this script will search current directory.', required=False)
	parser.add_argument('-O','--output', help='Name of output file. If no argument is given this script will produce/overwrite Flagstat_Summary.txt.', required=False)
	args = parser.parse_args()

# If a path is given then change to that directory

	if args.path:
		os.chdir(args.path)

# Put all coverage.txt files into an array.
# This assumes uncompressed fastqc folders with standard naming.

	Files = sorted(glob.glob('*coverage.txt'))

# If there are no files then give relevant error and suggestion.

	if not Files:
		print('There appear to be no flagstat files here. Either move this script to a directory with these files or specify path using -P (or --path) option.')
	else:

# If output file is specified then use that, otherwise use default

		if args.output:
			Output = open(args.output, 'w')
		else:
			Output = open('Coverage_Summary.txt', 'w')

# Write out the header for each category to the file

		Output.write('Sample	Coverage	Depth	Quality' + '\n')

# Loop through the coverage.txt files, make a name for sequence reads and then read in file
# Then define the relevant coverage information and output it in file

		for File in Files:

			Sample = File.replace('.coverage.txt','')
			CoverageFile = open(File).read().rstrip()
			Lines = CoverageFile.split('\n')
			Length = []
			for Line in Lines:
				Chromosome = Line.split('\t')
				try:
					Length.append(int(Chromosome[2]))
				except ValueError:
					NULL="NULL"
			Lengths = numpy.array(Length, numpy.int)
			TotalLength = numpy.sum(Lengths)

			CoverageFile = open(File, 'r')

			WeightedCov=[]
			WeightedDepth=[]
			WeightedQual=[]

			for Line in CoverageFile:			
				Line.strip('\n')
				Chromosome = Line.split('\t')
				try:
					Length = int(Chromosome[2])
					Percent = Length/TotalLength
					Coverage = float(Chromosome[5])
					Depth = float(Chromosome[6])
					BaseQ = float(Chromosome[7])
					WeightedCov.append(Percent * Coverage)
					WeightedDepth.append(Percent * Depth)
					WeightedQual.append(Percent * BaseQ)

				except ValueError:
					NULL="NULL"

			MeanCov = numpy.array(WeightedCov, numpy.float)
			MeanDepth = numpy.array(WeightedDepth, numpy.float)
			MeanQual = numpy.array(WeightedQual, numpy.float)

			Output.write('%s	%.2f	%.2f	%.2f' % (Sample, numpy.sum(MeanCov), numpy.sum(MeanDepth), numpy.sum(MeanQual)) + '\n')

		Output.close()

if __name__ == "__main__":
	main()
