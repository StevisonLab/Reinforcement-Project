#!/usr/bin/env python3

# Script written by Nick Bailey

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

# Put all flagstat.txt files into an array.
# This assumes uncompressed fastqc folders with standard naming.

	Files = sorted(glob.glob('*flagstat.txt'))

# If there are no files then give relevant error and suggestion.

	if not Files:
		print('There appear to be no flagstat files here. Either move this script to a directory with these files or specify path using -P (or --path) option.')
	else:

# If output file is specified then use that, otherwise use default

		if args.output:
			Output = open(args.output, 'w')
		else:
			Output = open('Flagstat_Summary.txt', 'w')

# Write out the header for each category to the file

		Output.write('Sample	Mapped	Duplicated' + '\n')

# Loop through the flagstat.txt files, make a name for sequence reads and then read in file
# Then define the relevant flagstat information and output it in file

		for File in Files:

			Sample = File.replace('.flagstat.txt','')
			Flagstat = open(File).read().rstrip()

			Total = float(obtain_value(re.findall('[0-9]* \+ 0 in total', Flagstat)))
			Duplicates = float(obtain_value(re.findall('[0-9]* \+ 0 duplicates', Flagstat)))
			DupPercent = re.findall('[0-9][0-9].[0-9][0-9]', str((Duplicates/Total)*100))[0]

			Mapped = re.findall('mapped \([0-9][0-9].[0-9][0-9]', Flagstat)[0].split(' ')[1].replace('(','')

			Output.write(Sample + '\t' + Mapped + '\t' + DupPercent + '\n')

		Output.close()

if __name__ == "__main__":
	main()
