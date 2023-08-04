#!/usr/bin/env python3

# Import relevant modules
# Notably, path was appended to use scikit-allel functions

from sys import path
import numpy
from pandas import read_csv
path.append('/home/npb0015/conda/pkgs/scikit-allel-1.3.5-py38h43a58ef_1/lib/python3.8/site-packages')
from allel import read_vcf
from allel import GenotypeArray
from allel import windowed_tajima_d

# Define input using pandas function

Populations = read_csv('../Reinforcement_Populations.txt', sep = '\t', names = ['samples','populations'])

# Open new file for output

Output = open('TajimaD.Output.txt', 'a')

# Loop through all chromosomes and all populations and begin by opening relevant chromosome spline defined windows file

for CHROM in tuple(range(1,21)) + tuple('XY'):

	for Population in Populations.populations.unique():

		Windows = open('../Reinforcement.chr' + str(CHROM) + '.MAF.SplineWindows.bed', 'r')

# For each window, read in data in format appropriate for scikit-allel and calculate size of window

		for Line in Windows:

			Line = Line.strip('\n')
			Fields = Line.split('\t')

			VCF = read_vcf('../Reinforcement.chr' + str(CHROM) + '.filtered.allsites.vcf.gz', samples = Populations[Populations.populations == Population ]['samples'], region = Fields[0] + ':' + Fields[1] + '-' + Fields[2])

			Size = int(Fields[2]) - int(Fields[1])

# Try to calculate Tajima's D  in window and output only the Tajima's D value

			try:

				D_calc = windowed_tajima_d(VCF['variants/POS'], GenotypeArray(VCF['calldata/GT']).count_alleles(), size = Size, start = int(Fields[1]), stop = int(Fields[2]))[0][0]

# If this fails, simply report the value as non-existent

			except TypeError:

				D_calc = 'nan'

# Write out the data to Output file defined near beginning and close out loop then close out output file

			Output.write(Line + '\t' + Population + '\t' + str(D_calc) + '\n')

	Windows.close()

Output.close()
