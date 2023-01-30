#!/usr/bin/env python3

# Read vcf file
# Delimit columns by tab
# Then delimit sample columns by colon
# Then isolate first value, which is the genotype
# Then replace phase bars with slashes
# Then split genotypes into individual alleles

import subprocess

VCF = open('Outlier.Genes.vcf', 'r')
subprocess.run('grep "^#" Outlier.Genes.vcf > ParaCyno.Unique.vcf', shell = True)
Output = open('ParaCyno.Unique.vcf', 'a')

for Line in VCF:
	Line = Line.strip('\n')
	Column = Line.split('\t')
	if "#" not in Column[0]:
		DRS139839 = Column[11].split(':')[0].replace('|','/').split('/')
		FascA12133 = Column[12].split(':')[0].replace('|','/').split('/')
		FascA15865 = Column[13].split(':')[0].replace('|','/').split('/')
		FascA17187 = Column[14].split(':')[0].replace('|','/').split('/')
		FascA18848 = Column[15].split(':')[0].replace('|','/').split('/')
		FascA19211 = Column[16].split(':')[0].replace('|','/').split('/')
		FascB00725 = Column[17].split(':')[0].replace('|','/').split('/')
		FascB00729 = Column[18].split(':')[0].replace('|','/').split('/')
		FascB00742 = Column[19].split(':')[0].replace('|','/').split('/')
		FascB2232 = Column[20].split(':')[0].replace('|','/').split('/')
		SRS117874 = Column[21].split(':')[0].replace('|','/').split('/')
		SRS212016 = Column[22].split(':')[0].replace('|','/').split('/')
		SRS2588651 = Column[23].split(':')[0].replace('|','/').split('/')
		SRS2588654 = Column[24].split(':')[0].replace('|','/').split('/')
		SRS2588667 = Column[25].split(':')[0].replace('|','/').split('/')
		SRS2588673 = Column[26].split(':')[0].replace('|','/').split('/')
		SRS3755433 = Column[27].split(':')[0].replace('|','/').split('/')
		SRS4048144 = Column[28].split(':')[0].replace('|','/').split('/')
		SRS4048145 = Column[29].split(':')[0].replace('|','/').split('/')
		SRS5035087 = Column[30].split(':')[0].replace('|','/').split('/')
		SRS5035100 = Column[31].split(':')[0].replace('|','/').split('/')
		SRS5035120 = Column[32].split(':')[0].replace('|','/').split('/')
		SRS5035173 = Column[33].split(':')[0].replace('|','/').split('/')
		SRS5035194 = Column[34].split(':')[0].replace('|','/').split('/')
		SRS674436 = Column[35].split(':')[0].replace('|','/').split('/')
		SRS674437 = Column[36].split(':')[0].replace('|','/').split('/')
		SRS693707 = Column[37].split(':')[0].replace('|','/').split('/')
		SRS693712 = Column[38].split(':')[0].replace('|','/').split('/')
		SRS693714 = Column[39].split(':')[0].replace('|','/').split('/')
		SRS748669 = Column[40].split(':')[0].replace('|','/').split('/')
		SRS778761 = Column[41].split(':')[0].replace('|','/').split('/')
		SRS791155 = Column[42].split(':')[0].replace('|','/').split('/')
		SRS823980 = Column[43].split(':')[0].replace('|','/').split('/')
		SRS823981 = Column[44].split(':')[0].replace('|','/').split('/')
		SRS837317 = Column[45].split(':')[0].replace('|','/').split('/')
		SRS886329 = Column[46].split(':')[0].replace('|','/').split('/')
		SRS886338 = Column[47].split(':')[0].replace('|','/').split('/')
		SRS893058 = Column[48].split(':')[0].replace('|','/').split('/')

# Define ParaCyno has a set of nonduplicate alleles present in ParaCyno individuals
# Define Others the same way but with all other individuals

		ParaCyno = set(FascA15865 + FascA17187 + FascA19211 + FascB00725 + FascB00729 + FascB00742 + SRS117874 + SRS4048144 + SRS4048145 + FascA15865)
		Others = set(DRS139839 + FascA12133 + FascA18848 + FascB2232 + SRS212016 + SRS2588651 + SRS2588654 + SRS2588667 + SRS2588673 + SRS3755433 + SRS5035087 + SRS5035100 + SRS5035120 + SRS5035173 + SRS5035194 + SRS674436 + SRS674437 + SRS693707 + SRS693712 + SRS693714 + SRS748669 + SRS778761 + SRS791155 + SRS823980 + SRS823981 + SRS837317 + SRS886329 + SRS886338 + SRS893058)

# If statement checks that the set of alleles in ParaCyno doesn't equal those shared between ParaCyno and Others
# That is, if there are alleles present in ParaCyno not in any individual outside that group

		if ParaCyno != ParaCyno & Others:

# Define specific alleles unique to ParaCyno set

			ParaCynoUniq = ParaCyno - Others

# If the "allele" unique to ParaCyno isn't just a missing symbol, then continue

			if ParaCynoUniq != {'.'}:

# For the ParaCyno population, count and sum the number of the unique allele
# Then output a line matching the input VCF but with a "Rank" appended with count of the allele
# Do NOT forget that python list slices are NOT inclusive of last number
# Rank should range within range of [1,18].
# 1 being a minimum because these should all have 1 allele unique to ParaCyno
# 18 being a maximum because there are 9 diploid individual ParaCynos
# Presently some are ranked 0 because there could be multiple unique alleles that are inappropriately concatenated together causing them to be treated as one allele not present, hopefully I will fix this soon and remove this comment
# "Fixed" the above error in that now counts can be done for sites with multiple ParaCyno Unique alleles but now if missing sites AND some actual allele are unique to ParaCyno, it'll count both
# Site chr1 110669495 is an example of this
# Pretty sure I fixed the above issues but now the script is convoluted so hopefully can edit later

				NumAlleles = 0
				for Individual in Column[13:15] + Column[16:20] + Column[21:22] + Column[28:30]:
					if len(ParaCynoUniq) == 1:
						GT = Individual.split(':')[0]
						NumAlleles += GT.count(''.join(ParaCynoUniq))
					elif len(ParaCynoUniq) > 1:
						if '.' in ParaCynoUniq:
							ParaCynoUniq = list(ParaCynoUniq)
							ParaCynoUniq.remove('.')
						AlleleCount = len(ParaCynoUniq)
						GT = Individual.split(':')[0]
						for Allele in range(0,AlleleCount):
							NumAlleles += GT.count(''.join(list(ParaCynoUniq)[Allele]))
				Output.write('\t'.join(Column[0:7]) + '\t' + Column[7] + ';RANK=' + str(NumAlleles) + '\t' + '\t'.join(Column[8:49]) + '\n')

VCF.close()
Output.close()
