#!/bin/bash

# This script exists to parse the Q file output of ADMIXTURE, which is represented as simply two columns of values with no labels

# First few steps as in ADMIXTURE.sh

for chrom in {1..20} X
do

STARTS=(`awk '{print $2}' Reinforcement.chr${chrom}.corrected.SplineWindows.bed`)

for start in ${STARTS[@]}
do

new_start=`expr ${start} + 1`
end=`awk -v START=${start} '$2==START {print $3}' Reinforcement.chr${chrom}.corrected.SplineWindows.bed`

# If ADMIXTURE output doesn't exist, end with error message

if [ ! -s chr${chrom}.${new_start}.${end}.2.Q ]
then
echo "ADMIXTURE output doesn't exist"
else

# If file containing individuals removed from analysis exists

if [ -s chr${chrom}.${new_start}.${end}.irem ]
then

# Then cut a single column out for new individual listing file

cut -f 1 chr${chrom}.${new_start}.${end}.irem > chr${chrom}.${new_start}.${end}.new.irem

# Exclude missing individuals from pop file and then search a random hard-coded rhesus

number=`grep -v -f chr${chrom}.${new_start}.${end}.new.irem autosome.ind.pop.txt | grep -n "SRS2588667" | awk -F: '{print $1}'`
else

# If no individuals removed from analysis, use random hard-coded rhesus

number=`grep -n "SRS2588667" autosome.ind.pop.txt | awk -F: '{print $1}'`
fi

# Then isolate first and second column of data file and multiply values by 100000 to get integers

first=`awk -v NUMBER=${number} 'NR == NUMBER {print $1}' chr${chrom}.${new_start}.${end}.2.Q`
second=`awk -v NUMBER=${number} 'NR == NUMBER {print $2}' chr${chrom}.${new_start}.${end}.2.Q`
First=`echo "${first} * 100000" | bc | awk -F. '{print $1}'`
Second=`echo "${second} * 100000" | bc | awk -F. '{print $1}'`

# Given that row, isolate which column has a higher value, that should represent rhesus ancestry proportion

if [ ${First} -gt ${Second} ]
then 
column=1
else
column=2
fi

# Isolate row numbers now for all ParaCyno individuals

if [ -s chr${chrom}.${new_start}.${end}.irem ]
then
numbers=(`grep -v -f chr${chrom}.${new_start}.${end}.new.irem autosome.ind.pop.txt | grep -n "ParaCyno" | awk -F: '{print $1}'`)
else
numbers=(`grep -n "ParaCyno" autosome.ind.pop.txt | awk -F: '{print $1}'`)
fi

# Loop through row numbers and isolate ParaCyno rows and column of rhesus ancestry

for number in ${numbers[@]}
do
awk -v COLUMN=${column} -v NUMBER=${number} 'NR == NUMBER {print $COLUMN}' chr${chrom}.${new_start}.${end}.2.Q >> chr${chrom}.${new_start}.${end}.ParaCyno.2.Q
done

# If this file is successfully generated, then average rhesus ancestry over all ParaCyno individuals

if [ -s chr${chrom}.${new_start}.${end}.ParaCyno.2.Q ]
then
rhesus_admix=`awk '{sum+=$1} END {print sum/NR}' chr${chrom}.${new_start}.${end}.ParaCyno.2.Q`
echo -e "chr${chrom}\t${new_start}\t${end}\t${rhesus_admix}" >> ParaCyno.admixture.txt
else
echo "No ParaCyno samples left"
fi

fi
done

done
