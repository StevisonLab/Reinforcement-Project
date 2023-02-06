#!/bin/bash

#SBATCH -J B-SIFT
#SBATCH -o %j
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=7G
#SBATCH -t 30:00:00
#SBATCH -p jro0014_amd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=npb0015@auburn.edu

# Define components of SIFT output file

CHROM=(`awk 'NR>1 {print $1}' ParaCyno.Unique.annotated_SIFTannotations.high.xls | sed 's/chr//'`)
POS=(`awk 'NR>1 {print $2}' ParaCyno.Unique.annotated_SIFTannotations.high.xls`)
REF=(`awk 'NR>1 {print $3}' ParaCyno.Unique.annotated_SIFTannotations.high.xls`)
TRANSCRIPT=(`awk 'NR>1 {print $5}' ParaCyno.Unique.annotated_SIFTannotations.high.xls`)
VarScore=(`awk 'NR>1 {print $13}' ParaCyno.Unique.annotated_SIFTannotations.high.xls`)
LineCount=`grep -c "^chr" ParaCyno.Unique.annotated_SIFTannotations.high.xls`

# Begin header line of output file

echo -e "Chromosome\tPosition\tTranscript\tB-SIFT" > B-SIFT.output.high.txt

# Loop through the line count of SIFT output file

for ((LINE=0; LINE<${LineCount}; LINE++))
do

# Do not bother with lines with NA SIFT scores

if [ ${VarScore[${LINE}]} != "NA" ]
then

# Isolate the reference allele SIFT score from the relevant chromosome file in the rheMac10 ref database

RefScore=`zless ../Mmul_10.99/${CHROM[${LINE}]}.gz | grep -m 1 "^${POS[${LINE}]}\s*${REF[${LINE}]}\s*${REF[${LINE}]}\s*${TRANSCRIPT[${LINE}]}" | awk '{print $11}'`

# Calculate difference of variant score and reference score, which is then a B-SIFT score
# Output that to file

Score=`echo ${VarScore[${LINE}]} - ${RefScore} | bc -l`

echo -e "chr${CHROM[${LINE}]}\t${POS[${LINE}]}\t${TRANSCRIPT[${LINE}]}\t${Score}" >> B-SIFT.output.high.txt

fi
done
