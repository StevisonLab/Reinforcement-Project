#!/bin/bash

cd SNP_vcfs

Files=(Reinforcement.chrX.filtered.allsites.vcf.gz Reinforcement.chrY.filtered.allsites.vcf.gz)

for File in ${Files[@]}
do

PREFIX=`echo ${File} | sed 's/.filtered.allsites.vcf.gz//'`

bcftools filter --threads 6 -i 'TYPE="snp" && QUAL!="."' ../${File} -O z -o ${PREFIX}.filtered.variant.vcf.gz 

bcftools +fill-tags ${PREFIX}.filtered.variant.vcf.gz -- -t MAF > ${PREFIX}.MAF.vcf

crabz -p 6 -l 9 ${PREFIX}.MAF.vcf > ${PREFIX}.MAF.vcf.gz

rm ${PREFIX}.MAF.vcf

bcftools query -f '%CHROM\t%POS\t%MAF\n' ${PREFIX}.MAF.vcf.gz > ${PREFIX}.MAF.txt

crabz -p 6 -l 9 ${PREFIX}.MAF.txt > ${PREFIX}.MAF.txt.gz

rm ${PREFIX}.MAF.txt

done
