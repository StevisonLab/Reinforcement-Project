#!/bin/bash

for chrom in {1..20} X
do

# Obtain start positions for each window on chromosome and loop through

STARTS=(`awk '{print $2}' Reinforcement.chr${chrom}.corrected.SplineWindows.bed`)

for start in ${STARTS[@]}
do

# Create new start matching VCF format (+1 of bed format)

new_start=`expr ${start} + 1`

# Isolate end coordinate for equivalent start coordinate

end=`awk -v START=${start} '$2==START {print $3}' Reinforcement.chr${chrom}.corrected.SplineWindows.bed`

# Isolate VCF for these coordinates

bcftools view --threads 8 -r chr${chrom}:${new_start}-${end} ../pca/Reinforcement.chr${chrom}.filtered.variant.ID.vcf.gz -O z -o chr${chrom}.${new_start}.${end}.vcf.gz

# Create appropriate plink formatted files for these

plink --threads 8 --vcf chr${chrom}.${new_start}.${end}.vcf.gz --keep Sample_List.txt --geno 0.999999 --mind 0.999999 --make-bed --out chr${chrom}.${new_start}.${end}

# Run ADMIXTURE with 8 threads and K value of 2

admixture -j8 chr${chrom}.${new_start}.${end}.bed 2 > chr${chrom}.${new_start}.${end}.log
done

done
