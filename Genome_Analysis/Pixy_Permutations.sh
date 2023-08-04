#!/bin/bash

#SBATCH -J Pixy_Permutations
#SBATCH -N 1
#SBATCH -n 15
#SBATCH --mem=75G
#SBATCH -t 5-00:00:00
#SBATCH -p jro0014_amd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=npb0015@auburn.edu

module load htslib/1.11

# Loop through all autosomes and chrX to calculate empirical Fst and Dxy
# This uses spline defined window files and constitutes the analysis

for x in {1..20} X
do 
~/conda/bin/pixy --n_cores 15 --stats dxy fst --vcf Reinforcement.chr${x}.filtered.allsites.vcf.gz --populations Reinforcement_Populations.txt --bed_file Reinforcement.chr${x}.MAF.SplineWindows.bed --output_prefix Reinforcement.SplineDefined.chr${x}
done

# Start a different loop to randomize sample list into different populations to generate permutations

for x in {1..100}
do
sort -R Sample_List.txt > Randomized_List.txt
paste Randomized_List.txt Reinforcement_Populations.txt | awk '{OFS="\t" ; print $1,$3}' > Randomized_Populations_${x}.txt

# Within this loop, conduct calculations again as above

for y in {1..20} X
do
~/conda/bin/pixy --n_cores 15 --stats dxy fst --vcf Reinforcement.chr${y}.filtered.allsites.vcf.gz --populations Randomized_Populations_${x}.txt --bed_file Reinforcement.chr${y}.MAF.SplineWindows.bed --output_prefix Randomized.SplineDefined.chr${y}.${x}
done
done
