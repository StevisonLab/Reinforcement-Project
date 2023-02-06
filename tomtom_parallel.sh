#!/bin/bash

#SBATCH -J parallel_tomtom
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=100G
#SBATCH -t 24:00:00
#SBATCH -p jro0014_amd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=npb0015@auburn.edu

module load gnu-parallel/20120222

# Define array of sequence names (represented simply as chr:start-end) in FASTA
# Names are the same in both Ref and Alt files

Seqs=(`grep ">" ParaCyno.Unique.Ref.fa | sed 's/>//'`)

# Loop through sequence names and grep the relevant sequences from Ref and Alt
# into individual sequence files

for Seq in ${Seqs[@]}
do 
grep -A 1 "${Seq}" ParaCyno.Unique.Ref.fa > tomtom_output/${Seq}.Ref.fa
grep -A 1 "${Seq}" ParaCyno.Unique.Alt.fa > tomtom_output/${Seq}.Alt.fa
done

# Define array of all FASTA files generated above

Files=(`ls tomtom_output/*.fa`)

# Parallelize converting the FASTA sequences to PWMs with meme-suite

# And depending on whether sequence is ref or alt, use rhesus or cyno database respectively for tomtom

parallel '
Prefix=`echo {} | sed "s/.fa//"`

~/meme/libexec/meme-5.4.1/rna2meme -dna {} > ${Prefix}.meme

if [[ ${Prefix} == *"Ref" ]]
then
~/meme/bin/tomtom -verbosity 4 ${Prefix}.meme ~/meme/motif_databases/CIS-BP_2.00/Macaca_mulatta.meme -o ${Prefix}_out
else
~/meme/bin/tomtom -verbosity 4 ${Prefix}.meme ~/meme/motif_databases/CIS-BP_2.00/Macaca_fascicularis.meme -o ${Prefix}_out
fi' ::: ${Files[@]}

