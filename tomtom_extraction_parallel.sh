#!/bin/bash

#SBATCH -J parallel_extraction
#SBATCH -N 1
#SBATCH -n 15
#SBATCH --mem=105G
#SBATCH -t 20:00:00
#SBATCH -p jro0014_amd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=npb0015@auburn.edu

module load snpeff/5.0
module load gnu-parallel/20120222

java -jar /tools/snpeff-5.0/SnpSift.jar extractFields ParaCyno.Unique.motif.vcf CHROM POS "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATUREID" "ANN[*].EFFECT" > Variant_Info.txt

# Create header for final output file.
# This tracks the order of chromosome, position, HGNC gene name, Ensembl gene name,
# Ensembl transcript name, SnpEff effect, and tomtom predicted transcription factor,
# but there can be multiple genes and other things so don't assume the field separation of genes matches
# that of the header. The header is simply a rough guide of the order
# and realistically the script should probably be parsed in rows 
# (e.g. a gene can be searched for by row)  

echo -e "CHROM\tSTART\tEND\tHGNC\tGENE\tTRANSCRIPT\tEFFECT\tTF" > Variant_Motifs.bed

# Define sequence names from FASTA file

Seqs=(`grep ">" ParaCyno.Unique.Ref.fa | sed 's/>//'`)

parallel '

# Define chromosome and site location from sequence name

chr=`echo {} | cut -d ":" -f 1`
site_inter=`echo {} | cut -d ":" -f 2 | cut -d "-" -f 1`
site=`expr ${site_inter} + 7`
pre_site=`expr ${site_inter} + 6`

# Find entry in VCF file and output gene name using SnpSift

Fields=(`grep "${chr}	${site}" Variant_Info.txt`)

# Define motifs in ref and alt tomtom output

ref_motifs=(`grep "^chr" tomtom_output/{}.Ref_out/tomtom.tsv | cut -f 2`)
alt_motifs=(`grep "^chr" tomtom_output/{}.Alt_out/tomtom.tsv | cut -f 2`)

# Identify motifs unique to ref or alt

ref_uniq=(`echo ${ref_motifs[@]} ${alt_motifs[@]} | sed "s/ /\n/g" | sort | uniq -d | xargs echo ${ref_motifs[@]} | sed "s/ /\n/g" | sort | uniq -u`)
alt_uniq=(`echo ${ref_motifs[@]} ${alt_motifs[@]} | sed "s/ /\n/g" | sort | uniq -d | xargs echo ${alt_motifs[@]} | sed "s/ /\n/g" | sort | uniq -u`)

Trans_count=$(((${#Fields[@]} - 2) / 4))

# Loop through motifs uniq to alt and find gene in meme database

	for motif in ${alt_uniq[@]}
	do
	Fields[${#Fields[@]}]=`grep "${motif}" ~/meme/motif_databases/CIS-BP_2.00/Macaca_fascicularis.meme | cut -d "(" -f 2 | sed "s/)_//"`
	done

HGNC=`printf "%s," "${Fields[@]:2:${Trans_count}}" | sed "s/,$//"`
GENE=`printf "%s," "${Fields[@]:2+${Trans_count}:${Trans_count}}" | sed "s/,$//"`
TRANSCRIPT=`printf "%s," "${Fields[@]:2+2*${Trans_count}:${Trans_count}}" | sed "s/,$//"`
EFFECT=`printf "%s," "${Fields[@]:2+3*${Trans_count}:${Trans_count}}" | sed "s/,$//"`
TF=`printf "%s," "${Fields[@]:2+4*${Trans_count}:${#alt_uniq[@]}}" | sed "s/,$//"`

echo -e "${Fields[0]}\t${pre_site}\t${Fields[1]}\t${HGNC}\t${GENE}\t${TRANSCRIPT}\t${EFFECT}\t${TF}" >> Variant_Motifs.txt' ::: ${Seqs[@]}

sort -k1,1 Variant_Motifs.txt >> Variant_Motifs.bed
rm Variant_Motifs.txt
