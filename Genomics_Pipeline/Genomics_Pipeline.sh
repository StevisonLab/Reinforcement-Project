#!/bin/bash

#SBATCH -J Genomics_Pipeline
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=64G
#SBATCH -t 2-00:00:00
#SBATCH -p general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=npb0015@auburn.edu

# Load all relevent modules

module load fastqc/0.11.9
module load samtools/1.11
module load gatk/4.1.9.0
module load R/4.0.3
module load bwa/0.7.17
module load sratoolkit/2.11.0
module load bcftools/1.11
module load python/3.9.2

# Define arrays for the samples file of interest and then extract reads from this file

Sample_Data=~/reinforcement_project/Reference_files/Master_SraRunInfo.csv
my_reads=(`awk -F, '{print $1}' ${Sample_Data}`)

# Loop through all reads then download and compress them if they don't already exist

cd fastq_files/

for i in ${my_reads[@]}
do
	if [ -s "${i}_1.fastq" ] && [ -s "${i}_2.fastq" ] || [ -s "${i}_1.fastq.gz" ] && [ -s "${i}_2.fastq.gz" ]
	then
	echo "Download of ${i} has already occurred"
	else
	fasterq-dump -S -e 8 ${i}
	~/pigz-2.7/pigz -9 -p 8 ${i}_1.fastq
	~/pigz-2.7/pigz -9 -p 8 ${i}_2.fastq
	fi
done

	if [ -s "${i}_1.fastq.gz" ] && [ -s "${i}_2.fastq.gz" ] && [ ! -s ${i}_1_fastqc.zip ] && [ ! -s ${i}_2_fastqc.zip ] # If files were successfully downloaded and zipped then run a quality check
	then
	fastqc -t 8 ${i}_1.fastq.gz
	fastqc -t 8 ${i}_2.fastq.gz
	else
	echo "There are no gzipped fastq files to check or fastqc is already complete"
	fi

cd ..

# Add reads from our newly sequenced samples

my_reads+=(`ls ~/reinforcement_project/Novogene_Samples/usftp21.novogene.com/raw_data/Fasc*/*fq.gz | awk -F/ '{print $9}' | sort -u | sed "s/_[1-2].fq.gz//"`)

for i in ${my_reads[@]}
do

	sample=`awk -F, -v SRR=${i} '$1==SRR {print $25}' ${Sample_Data}` # Define variables necessary for bwa-mem alignment
	PU=`awk -F, -v SRR=${i} '$1==SRR {print $26}' ${Sample_Data}`
	LB=`awk -F, -v SRR=${i} '$1==SRR {print $12}' ${Sample_Data}`


if [ -z ${sample} ] && [ -z ${PU} ] && [ -z ${LB} ] # these paremeters should be missing for newly sequences samples
then

# These commands are necessary to extract data from the files according to paremeters given by Table 1 of the Novogene Data Summary File

	sample=`echo ${i} | awk -F_ '{print $1}'`
	PU=`zcat ~/reinforcement_project/Novogene_Samples/usftp21.novogene.com/raw_data/${sample}/${i}_1.fq.gz | head -n 1 | awk -F: '{OFS="." ; print $3,$4}'`
	LB=`echo ${i}_1.fq.gz | awk -F_ '{OFS="_" ; print $2,$4}'`

fi

	if [ -z ${LB} ] && [ ! -z ${PU} ] && [ ! -z ${PU} ] # If Library doesn't exist in the RunInfo file, then just make it "unknown" so it can still be aligned
	then
	LB="unknown"
	fi

	if [ -s "fastq_files/${i}_1.fastq.gz" ] && [ -s "fastq_files/${i}_2.fastq.gz" ] && [ ! -s "samtools_files/${sample}.${i}.sam" ] # If the fastq.gz files exist and they haven't been aligned to a reference genome, then do that
	then
	bwa mem -M -v 2 -t 8 -R "@RG\tID:${i}\tSM:${sample}\tPU:${PU}\tPL:Illumina\tLB:${LB}" ~/reinforcement_project/Reference_files/rheMac10.masked.female.fa fastq_files/${i}_1.fastq.gz fastq_files/${i}_2.fastq.gz > samtools_files/${sample}.${i}.sam
	elif [ -s ~/reinforcement_project/Novogene_Samples/usftp21.novogene.com/raw_data/${sample}/${i}_1.fq.gz ] && [ -s ~/reinforcement_project/Novogene_Samples/usftp21.novogene.com/raw_data/${sample}/${i}_2.fq.gz ] && [ ! -s samtools_files/${sample}.${i}.sam ] # same as above but for newly sequenced samples
	then
	bwa mem -M -v 2 -t 8 -R "@RG\tID:${i}\tSM:${sample}\tPU:${PU}\tPL:Illumina\tLB:${LB}" ~/reinforcement_project/Reference_files/rheMac10.masked.female.fa ~/reinforcement_project/Novogene_Samples/usftp21.novogene.com/raw_data/${sample}/${i}_1.fq.gz ~/reinforcement_project/Novogene_Samples/usftp21.novogene.com/raw_data/${sample}/${i}_2.fq.gz > samtools_files/${sample}.${i}.sam
	else
	echo "${i} has already been aligned to reference"
	fi

	if [ -s "samtools_files/${sample}.${i}.sam" ] && [ ! -s "samtools_files/${sample}.${i}.bam" ] # If alignment was successful, compress the files to bam
	then
	samtools view -b -@ 8 samtools_files/${sample}.${i}.sam -o samtools_files/${sample}.${i}.bam
	else
	echo "${i} has already been compressed to bam"
	fi

	if [ -s "samtools_files/${sample}.${i}.bam" ] && [ ! -s "samtools_files/sorted.${sample}.${i}.bam" ] # If compression was successful, sort the bam files
	then
	samtools sort -@ 8 samtools_files/${sample}.${i}.bam -o samtools_files/sorted.${sample}.${i}.bam
	samtools index -@ 8 samtools_files/sorted.${sample}.${i}.bam
	else
	echo "bam for ${i} has already been sorted"
	fi

	if [ -s "samtools_files/sorted.${sample}.${i}.bam" ] && [ ! -s "samtools_files/sorted.${sample}.${i}.chrX.bam" ] # This step is necessary only when working with sex chromosomes
	then
        samtools view -@ 8 samtools_files/sorted.${sample}.${i}.bam chrX -o samtools_files/sorted.${sample}.${i}.chrX.bam # this can be done with chrX or chrY
	else
	echo "chrX bam for ${i} has already been generated"
	fi

done

Samples=(awk '{print $1}' Sample_List.txt) # Define SRS samples and then merge any reads that fall under a common sample
for x in ${Samples[@]}
do

	ls samtools_files/sorted.${x}*.bam > bam_list.txt
	Lines=`wc -l bam_list.txt | awk '{print $1}'`
	if [ ${Lines} -le 1 ] || [ -s sorted.${x}.combined.bam ] # Check that there are actually multiple reads for a given sample
	then
        echo "Sample ${x} has one or zero reads, or is already merged, no need for merging."
        else
        echo "Sample ${x} will be merged now."
	samtools merge samtools_files/sorted.${x}.combined.bam -@ 8 -b bam_list.txt
	fi

done

cd samtools_files/

sample_bams=(`ls sorted*.bam | sort -t . -k2,2 -u`) # Define bam files including he combined ones generated above and excluding ones where reads were separate
# The sort command depends on the combined files alphabetically coming before others (e.g. "combined" always come before "SRR") but would likely pick an incorrect file to sort unique for if this was not the case

# Run loops to build indices for each bam and mark duplicates

for f in ${sample_bams[@]}
do
prefix=`echo ${f} | sed 's/.bam//'`
if [ -s ${prefix}.bai ]
then
echo "Index for ${f} already exists"
else
gatk BuildBamIndex -I ${f} -R ~/reinforcement_project/Reference_files/rheMac10.masked.female.fa
fi
done

for f in ${sample_bams[@]}
do
sample=`echo ${f} | awk -F. '{print $2}'`
if [ -s ${sample}.markdup.bam ]
then
echo "${f} has already had duplicates marked"
else
gatk MarkDuplicates --java-options "-Xmx8G" --COMPRESSION_LEVEL 9 -I ${f} -O ${sample}.markdup.bam -M ${sample}.markdup.txt
fi
done

# Use GNU parallel to run BQSR steps for all samples

parallel '
if [ -s {}.based.table ]
then
echo "Bases have been recalibrated for {}"
else
gatk BaseRecalibrator --java-options "-Xmx8G" -I {}.markdup.bam -R ~/reinforcement_project/Reference_files/rheMac10.masked.fa --known-sites ~/reinforcement_project/Reference_files/mGAP.sorted.biallelic.rheMac10.vcf.gz --known-sites ~/reinforcement_project/Reference_files/81Chinese.sorted.rheMac10.vcf.gz -O {}.based.table
fi' ::: ${Samples[@]}

parallel '
if [ -s {}.postBQSR.bam ]
then
echo "BQSR has been applied to {}"
else
gatk ApplyBQSR --java-options "-Xmx8G" -R ~/reinforcement_project/Reference_files/rheMac10.masked.fa -I {}.markdup.bam --bqsr-recal-file {}.based.table -O {}.postBQSR.bam
fi' ::: ${Samples[@]}

parallel '
if [ -s {}.postBQSR.table ]
then
echo "postBQSR bases have been recalibrated for {}"
else
gatk BaseRecalibrator --java-options "-Xmx8G" -I {}.postBQSR.bam -R ~/reinforcement_project/Reference_files/rheMac10.masked.fa --known-sites ~/reinforcement_project/Reference_files/mGAP.sorted.biallelic.rheMac10.vcf.gz --known-sites ~/reinforcement_project/Reference_files/81Chinese.sorted.rheMac10.vcf.gz -O {}.postBQSR.table
fi' ::: ${Samples[@]}

# Run quality checks for bam files at each stage above, namely samtools flagstat and coverage

Files=(`ls *.bam | grep -E "postBQSR|markdup|sorted" | sed 's/.bam//'`)

parallel '
if [ -s {}.coverage.txt ]
then
echo "Coverage file for {} has already been generated"
else
samtools coverage {}.bam -o {}.coverage.txt
fi' ::: ${Files[@]}

for x in ${Files[@]}
do
if [ -s ${x}.flagstat.txt ]
then
echo "Flagstat file for ${x} has already been generated"
else
samtools flagstat -@ 8 ${x}.bam > ${x}.flagstat.txt
fi
done

# Run AnalyzeCovariates in a loop since I've had trouble parallelizing in the past and it goes pretty quickly

for a in ${Samples[@]}
do
if [ -s ${a}.postBQSR.pdf ]
then
echo "AnalyzeCovariates file for ${a} has already been generated"
else
gatk AnalyzeCovariates --java-options "-Xmx8G" -before ${a}.based.chrX.table -after ${a}.postBQSR.chrX.table -plots ${a}.postBQSR.chrX.pdf
fi
done

# Use GNU parallel for HaplotypeCaller and output MD5 text files for each

parallel '
if [ -s ../variant_files/{}.chrX.g.vcf.gz ]
then
echo "HaplotypeCaller has already been run on this file"
else
gatk --java-options "-Xmx8G" HaplotypeCaller -ploidy 2 -R ~/reinforcement_project/Reference_files/rheMac10.masked.female.fa -I {}.postBQSR.chrX.bam -O ../variant_files/{}.chrX.g.vcf.gz -OVM true -ERC GVCF
fi' ::: ${Samples[@]}

cd ../variant_files/

# Create sample map in appropriate format for GenomicsDBImport, this involves an awk and sed combo to remove suffixes from the first column

for y in ${Samples[@]}
do
echo ${y}.chrX.g.vcf.gz | awk '{OFS="\t" ; print $1"remove",$1}' | sed 's/.chrX.g.vcf.gzremove//' >> Reinforcement.sample_map
done

# Create Genomics Database for all chromosomes

gatk GenomicsDBImport --java-options "-Xmx8G" --genomicsdb-workspace-path ~/reinforcement_project/GenomicsDB --sample-name-map Reinforcement.sample_map --reader-threads 8 -L chr.list

rm Reinforcement.sample_map

# Define array of chromosomes and then run GenotypeGVCFs for each

chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20)

# Joint Genotype GVCFS represented in Genomics Database generated above for all chromosomes
# Ploidy should be switched to 1 for X and Y chromosomes in males

parallel 'gatk GenotypeGVCFs -ploidy 2 -R ~/reinforcement_project/Reference_files/rheMac10.masked.fa -V gendb://../../../../home/npb0015/reinforcement_project/GenomicsDB/ -all-sites -L {} -O Reinforcement.{}.female.jointcalls.vcf.gz' ::: ${chromosomes[@]}

# Combine all autosomal jointcall files (unnecessary when working with sex chromosomes)

bcftools concat -threads 8 Reinforcement.chr*.jointcalls.vcf.gz -O z -o Reinforcement.jointcalls.vcf.gz

# Steps below will filter invariant and variant sites by different standards

if [ ! -s Reinforcement.invariant.vcf.gz ]
then
bcftools filter --threads 8 -O z -i 'TYPE="ref" & FORMAT/DP>5 & FORMAT/RGQ>20' Reinforcement.jointcalls.vcf.gz -o Reinforcement.invariant.vcf.gz
else
echo "Invariant file already exists"
fi

if [ ! -s Reinforcement.variant.vcf.gz ]
then
bcftools filter --threads 8 -O z -i 'TYPE="snp"' Reinforcement.jointcalls.vcf.gz -o Reinforcement.variant.vcf.gz
else
echo "Variant file already exists"
fi

if [ -s Reinforcement.filtered.variant.vcf ]
then
echo "VQSR has already been conducted"
else

# VQSR is applied for variant filtering

gatk IndexFeatureFile -I Reinforcement.variant.vcf.gz

gatk VariantRecalibrator \
   -R ~/reinforcement_project/Reference_files/rheMac10.masked.fa \
   -V Reinforcement.variant.vcf.gz \
   --resource:mGAP,known=false,training=true,truth=true,prior=12.0 ~/reinforcement_project/Reference_files/mGAP.sorted.biallelic.rheMac10.vcf.gz \
   --resource:Xue,known=false,training=true,truth=false,prior=10.0 ~/reinforcement_project/Reference_files/81Chinese.sorted.rheMac10.vcf.gz \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
   --output Reinforcement.recal \
   --tranches-file Reinforcement.tranches \
   --rscript-file Reinforcement.VQSR.R

gatk ApplyVQSR \
   -R ~/reinforcement_project/Reference_files/rheMac10.masked.fa \
   -V Reinforcement.variant.vcf.gz \
   -O Reinforcement.filtered.variant.vcf \
   --recal-file Reinforcement.recal \
   --tranches-file Reinforcement.tranches \
   --truth-sensitivity-filter-level 99.0 \
   --add-output-vcf-command-line true \
   --verbosity INFO \
   --exclude-filtered true \
   -mode SNP
fi

# Combine and sort invariant and variant files

if [ -s Reinforcement.filtered.allsites.vcf.gz ]
then
echo "Invariant and variant files have already been combined and sorted"
else
bcftools concat --threads 8 Reinforcement.invariant.vcf.gz Reinforcement.filtered.variant.vcf -O z -o Reinforcement.allsites.vcf.gz
bcftools sort Reinforcement.allsites.vcf.gz -O z -o Reinforcement.filtered.allsites.vcf.gz
fi

# Compress bam files to cram for permanent storage
# Check these by using https://github.com/wtsi-hgi/bam2cram-check

for x in ${Samples[@]}
do

samtools view -T ~/reinforcement_project/Reference_files/rheMac10.masked.female.fa -C -@ 8 ../samtools_files/${x}.postBQSR.chrX.bam -o ../samtools_files/${x}.postBQSR.chrX.cram

python3 ~/bam2cram-check/main.py -b ../samtools_files/${x}.postBQSR.chrX.bam -c ../samtools_files/${x}.postBQSR.chrX.cram -e ../samtools_files/${x}.postBQSR.chrX.error --log ../samtools_files/${x}.postBQSR.chrX.log

done
