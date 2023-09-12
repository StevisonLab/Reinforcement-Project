# Genomic analysis of the rhesus macaque (_Macaca mulatta_) and the cynomolgus macaque (_Macaca fascicularis_) uncover polygenic signatures of reinforcement speciation

This repository contains data files and scripts to reproduce the analysis for the preprint titled: "Genomic analysis of the rhesus macaque (_Macaca mulatta_) and the cynomolgus macaque (_Macaca fascicularis_) uncover polygenic signatures of reinforcement speciation"

# Abstract

Speciation can involve phases of divergent adaptation in allopatry and ecological/reproductive character displacement in sympatry or parapatry. Reproductive character displacement can result as a means of preventing hybridization, a process known as reinforcement speciation. In this study, we use whole-genome sequencing (WGS) of two closely related primate species that have experienced introgression in their history, the rhesus (Macaca mulatta) and cynomolgus (M. fascicularis) macaques, to identify genes exhibiting reproductive character displacement and other patterns consistent with reinforcement speciation. Using windowed scans of various population genetic statistics to identify signatures of reinforcement, we find 184 candidate genes associated with a variety of functions, including an overrepresentation of multiple neurological functions and several genes involved in sexual development and gametogenesis. These results are consistent with a variety of genes acting in a reinforcement process between these species. We also find signatures of introgression of the Y-chromosome that confirm previous studies suggesting male-driven introgression of _M. mulatta_ into _M. fascicularis_ populations. This study uses whole-genome sequencing to find evidence of the process of reinforcement in primates that have medical and conservation relevance.

# Subdirectories

The subdirectories are organized in order of actual usage (as best as possible) and scale of data processing (from largest to smallest); 1) [Genomics_Pipeline](https://github.com/StevisonLab/Reinforcement-Project/tree/main/Genomics_Pipeline); 2) [Genome_Analysis](https://github.com/StevisonLab/Reinforcement-Project/tree/main/Genome_Analysis); and 3) [Gene_Analysis](https://github.com/StevisonLab/Reinforcement-Project/tree/main/Gene_Analysis).

## Genomics Pipeline

The primary script here is Genomics_Pipeline.sh. This script will conduct a complete germline variant calling pipeline essentially as described in the book Genomics in the Cloud (Van der Auwera and O'Connor 2020) and constitutes a GATK Best Practices pipeline. This involves in brief; 1) downloading FASTQ files from NCBI SRA using fasterq-dump; 2) compressing these files using pigz; 3) aligning these files to a reference genome using bwa-mem; 4) compressing sam files to bam using samtools; 5) sorting bam files; 6) merging bam files common to a single invididual organism when necessary; 7) checking coverage of bam files; 8) checking other stats of bam files; 9) marking duplicates using GATK; 10) base quality score recalibration; 11) running checks on this output; 12) call variants; 13) create genomics database for joint genotyping; 14) run joint genotyping and retain invariant sites; 15) filter invariant sites with a hard-coded filter; 16) separately filter variant sites with variant quality score recalibration; 17) combine sites back together; 18) compress all bam files to cram; 19) use a publicly available script to check if cram has been generated properly.

**This script is heavily hard-coded.** We have run into numerous issues while running this analysis and inserted conditional statements for multiple possibilities and checks that different stages of the pipeline have already been run. However, several actual file names and directories and SLURM parameters are given as actually used for this study. The script can extract most relevant sample data from a standard NCBI SRA csv file, though there are sections hard-coded for extracting data from our Novogene sequenced files (though may be somewhat useful for Novogene sequences generally if modified). Some details of the script are currently set for calling X-chromosomes from female samples and comments throughout (and in the methods of the associated manuscript) should detail sections that may need to be changed for sex chromosomes. 

This directory also includes data generated using MultiQC and custom python scripts to compiling data from various file checks above (specifically, samtools coverage, samtools flagstat, and bcftools stats). These were used to generate relevant parts of Table 1 of the associated manuscript, and should be relatively flexible and user friendly as they have explicit help options and are relatively soft-coded. 

## Genome Analysis

These scripts primarily constitute different parts of a windowed genome-scale analysis of the samples generated from the above pipeline. PCA_Plot.R is the one exception to this as it conducts PCA over all autosomes joined together and X and Y chromosomes separately (and X chromosomes separate in females from males), and not in windows. It uses PLINK2 commands (commented out currently) directly in the same and then plots the output from these commands. It should probably be the first script used as it may provide basic information on samples to include or exclude from downstream windowed analyses (we excluded two samples based on our first PCA run).

Once beginning the window analysis, the first script to run is bcftools_MAF.sh, which is a short script necessary for generating the minor allele frequency data for our samples used by GenWin to generate genomic window intervals for analysis. The shell script is hard-coded to some degree but the essential bcftools components for a single file are like so:

```bcftools filter --threads 6 -i 'TYPE="snp" && QUAL!="."' all_sites.vcf.gz -O z -o variants_only.vcf.gz```

```bcftools +fill-tags variants_only.vcf.gz -- -t MAF > variants_only.MAF.vcf```

```gzip variants_only.MAF.vcf # optional compression step```

```bcftools query -f '%CHROM\t%POS\t%MAF\n' variants_only.MAF.vcf.gz > variants_only.MAF.txt```

GenWin_SNPs.R must then be run with this data (for all chromosomes) as input. It will output text files "*Spline.txt" of defined intervals. TajimaD.py, pixy_permutations.sh, and ADMIXTURE.sh run different summary statistics over the same data so all should use the VCF files and Spline.txt files for input and all don't need to be run in any particular order relative to each other. ADMIXTURE_parse.sh must be run after ADMIXTURE.sh as the output from the latter is confusing and doesn't directly connect data to samples. 

## Gene Analysis

These final scripts constitute some analysis of candidate genes extracted from the above data (process detailed more in associated manuscript). ParaCyno.Unique.py should be run first to isolate variants that are unique to parapatric (Para) _M. fascicularis_ (Cyno). This data can then be input into SnpEff like so:

```/tools/snpeff-5.0/scripts/snpEff eff -c /home/npb0015/Arctoides_unique/snpEff.config Mmul_10.99 -csvStats ParaCyno.Unique.stats.csv ParaCyno.Unique.vcf > ParaCyno.Unique.annotated.vcf```

The annotated.vcf file resulting from this can then be input into SnpSift (part of SnpEff, not the same as SIFT also used in this study!) to extract specific kinds of sites. For example:

```java -jar /tools/snpeff-5.0/SnpSift.jar filter "(ANN[*].EFFECT has 'upstream_gene_variant') | (ANN[*].EFFECT has 'intron_variant') | (ANN[*].EFFECT has 'downstream_gene_variant') | (ANN[*].EFFECT has '5_prime_UTR_variant') | (ANN[*].EFFECT has '3_prime_UTR_variant')" -f ../SnpEff_data/ParaCyno.Unique.annotated.vcf > ParaCyno.Unique.motif.vcf```

This produces variant that can be subject to motif analysis through tomtom as it includes a variant of types of sites (e.g. upstream variants or intron variants) that are possible binding sites for transcription factors. The above command specifically does not include coding variants, but an equivalent command including those can generate output usable by SIFT. 

So SIFT can be run on that data and output from that can be used for B-SIFT.sh. Similarly the output from the above command is necessary for running tomtom_parallel.sh. For that script, it is also necessary to generate reference and alternate FASTA files using something like bcftools consensus, which is also used in our previous study [Mitonuclear-Analysis-Project](https://github.com/StevisonLab/Mitonuclear-Analysis-Project). tomtom_extraction_parallel.sh must be run after this script to extract data in a neat format. This is analogous to ADMIXTURE_parse.sh above. 

We neglect the details of generating or obtaining databases for SnpEff, SIFT, and tomtom as these are partly described in the associated manuscript for this study (and the respective documentation of these programs) and are too detailed to include here. 
