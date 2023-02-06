setwd("~/Desktop/Reinforcement_Project/permutations/")

library(GenomicRanges)
library(plyranges)
library(rtracklayer)
library(dplyr)
require(data.table)

# Load in annotation files for rheMac10 (including gene and repeat element files)

GTF = import('../rheMac10.ensGene.gtf.gz')
Transcripts = subset(GTF, GTF$type == 'transcript')
EnsGenetoHGNC = fread('../ensemblToGeneName.txt', header = FALSE)
TE_file = import('../rheMac10.TE.bed.gz', extraCols = c(Name = 'character', Class = 'character', Family = 'character'))
#ensembl = useEnsembl(biomart="genes", dataset="mmulatta_gene_ensembl")

# Read in data files for ADMIXTURE and Tajima's D calculations and create genomic range objects for them

ADMIXTURE = fread('../ADMIXTURE/ParaCyno.admixture.txt')
ADMIXTURE_GRange = GRanges(seqname = ADMIXTURE$V1, IRanges(start = as.numeric(ADMIXTURE$V2) + 1, 
      end = as.numeric(ADMIXTURE$V3)), admixture = ADMIXTURE$V4)

TajimaD = fread('../TajimaD.Output.txt')

# Sanity-check that Tajima's D window lengths are positive values

TajimaD_windows = subset(TajimaD, V2 < V3)

# Subset populations

TajimaD_ParaCyno = subset(TajimaD_windows, V4 == 'ParaCyno')
TajimaD_AlloCyno = subset(TajimaD_windows, V4 == 'AlloCyno')
TajimaD_ParaRhe = subset(TajimaD_windows, V4 == 'ParaRhe')
TajimaD_AlloRhe = subset(TajimaD_windows, V4 == 'AlloRhe')

# Perform calculation to normalize ParaCyno Tajima's D values compared to neighboring populations

TajimaD_CynoDiff = TajimaD_ParaCyno$V5 - TajimaD_AlloCyno$V5
TajimaD_ParaDiff = TajimaD_ParaCyno$V5 - TajimaD_ParaRhe$V5
TajimaD_AvgDiff = (TajimaD_CynoDiff + TajimaD_ParaDiff)/2

TajimaD_all_GRange = GRanges(seqname = TajimaD_ParaCyno$V1, IRanges(start = as.numeric(TajimaD_ParaCyno$V2) + 1, 
          end = as.numeric(TajimaD_ParaCyno$V3)), ParaCyno = TajimaD_ParaCyno$V5, 
          AlloCyno = TajimaD_AlloCyno$V5, ParaRhe = TajimaD_ParaRhe$V5, 
          AlloRhe = TajimaD_AlloRhe$V5, TajimaDiff = TajimaD_AvgDiff)

# Run loop through all chromosomes to extract data for Dxy and Fst calculations

combined_stats = NULL
for (CHROM in c(1:20,'X','Y'))
{
# Extract data from empirical Dxy and Fst files
  
  result_dxy = paste("Reinforcement.SplineDefined.chr", CHROM, "_dxy.txt", sep = "")
  dxy_file = fread(result_dxy)
  dxy_file = subset(dxy_file, window_pos_1 < window_pos_2)
  result_fst = paste("Reinforcement.SplineDefined.chr", CHROM, "_fst.txt", sep = "")
  fst_file = fread(result_fst)
  fst_file = subset(fst_file, window_pos_1 < window_pos_2)

# Subset Dxy and Fst calculations for all pairwise comparisons
  
  ParaChromDxy = subset(dxy_file, pop1 == "ParaCyno" & pop2 == "ParaRhe")
  AlloChromDxy = subset(dxy_file, pop1 == "AlloCyno" & pop2 == "AlloRhe")
  CynoChromDxy = subset(dxy_file, pop1 == "ParaCyno" & pop2 == "AlloCyno")
  RheChromDxy = subset(dxy_file, pop1 == "AlloRhe" & pop2 == "ParaRhe")
  ParaChromFst = subset(fst_file, pop1 == "ParaCyno" & pop2 == "ParaRhe")
  AlloChromFst = subset(fst_file, pop1 == "AlloCyno" & pop2 == "AlloRhe")
  CynoChromFst = subset(fst_file, pop1 == "ParaCyno" & pop2 == "AlloCyno")
  RheChromFst = subset(fst_file, pop1 == "AlloRhe" & pop2 == "ParaRhe")
 
# Convert negative Fst values to 0
  
  ParaChromFst$avg_wc_fst[ParaChromFst$avg_wc_fst < 0] = 0 
  AlloChromFst$avg_wc_fst[AlloChromFst$avg_wc_fst < 0] = 0 
  CynoChromFst$avg_wc_fst[CynoChromFst$avg_wc_fst < 0] = 0 
  RheChromFst$avg_wc_fst[RheChromFst$avg_wc_fst < 0] = 0 

# Take the difference of Dxy and Fst between the parapatric and allopatric populations
# These represent genomic signature of character displacement
  
  DxyDiff = ParaChromDxy$avg_dxy - AlloChromDxy$avg_dxy
  FstDiff = ParaChromFst$avg_wc_fst - AlloChromFst$avg_wc_fst

# Create genomic range objects incorporating all the above calculations
# Also use subsetByOverlaps to obtain windows that have Dxy and Fst values
  
  fst_stats = data.frame(cbind(ParaChromFst$chromosome,ParaChromFst$window_pos_1,ParaChromFst$window_pos_2))
  dxy_stats_GRange = GRanges(seqname=ParaChromDxy$chromosome, IRanges(start=ParaChromDxy$window_pos_1 + 1, end=ParaChromDxy$window_pos_2), 
                             DxyDiff = DxyDiff, ParaDxy = ParaChromDxy$avg_dxy, 
                             AlloDxy = AlloChromDxy$avg_dxy, CynoDxy = CynoChromDxy$avg_dxy,
                             RheDxy = RheChromDxy$avg_dxy)
  fst_stats_GRange = GRanges(seqname=ParaChromFst$chromosome, IRanges(start=ParaChromFst$window_pos_1 + 1, end=ParaChromFst$window_pos_2))
  dxy_overlaps = subsetByOverlaps(dxy_stats_GRange, fst_stats_GRange)
  dxy_fst_combined = cbind(fst_stats,FstDiff,ParaChromFst$avg_wc_fst, 
                           AlloChromFst$avg_wc_fst, CynoChromFst$avg_wc_fst, RheChromFst$avg_wc_fst,  
                           dxy_overlaps$DxyDiff, dxy_overlaps$ParaDxy, 
                           dxy_overlaps$AlloDxy, dxy_overlaps$CynoDxy, dxy_overlaps$RheDxy)
  names(dxy_fst_combined) = c('CHROM','START','END','FstDiff_original',
                              'ParaFst_original','AlloFst_original','CynoFst_original','RheFst_original',
                              'DxyDiff_original','ParaDxy_original','AlloDxy_original','CynoDxy_original',
                              'RheDxy_original')

# Loop through 100 permutations to extract data for Dxt and Fst calculations
# This largely follows the same structure as above
  
for (permutation in 1:100)
{
  result_dxy = paste('Randomized.SplineDefined.chr', CHROM, '.', permutation, '_dxy.txt.gz', sep = '')
  dxy_file = fread(result_dxy)
  dxy_file = subset(dxy_file, window_pos_1 < window_pos_2)
  result_fst = paste('Randomized.SplineDefined.chr', CHROM, '.', permutation, '_fst.txt.gz', sep = '')
  fst_file = fread(result_fst)
  fst_file = subset(fst_file, window_pos_1 < window_pos_2)
  ParaChromDxy = subset(dxy_file, pop1 == 'ParaCyno' & pop2 == 'ParaRhe')
  AlloChromDxy = subset(dxy_file, pop1 == 'AlloCyno' & pop2 == 'AlloRhe')
  ParaChromFst = subset(fst_file, pop1 == 'ParaCyno' & pop2 == 'ParaRhe')
  AlloChromFst = subset(fst_file, pop1 == 'AlloCyno' & pop2 == 'AlloRhe')
  ParaChromFst$avg_wc_fst[ParaChromFst$avg_wc_fst < 0] = 0 
  AlloChromFst$avg_wc_fst[AlloChromFst$avg_wc_fst < 0] = 0 
  DxyDiff = ParaChromDxy$avg_dxy - AlloChromDxy$avg_dxy
  FstDiff = ParaChromFst$avg_wc_fst - AlloChromFst$avg_wc_fst
  dxy_stats = data.frame(cbind(ParaChromDxy$chromosome,ParaChromDxy$window_pos_1,ParaChromDxy$window_pos_2,DxyDiff))
  fst_stats = data.frame(cbind(ParaChromFst$chromosome,ParaChromFst$window_pos_1,ParaChromFst$window_pos_2,FstDiff))
  dxy_stats_GRange = GRanges(seqname=ParaChromDxy$chromosome, IRanges(start=ParaChromDxy$window_pos_1 + 1, end=ParaChromDxy$window_pos_2), data_id = DxyDiff)
  fst_stats_GRange = GRanges(seqname=ParaChromFst$chromosome, IRanges(start=ParaChromFst$window_pos_1 + 1, end=ParaChromFst$window_pos_2))
  dxy_overlaps = subsetByOverlaps(dxy_stats_GRange, fst_stats_GRange)
  dxy_fst_combined = cbind(dxy_fst_combined,FstDiff,dxy_overlaps$data_id)

# Create column names for all permutations in an automated fashion
  
  names(dxy_fst_combined)[permutation*2 + 12] = paste('FstDiff_', permutation, sep = '')
  names(dxy_fst_combined)[permutation*2 + 13] = paste('DxyDiff_', permutation, sep = '')
}
# Append empirical and permutated calculations for chromosome to a data frame
  
   combined_stats = rbind(combined_stats, dxy_fst_combined)
}

# Reorder the combined_stats file generated above to place Dxy/Fst permutations in proper order then create genomic ranges object

reordered_stats = combined_stats[ , c(1:13, seq(14,212,2), seq(15,213,2))]
reordered_stats_GRange = GRanges(seqname = reordered_stats$CHROM, IRanges(start = as.numeric(reordered_stats$START) + 1, 
                                end = as.numeric(reordered_stats$END)), FST = reordered_stats$FstDiff_original, DXY = reordered_stats$DxyDiff_original)

# Merge Tajima's D and ADMIXTURE ranges, then merge these both to reordered_stats
# This results in aligned data for all stats of interest

Tajima_ADMIXTURE_dframe = mergeByOverlaps(TajimaD_all_GRange, ADMIXTURE_GRange)
Tajima_Divergence_dframe = mergeByOverlaps(Tajima_ADMIXTURE_dframe$TajimaD_all_GRange, reordered_stats_GRange)
ADMIXTURE_Divergence_dframe = mergeByOverlaps(Tajima_ADMIXTURE_dframe$ADMIXTURE_GRange, reordered_stats_GRange)

# Loop through all chromosomes to make a vector of chromosomes aligned to stats
# The replicate function is using two parameters
# These are the "lengths" or numbers of times each chromosome is associated with a window in the file
# and the then this number is used to replicate the chromosome name or "value" that many times

all_stats_chromosomes = NULL
for (CHROM in c(1:22)) 
{
  all_stats_chromosomes = append(all_stats_chromosomes, 
  replicate(ADMIXTURE_Divergence_dframe$`Tajima_ADMIXTURE_dframe$ADMIXTURE_GRange`@seqnames@lengths[CHROM],
  toString(ADMIXTURE_Divergence_dframe$`Tajima_ADMIXTURE_dframe$ADMIXTURE_GRange`@seqnames@values[CHROM])))
}

# All the above data is combined to generate a bed file for all stats for all windows
# This data is then renamed and made into a genomic ranges object
# Then it is output into a file

all_stats_bed = data.frame(cbind(all_stats_chromosomes,
           ADMIXTURE_Divergence_dframe$`Tajima_ADMIXTURE_dframe$ADMIXTURE_GRange`@ranges@start,
           ADMIXTURE_Divergence_dframe$`Tajima_ADMIXTURE_dframe$ADMIXTURE_GRange`@ranges@start + ADMIXTURE_Divergence_dframe$`Tajima_ADMIXTURE_dframe$ADMIXTURE_GRange`@ranges@width - 1,
           ADMIXTURE_Divergence_dframe$admixture,
           ADMIXTURE_Divergence_dframe$FST, ADMIXTURE_Divergence_dframe$DXY, 
           Tajima_Divergence_dframe$ParaCyno, Tajima_Divergence_dframe$AlloCyno,
           Tajima_Divergence_dframe$ParaRhe, Tajima_Divergence_dframe$AlloRhe,
           Tajima_Divergence_dframe$TajimaDiff))
names(all_stats_bed) = c('CHROM', 'START', 'END', 'ParaCyno_ADMIXTURE', 
                         'FstDiff', 'DxyDiff', 'ParaCyno_D', 'AlloCyno_D', 
                         'ParaRhe_D', 'AlloRhe_D', 'TajimaDiff')

all_stats_GRange = GRanges(seqname = all_stats_chromosomes, 
                   IRanges(start = ADMIXTURE_Divergence_dframe$`Tajima_ADMIXTURE_dframe$ADMIXTURE_GRange`@ranges@start, 
                   end = ADMIXTURE_Divergence_dframe$`Tajima_ADMIXTURE_dframe$ADMIXTURE_GRange`@ranges@start + ADMIXTURE_Divergence_dframe$`Tajima_ADMIXTURE_dframe$ADMIXTURE_GRange`@ranges@width - 1),
                   ParaCyno_ADMIXTURE = ADMIXTURE_Divergence_dframe$admixture, Fst = ADMIXTURE_Divergence_dframe$FST,
                   Dxy = ADMIXTURE_Divergence_dframe$DXY,ParaCyno_D = Tajima_Divergence_dframe$ParaCyno, 
                   AlloCyno_D = Tajima_Divergence_dframe$AlloCyno, ParaRhe_D = Tajima_Divergence_dframe$ParaRhe,
                   AlloRhe_D = Tajima_Divergence_dframe$AlloRhe, TajimaDiff = Tajima_Divergence_dframe$TajimaDiff)

write.table(all_stats_bed, file = "../Genome_Stats.bed", quote = FALSE,
            row.names = FALSE, sep = '\t', col.names = TRUE)

# Merge rheMac10 transcripts annotation to genomic range above

transcripts_overlap_all = mergeByOverlaps(Transcripts, all_stats_GRange)

# Run function like above to get chromosomes for transcripts the proper amount of times

transcript_chromosomes_all = NULL
for (CHROM in c(1:22)) 
{
  transcript_chromosomes_all = append(transcript_chromosomes_all, 
  replicate(transcripts_overlap_all$all_stats_GRange@seqnames@lengths[CHROM], 
  toString(transcripts_overlap_all$all_stats_GRange@seqnames@values[CHROM])))
}

# Make objects with all statistics and Ensembl gene and transcript names

transcripts_bed_all = data.frame(cbind(transcript_chromosomes_all,
                  transcripts_overlap_all$all_stats_GRange@ranges@start,
                  transcripts_overlap_all$all_stats_GRange@ranges@start + transcripts_overlap_all$all_stats_GRange@ranges@width - 1,
                  transcripts_overlap_all$gene_id, transcripts_overlap_all$transcript_id,
                  transcripts_overlap_all$ParaCyno_ADMIXTURE,
                  transcripts_overlap_all$Fst, transcripts_overlap_all$Dxy, 
                  transcripts_overlap_all$ParaCyno_D, transcripts_overlap_all$AlloCyno_D,
                  transcripts_overlap_all$ParaRhe_D, transcripts_overlap_all$AlloRhe_D,
                  transcripts_overlap_all$TajimaDiff))
names(transcripts_bed_all) = c('CHROM', 'START', 'END', 'Gene',
                               'Transcript', 'ParaCyno_ADMIXTURE', 'FstDiff', 
                               'DxyDiff', 'ParaCyno_D', 'AlloCyno_D', 'ParaRhe_D', 
                               'AlloRhe_D', 'TajimaDiff')

# Do all the above with a repeat element annotation file instead of transcripts

transposons_overlap_all = mergeByOverlaps(TE_file, all_stats_GRange)

transposon_chromosomes_all = NULL
for (CHROM in c(1:22)) 
{
  transposon_chromosomes_all = append(transposon_chromosomes_all, 
  replicate(transposons_overlap_all$all_stats_GRange@seqnames@lengths[CHROM], 
  toString(transposons_overlap_all$all_stats_GRange@seqnames@values[CHROM])))
}

transposons_bed_all = data.frame(cbind(transposon_chromosomes_all,
                      transposons_overlap_all$all_stats_GRange@ranges@start,
                      transposons_overlap_all$all_stats_GRange@ranges@start + transposons_overlap_all$all_stats_GRange@ranges@width - 1,
                      transposons_overlap_all$Name, transposons_overlap_all$Class,transposons_overlap_all$Family,
                      transposons_overlap_all$ParaCyno_ADMIXTURE,
                      transposons_overlap_all$Fst, transposons_overlap_all$Dxy, 
                      transposons_overlap_all$ParaCyno_D, transposons_overlap_all$AlloCyno_D,
                      transposons_overlap_all$ParaRhe_D, transposons_overlap_all$AlloRhe_D,
                      transposons_overlap_all$TajimaDiff))
names(transposons_bed_all) = c('CHROM', 'START', 'END', 'Name','Class','Family', 'ParaCyno_ADMIXTURE', 'FstDiff', 'DxyDiff', 'ParaCyno_D', 'AlloCyno_D', 'ParaRhe_D', 'AlloRhe_D', 'TajimaDiff')

# Output bed files for transcripts and repeat elements separately

write.table(transcripts_bed_all, file = "../Transcripts_Stats.bed", quote = FALSE,
            row.names = FALSE, sep = '\t', col.names = TRUE)

write.table(transposons_bed_all, file = "../Transposons_Stats.bed", quote = FALSE,
            row.names = FALSE, sep = '\t', col.names = TRUE)


fst_ranking = reordered_stats[ , c(4, 10:109)]
dxy_ranking = reordered_stats[ , c(7, 110:209)]

Fst_P5 = NULL
for (row in 1:NROW(fst_ranking))
{
  if (rank(fst_ranking[row,], na.last = FALSE)[1] > 96)
    Fst_P5 = rbind(Fst_P5, reordered_stats[row, c(1:4)])
}

Fst_P5_GRange = GRanges(seqname = Fst_P5$CHROM, IRanges(start = as.numeric(Fst_P5$START) + 1, end = as.numeric(Fst_P5$END)), data_id = Fst_P5$FstDiff_original)

Dxy_P5 = NULL
for (row in 1:NROW(dxy_ranking))
{
  if (rank(dxy_ranking[row,], na.last = FALSE)[1] > 96 )
    Dxy_P5 = rbind(Dxy_P5, reordered_stats[row, c(1:3,7)])
  
}

Dxy_P5_GRange = GRanges(seqname=Dxy_P5$CHROM, IRanges(start = as.numeric(Dxy_P5$START) + 1, end = as.numeric(Dxy_P5$END)), data_id = Dxy_P5$DxyDiff_original)

TajimaD_bottom = subset(all_stats_GRange, as.numeric(TajimaDiff) < -1)

ADMIXTURE_2.5_percent = subset(all_stats_GRange, as.numeric(ParaCyno_ADMIXTURE) < 0.025) 

Top_Dxy_overlaps = mergeByOverlaps(all_stats_GRange, Dxy_P5_GRange)
Top_Fst_overlaps = mergeByOverlaps(all_stats_GRange, Fst_P5_GRange)
Bottom_TajimaD_overlaps = mergeByOverlaps(all_stats_GRange, TajimaD_bottom)
Bottom_ADMIXTURE_overlaps = mergeByOverlaps(all_stats_GRange, ADMIXTURE_2.5_percent)
TajimaD_ADMIXTURE_overlaps = mergeByOverlaps(Bottom_TajimaD_overlaps$TajimaD_bottom, Bottom_ADMIXTURE_overlaps$ADMIXTURE_2.5_percent)
Fst_TajimaD_overlap = mergeByOverlaps(Top_Fst_overlaps$Fst_P5_GRange, TajimaD_ADMIXTURE_overlaps$`Bottom_TajimaD_overlaps$TajimaD_bottom`)
Dxy_TajimaD_overlap = mergeByOverlaps(Top_Dxy_overlaps$Dxy_P5_GRange, TajimaD_ADMIXTURE_overlaps$`Bottom_TajimaD_overlaps$TajimaD_bottom`)
colnames(Fst_TajimaD_overlap)[colnames(Fst_TajimaD_overlap) == "TajimaD_ADMIXTURE_overlaps$`Bottom_TajimaD_overlaps$TajimaD_bottom`"] <- "final_data" 
colnames(Dxy_TajimaD_overlap)[colnames(Dxy_TajimaD_overlap) == "TajimaD_ADMIXTURE_overlaps$`Bottom_TajimaD_overlaps$TajimaD_bottom`"] <- "final_data" 
Fst_unique = subsetByOverlaps(Fst_TajimaD_overlap$final_data, Dxy_TajimaD_overlap$final_data, invert = TRUE)
Fst_unique$outlier_type = "Fst"
Dxy_unique = subsetByOverlaps(Dxy_TajimaD_overlap$final_data, Fst_TajimaD_overlap$final_data, invert = TRUE)
Dxy_unique$outlier_type = "Dxy"
Fst_Dxy_overlap = mergeByOverlaps(Fst_TajimaD_overlap$final_data, Dxy_TajimaD_overlap$final_data)
final_stats = c(Fst_unique, Dxy_unique, Fst_Dxy_overlap$`Dxy_TajimaD_overlap$final_data`)
transcripts_overlap_outliers = mergeByOverlaps(Transcripts, final_stats)
transposons_overlap_outliers = mergeByOverlaps(TE_file, final_stats)
length(unique(transcripts_overlap_outliers$Transcripts@seqnames))

transcript_chromosomes_outliers = NULL
for (CHROM in c(1:20)) 
{
  transcript_chromosomes_outliers = append(transcript_chromosomes_outliers, 
  replicate(transcripts_overlap_outliers$final_stats@seqnames@lengths[CHROM], 
  toString(transcripts_overlap_outliers$final_stats@seqnames@values[CHROM])))
}

transposon_chromosomes_outliers = NULL
for (CHROM in c(1:20)) 
{
  transposon_chromosomes_outliers = append(transposon_chromosomes_outliers, 
  replicate(transposons_overlap_outliers$final_stats@seqnames@lengths[CHROM], 
  toString(transposons_overlap_outliers$final_stats@seqnames@values[CHROM])))
}

transcripts_bed_outliers = data.frame(cbind(transcript_chromosomes_outliers,
           transcripts_overlap_outliers$final_stats@ranges@start,
           transcripts_overlap_outliers$final_stats@ranges@start + transcripts_overlap_outliers$final_stats@ranges@width - 1,
           transcripts_overlap_outliers$gene_id, transcripts_overlap_outliers$transcript_id,
           transcripts_overlap_outliers$ParaCyno_ADMIXTURE,
           transcripts_overlap_outliers$Fst, transcripts_overlap_outliers$Dxy, 
           transcripts_overlap_outliers$ParaCyno_D, transcripts_overlap_outliers$AlloCyno_D,
           transcripts_overlap_outliers$ParaRhe_D, transcripts_overlap_outliers$AlloRhe_D, 
           transcripts_overlap_outliers$TajimaDiff, transcripts_overlap_outliers$outlier_type))
names(transcripts_bed_outliers) = c('CHROM', 'START', 'END', 'Gene',
                                    'Transcript', 'ParaCyno_ADMIXTURE', 'FstDiff', 
                                    'DxyDiff', 'ParaCyno_D', 'AlloCyno_D', 'ParaRhe_D', 
                                    'AlloRhe_D', 'TajimaDiff', 'Outlier_type')

ADMIXTURE_ranked = rank(transcripts_bed_outliers$ParaCyno_ADMIXTURE, ties.method = "min")
Fst_ranked = rank(transcripts_bed_outliers$FstDiff, ties.method = "max")
Dxy_ranked = rank(transcripts_bed_outliers$DxyDiff, ties.method = "max")
Tajima_ranked = rank(transcripts_bed_outliers$TajimaDiff, ties.method = "min")

transcripts_bed_ranked = data.frame(cbind(transcript_chromosomes_outliers,
                           transcripts_overlap_outliers$final_stats@ranges@start,
                           transcripts_overlap_outliers$final_stats@ranges@start + transcripts_overlap_outliers$final_stats@ranges@width - 1,
                           transcripts_overlap_outliers$gene_id, transcripts_overlap_outliers$transcript_id,
                           ADMIXTURE_ranked, Fst_ranked, Dxy_ranked,Tajima_ranked, ADMIXTURE_ranked + Fst_ranked + Dxy_ranked + Tajima_ranked,
                           transcripts_overlap_outliers$ParaCyno_ADMIXTURE,
                           transcripts_overlap_outliers$Fst, transcripts_overlap_outliers$Dxy, 
                           transcripts_overlap_outliers$ParaCyno_D, transcripts_overlap_outliers$AlloCyno_D, 
                           transcripts_overlap_outliers$ParaRhe_D, transcripts_overlap_outliers$AlloRhe_D, 
                           transcripts_overlap_outliers$TajimaDiff, transcripts_overlap_outliers$outlier_type))
names(transcripts_bed_ranked) = c('CHROM', 'START', 'END', 'Gene',
                                  'Transcript', 'ADMIXTURE_rank', 'Fst_rank',
                                  'Dxy_rank', 'Tajima_rank', 'Rank','ParaCyno_ADMIXTURE', 
                                  'FstDiff', 'DxyDiff', 'ParaCyno_D', 'AlloCyno_D', 
                                  'ParaRhe_D', 'AlloRhe_D', 'TajimaDiff', 'Outlier_type')

transposons_bed_outliers = data.frame(cbind(transposon_chromosomes_outliers,
           transposons_overlap_outliers$final_stats@ranges@start,
           transposons_overlap_outliers$final_stats@ranges@start + transposons_overlap_outliers$final_stats@ranges@width - 1,
           transposons_overlap_outliers$Name, transposons_overlap_outliers$Class, transposons_overlap_outliers$Family,
           transposons_overlap_outliers$ParaCyno_ADMIXTURE,
           transposons_overlap_outliers$Fst, transposons_overlap_outliers$Dxy, 
           transposons_overlap_outliers$ParaCyno_D, transposons_overlap_outliers$AlloCyno_D,
           transposons_overlap_outliers$ParaRhe_D, transposons_overlap_outliers$AlloRhe_D, 
           transposons_overlap_outliers$TajimaDiff, transposons_overlap_outliers$outlier_type))
names(transposons_bed_outliers) = c('CHROM', 'START', 'END', 'Name',
                                    'Class','Family', 'ParaCyno_ADMIXTURE', 
                                    'FstDiff', 'DxyDiff', 'ParaCyno_D', 'AlloCyno_D', 
                                    'ParaRhe_D', 'AlloRhe_D', 'TajimaDiff', 'Outlier_type')

write.table(transcripts_bed_outliers, file="../Outlier_Transcripts.bed", quote=FALSE,
            row.names = FALSE, sep = '\t', col.names = TRUE)

write.table(transcripts_bed_ranked, file="../Outlier_Transcripts_Ranked.bed", quote=FALSE,
            row.names = FALSE, sep = '\t', col.names = TRUE)

write.table(transposons_bed_outliers, file="../Outlier_Transposons.bed", quote=FALSE,
            row.names = FALSE, sep = '\t', col.names = TRUE)

write.table(transcripts_bed_outliers$Transcript, file="../Outlier_Transcripts.txt", quote=FALSE,
            row.names = FALSE, sep = '\t', col.names = TRUE)

write.table(transcripts_bed_outliers$Gene, file="../Outlier_Genes.txt", quote=FALSE,
            row.names = FALSE, sep = '\t', col.names = TRUE)


png("../Fst_Histogram.png", res=700, height = 5.29, width = 8, units="in")
plot(ParaFstHist, col = rgb(1,0,0,0.4), xlab = 'Genetic Divergence (Fst)', main = NA)
plot(AlloFstHist, col = rgb(0,0,1,0.4), add = TRUE)
legend(0.5, 50000, c(paste('Parapatric Fst: Mean = ', 
      signif(mean(combined_stats$ParaFst_original, na.rm = TRUE), digits = 2)),
        paste('Allopatric Fst: Mean = ', signif(mean(combined_stats$AlloFst_original, na.rm = TRUE), digits = 2))), 
       fill = c(rgb(1,0,0,0.4), rgb(0,0,1,0.4)), cex = 1, bty = 'n')
text(0.8, 25000, "SNP Count = 31,639,415\nSite Count = 1,333,185,598\nSNP/Sites = 2.37%", cex = 1)
dev.off()

plot(ParaDxyHist, col = rgb(1,0,0,0.4), xlab = 'Genetic Divergence (Fst)', main = NA)
plot(AlloDxyHist, col = rgb(0,0,1,0.4), add = TRUE)
legend(0.5, 50000, c(paste('Parapatric Dxy: Mean = ', 
                           signif(mean(combined_stats$ParaDxy_original, na.rm = TRUE), digits = 2)),
                     paste('Allopatric Dxy: Mean = ', signif(mean(combined_stats$AlloDxy_original, na.rm = TRUE), digits = 2))), 
       fill = c(rgb(1,0,0,0.4), rgb(0,0,1,0.4)), cex = 1, bty = 'n')
text(0.8, 25000, "SNP Count = 31,639,415\nSite Count = 1,333,185,598\nSNP/Sites = 2.37%", cex = 1)


