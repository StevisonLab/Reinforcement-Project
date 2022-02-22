setwd("~/Desktop/Reinforcement_Project/Spline_Analysis/")

library(dplyr)
require(data.table)
combined_stats = NULL
for (CHROM in 1:20)
{
  result_dxy = paste("Reinforcement.SplineDefined.chr", CHROM, "_dxy.txt", sep = "")
  dxy_file = fread(result_dxy)
  result_fst = paste("Reinforcement.SplineDefined.chr", CHROM, "_fst.txt", sep = "")
  fst_file = fread(result_fst)
  ParaChromDxy = subset(dxy_file, pop1 == "ParaCyno" & pop2 == "ParaRhe")
  AlloChromDxy = subset(dxy_file, pop1 == "AlloCyno" & pop2 == "AlloRhe")
  ParaChromFst = subset(fst_file, pop1 == "ParaCyno" & pop2 == "ParaRhe")
  AlloChromFst = subset(fst_file, pop1 == "AlloCyno" & pop2 == "AlloRhe")
  ParaChromFst$avg_wc_fst[ParaChromFst$avg_wc_fst < 0] = 0 
  AlloChromFst$avg_wc_fst[AlloChromFst$avg_wc_fst < 0] = 0 
  DxyDiff = ParaChromDxy$avg_dxy - AlloChromDxy$avg_dxy
  FstDiff = ParaChromFst$avg_wc_fst - AlloChromFst$avg_wc_fst
  stat_name = paste("../VCFS/filtered_vcfs/SNP_vcfs/Reinforcement.chr", CHROM, ".MAF.Spline.txt", sep = "")
  stat_file = fread(stat_name)
  res_stats = data.frame(cbind(ParaChromDxy$chromosome,AlloChromFst$window_pos_1,ParaChromFst$window_pos_2,DxyDiff,FstDiff,stat_file$Wstat))
  res_stats = subset(res_stats, !is.na(res_stats$V6))
  combined_stats = rbind(combined_stats, res_stats)
}
write.table(combined_stats, file="Reinforcement.ParaAlloDiff.bed", quote=FALSE,
            row.names = FALSE, sep = '\t', col.names = FALSE)

  Dxy_ordered_stats = combined_stats[order((as.numeric(combined_stats$DxyDiff)), decreasing=TRUE),]
  Fst_ordered_stats = combined_stats[order((as.numeric(combined_stats$FstDiff)), decreasing=TRUE),]
  Dxy_top = slice_head(Dxy_ordered_stats, prop = 0.01)
  Fst_top = slice_head(Fst_ordered_stats, prop = 0.01)
  write.table(Dxy_top, file="Reinforcement.top1percent.Dxy.bed", quote=FALSE,
              row.names = FALSE, sep = '\t', col.names = FALSE)
  write.table(Fst_top, file="Reinforcement.top1percent.Fst.bed", quote=FALSE,
              row.names = FALSE, sep = '\t', col.names = FALSE)
  
