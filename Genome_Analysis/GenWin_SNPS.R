library('GenWin')
setwd('~/Desktop/Reinforcement_Project/VCFS/filtered_vcfs/SNP_vcfs')
require(data.table)

for (CHROM in 1:20)
{
  Filename = paste('Reinforcement.chr', CHROM, '.MAF.txt.gz', sep = "")
  File = fread(Filename)
  Spline = paste('Reinforcement.chr', CHROM, '.MAF.Spline.txt', sep = "")
  SplineList <- assign(Spline, splineAnalyze(Y = File$V3, smoothness = 100,
                       map = File$V2, plotRaw = FALSE,
                       plotWindows = FALSE, method = 4))
  write.table(SplineList$windowData, file=Spline, quote=FALSE,
              row.names = FALSE, sep = '\t')
}
