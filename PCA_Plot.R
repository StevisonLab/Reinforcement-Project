setwd('~/Desktop/Reinforcement_Project/pca')

library(tidyverse)
library(ggpubr)

# Run plink2 commands directly on command line to conduct PCA
# File names should be changed for autosome, chrY and chrX

#system("plink2 --threads 6 --keep Samples.txt --vcf Reinforcement.chr21.recombined.allsites.ID.vcf.gz --make-pgen --out Reinforcement.chr21.recombined.allsites.ID")
#system("mv Reinforcement.chr21.recombined.allsites.ID.log Reinforcement.chr21.recombined.allsites.ID.pfiles.log")
#system("echo \"#IID\tSID\tPAT\tMAT\tSEX\" > Reinforcement.chr21.recombined.allsites.ID.psam")
#system("paste Sample_List.plink.txt Populations.plink.txt | awk '{OFS=\"\t\" ; print $1,$4,0,0,$2}' >> Reinforcement.chr21.recombined.allsites.ID.psam")
#system("sed -i .backup \"s/M$/1/\" Reinforcement.chr21.recombined.allsites.ID.psam")
#system("sed -i .backup \"s/F$/2/\" Reinforcement.chr21.recombined.allsites.ID.psam")
#system("plink2 --threads 6 --psam Reinforcement.chr21.recombined.allsites.ID.psam --pvar Reinforcement.chr21.recombined.allsites.ID.pvar --pgen Reinforcement.chr21.recombined.allsites.ID.pgen --freq --out Reinforcement.chr21.recombined.allsites.ID")
#system("mv Reinforcement.chr21.recombined.allsites.ID.log Reinforcement.chr21.recombined.allsites.ID.afreq.log")
#system("plink2 --threads 6 --psam Reinforcement.chr21.recombined.allsites.ID.psam --pvar Reinforcement.chr21.recombined.allsites.ID.pvar --pgen Reinforcement.chr21.recombined.allsites.ID.pgen --read-freq Reinforcement.chr21.recombined.allsites.ID.afreq --pca --out Reinforcement.chr21.recombined.allsites.ID")
#system("mv Reinforcement.chr21.recombined.allsites.ID.log Reinforcement.chr21.recombined.allsites.ID.pca.log")

# Read in autosome and chrY eigenfiles respectively

AutoEigenValues <- read_delim("Reinforcement.autosomes.MAF.ID.eigenval", delim = "\t")
AutoEigenVectors <- read_delim("Reinforcement.autosomes.MAF.ID.eigenvec", delim = "\t")
AutoEigenPercent <- round((AutoEigenValues / (sum(AutoEigenValues))*100), 2)

YEigenValues <- read_delim("Reinforcement.chr22.filtered.variant.ID.eigenval", delim = "\t")
YEigenVectors <- read_delim("Reinforcement.chr22.filtered.variant.ID.eigenvec", delim = "\t")
YEigenPercent <- round((YEigenValues / (sum(YEigenValues))*100), 2)

# Create plots for autosome and chrY PCA respectively
AutoPlot = ggplot(data = AutoEigenVectors) +
  geom_point(mapping = aes(x = PC1, y = PC2, color = SID, shape = SID), 
             size = 3, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "Autosomes", x = paste0("PC1 (",AutoEigenPercent[1,1]," %)"),
       y = paste0("PC2 (",AutoEigenPercent[2,1]," %)"),
       colour = "Populations", shape = "Populations") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

YPlot = ggplot(data = YEigenVectors) +
  geom_point(mapping = aes(x = PC1, y = PC2, color = SID, shape = SID), 
             size = 3, show.legend = TRUE) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "Y-chromosome", x = paste0("PC1 (",YEigenPercent[1,1]," %)"),
       y = paste0("PC2 (",YEigenPercent[2,1]," %)"),
       colour = "Populations", shape = "Populations") +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.85), 
        legend.box.background = element_rect(color = "black", size = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

# Arrange both plots into a two panel image with some space in the middle

PCA_plot = ggarrange(AutoPlot, NULL, YPlot,
                    labels = c('B','', 'C'), font.label = list(size = 25),
                    ncol = 3, widths = c(1, 0.05, 1))

# Save image as PNG file

ggsave("PCA.png",
  PCA_plot, units = 'mm',
  width = 300,
  height = 150,
  bg = 'white',
  dpi = 800)
