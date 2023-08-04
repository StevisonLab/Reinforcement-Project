setwd('~/Desktop/Reinforcement_Project/pca')

library(tidyverse)
library(ggpubr)

# Run plink2 commands directly on command line to conduct PCA for autosomes, chrY, chrX for males, and chrX for females

#system("plink2 --threads 6 --keep Samples.txt --vcf Reinforcement.autosomes.MAF.ID.vcf.gz --make-pgen --out Reinforcement.autosomes.MAF.ID")
#system("mv Reinforcement.autosomes.MAF.ID.log Reinforcement.autosomes.MAF.ID.pfiles.log")
#system("echo \"#IID\tSID\tPAT\tMAT\tSEX\" > Reinforcement.autosomes.MAF.ID.psam")
#system("paste Sample_List.plink.txt Populations.plink.txt | awk '{OFS=\"\t\" ; print $1,$4,0,0,$2}' >> Reinforcement.autosomes.MAF.ID.psam")
#system("sed -i .backup \"s/M$/1/\" Reinforcement.autosomes.MAF.ID.psam")
#system("sed -i .backup \"s/F$/2/\" Reinforcement.autosomes.MAF.ID.psam")
#system("plink2 --threads 6 --psam Reinforcement.autosomes.MAF.ID.psam --pvar Reinforcement.autosomes.MAF.ID.pvar --pgen Reinforcement.autosomes.MAF.ID.pgen --freq --out Reinforcement.autosomes.MAF.ID")
#system("mv Reinforcement.autosomes.MAF.ID.log Reinforcement.autosomes.MAF.ID.afreq.log")
#system("plink2 --threads 6 --psam Reinforcement.autosomes.MAF.ID.psam --pvar Reinforcement.autosomes.MAF.ID.pvar --pgen Reinforcement.autosomes.MAF.ID.pgen --read-freq Reinforcement.autosomes.MAF.ID.afreq --pca 37 --out Reinforcement.autosomes.MAF.ID")
#system("mv Reinforcement.autosomes.MAF.ID.log Reinforcement.autosomes.MAF.ID.pca.log")

#system("plink2 --threads 6 --keep Males.txt --vcf Reinforcement.chr22.filtered.variant.ID.vcf.gz --make-pgen --out Reinforcement.chr22.filtered.variant.ID")
#system("mv Reinforcement.chr22.filtered.variant.ID.log Reinforcement.chr22.filtered.variant.ID.pfiles.log")
#system("echo \"#IID\tSID\tPAT\tMAT\tSEX\" > Reinforcement.chr22.filtered.variant.ID.psam")
#system("paste Male_List.txt Populations.male.plink.txt | awk '{OFS=\"\t\" ; print $1,$4,0,0,$2}' >> Reinforcement.chr22.filtered.variant.ID.psam")
#system("sed -i .backup \"s/M$/1/\" Reinforcement.chr22.filtered.variant.ID.psam")
#system("sed -i .backup \"s/F$/2/\" Reinforcement.chr22.filtered.variant.ID.psam")
#system("plink2 --threads 6 --psam Reinforcement.chr22.filtered.variant.ID.psam --pvar Reinforcement.chr22.filtered.variant.ID.pvar --pgen Reinforcement.chr22.filtered.variant.ID.pgen --freq --out Reinforcement.chr22.filtered.variant.ID")
#system("mv Reinforcement.chr22.filtered.variant.ID.log Reinforcement.chr22.filtered.variant.ID.afreq.log")
#system("plink2 --threads 6 --psam Reinforcement.chr22.filtered.variant.ID.psam --pvar Reinforcement.chr22.filtered.variant.ID.pvar --pgen Reinforcement.chr22.filtered.variant.ID.pgen --read-freq Reinforcement.chr22.filtered.variant.ID.afreq --pca 22 --out Reinforcement.chr22.filtered.variant.ID")
#system("mv Reinforcement.chr22.filtered.variant.ID.log Reinforcement.chr22.filtered.variant.ID.pca.log")

#system("plink2 --threads 6 --keep Males.txt --vcf Reinforcement.chr21.variants.ID.vcf.gz --make-pgen --out Reinforcement.chr21.male")
#system("mv Reinforcement.chr21.male.log Reinforcement.chr21.male.pfiles.log")
#system("echo \"#IID\tSID\tPAT\tMAT\tSEX\" > Reinforcement.chr21.male.psam")
#system("paste Male_List.txt Populations.male.plink.txt | awk '{OFS=\"\t\" ; print $1,$4,0,0,$2}' >> Reinforcement.chr21.male.psam")
#system("sed -i .backup \"s/M$/1/\" Reinforcement.chr21.male.psam")
#system("sed -i .backup \"s/F$/2/\" Reinforcement.chr21.male.psam")
#system("plink2 --threads 6 --psam Reinforcement.chr21.male.psam --pvar Reinforcement.chr21.male.pvar --pgen Reinforcement.chr21.male.pgen --freq --out Reinforcement.chr21.male")
#system("mv Reinforcement.chr21.male.log Reinforcement.chr21.male.afreq.log")
#system("plink2 --threads 6 --psam Reinforcement.chr21.male.psam --pvar Reinforcement.chr21.male.pvar --pgen Reinforcement.chr21.male.pgen --read-freq Reinforcement.chr21.male.afreq --pca 22 --out Reinforcement.chr21.male")
#system("mv Reinforcement.chr21.male.log Reinforcement.chr21.male.pca.log")

#system("plink2 --threads 6 --keep Females.txt --vcf Reinforcement.chr21.variants.ID.vcf.gz --make-pgen --out Reinforcement.chr21.female")
#system("mv Reinforcement.chr21.female.log Reinforcement.chr21.female.pfiles.log")
#system("echo \"#IID\tSID\tPAT\tMAT\tSEX\" > Reinforcement.chr21.female.psam")
#system("paste Female_List.txt Populations.female.plink.txt | awk '{OFS=\"\t\" ; print $1,$4,0,0,$2}' >> Reinforcement.chr21.female.psam")
#system("sed -i .backup \"s/M$/1/\" Reinforcement.chr21.female.psam")
#system("sed -i .backup \"s/F$/2/\" Reinforcement.chr21.female.psam")
#system("plink2 --threads 6 --psam Reinforcement.chr21.female.psam --pvar Reinforcement.chr21.female.pvar --pgen Reinforcement.chr21.female.pgen --freq --out Reinforcement.chr21.female")
#system("mv Reinforcement.chr21.female.log Reinforcement.chr21.female.afreq.log")
#system("plink2 --threads 6 --psam Reinforcement.chr21.female.psam --pvar Reinforcement.chr21.female.pvar --pgen Reinforcement.chr21.female.pgen --read-freq Reinforcement.chr21.female.afreq --pca 14 --out Reinforcement.chr21.female")
#system("mv Reinforcement.chr21.female.log Reinforcement.chr21.female.pca.log")

# Read in autosome, chrY, chrX for males, and chrX for females eigenfiles respectively

AutoEigenValues <- read_delim("Reinforcement.autosomes.MAF.ID.eigenval", delim = "\t", col_names = FALSE)
AutoEigenVectors <- read_delim("Reinforcement.autosomes.MAF.ID.eigenvec", delim = "\t")
AutoEigenPercent <- round((AutoEigenValues / (sum(AutoEigenValues))*100), 2)

YEigenValues <- read_delim("Reinforcement.chr22.filtered.variant.ID.eigenval", delim = "\t", col_names = FALSE)
YEigenVectors <- read_delim("Reinforcement.chr22.filtered.variant.ID.eigenvec", delim = "\t")
YEigenPercent <- round((YEigenValues / (sum(YEigenValues))*100), 2)

XMEigenValues <- read_delim("Reinforcement.chr21.male.eigenval", delim = "\t", col_names = FALSE)
XMEigenVectors <- read_delim("Reinforcement.chr21.male.eigenvec", delim = "\t")
XMEigenPercent <- round((XMEigenValues / (sum(XMEigenValues))*100), 2)

XFEigenValues <- read_delim("Reinforcement.chr21.female.eigenval", delim = "\t", col_names = FALSE)
XFEigenVectors <- read_delim("Reinforcement.chr21.female.eigenvec", delim = "\t")
XFEigenPercent <- round((XFEigenValues / (sum(XFEigenValues))*100), 2)

# Create plots for autosome, chrY, male chrX, and female chrX PCA respectively

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

XMPlot = ggplot(data = XMEigenVectors) +
  geom_point(mapping = aes(x = PC1, y = PC2, color = SID, shape = SID), 
             size = 3, show.legend = TRUE) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "X-chromosome Males", x = paste0("PC1 (",XMEigenPercent[1,1]," %)"),
       y = paste0("PC2 (",XMEigenPercent[2,1]," %)"),
       colour = "Populations", shape = "Populations") +
  theme_minimal() +
  theme(legend.position = c(0.15, 0.85), 
        legend.box.background = element_rect(color = "black", size = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

XFPlot = ggplot(data = XFEigenVectors) +
  geom_point(mapping = aes(x = PC1, y = PC2, color = SID, shape = SID), 
             size = 3, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "X-chromosome Females", x = paste0("PC1 (",XFEigenPercent[1,1]," %)"),
       y = paste0("PC2 (",XFEigenPercent[2,1]," %)"),
       colour = "Populations", shape = "Populations") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

# Arrange autosome and chrY plots into a two panel image with some space in the middle to make main figure
# Do the same for male and female chrX to make supplemental figure

PCA_plot = ggarrange(AutoPlot, NULL, YPlot,
                    labels = c('B','', 'C'), font.label = list(size = 25),
                    ncol = 3, widths = c(1, 0.05, 1))

Supp_PCA_plot = ggarrange(XMPlot, NULL, XFPlot,
                    labels = c('A','', 'B'), font.label = list(size = 25),
                    ncol = 3, widths = c(1, 0.05, 1))

# Save images as PNG files

ggsave("PCA.png",
  PCA_plot, units = 'mm',
  width = 300,
  height = 150,
  bg = 'white',
  dpi = 800)

ggsave("Figure1_SuppInfo.png",
       Supp_PCA_plot, units = 'mm',
       width = 300,
       height = 150,
       bg = 'white',
       dpi = 800)
