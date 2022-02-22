setwd('~/Desktop/Reinforcement_Project')

library(tidyverse)

system("plink2 --vcf Reinforcement.chr21.recombined.allsites.ID.vcf.gz --out Reinforcement.chr21.recombined.allsites.ID")
system("mv Reinforcement.chr21.recombined.allsites.ID.log Reinforcement.chr21.recombined.allsites.ID.pfiles.log")
system("echo \"#IID\tSID\tPAT\tMAT\tSEX\" > Reinforcement.chr21.recombined.allsites.ID.psam")
system("paste Sample_List.txt VCFS/filtered_vcfs/sorted.Population_File_BCFtools.txt | awk '{OFS=\"\t\" ; print $1,$5,0,0,$2}' >> Reinforcement.chr21.recombined.allsites.ID.psam")
system("sed -i .backup \"s/M$/1/\" Reinforcement.chr21.recombined.allsites.ID.psam")
system("sed -i .backup \"s/F$/2/\" Reinforcement.chr21.recombined.allsites.ID.psam")
system("plink2 --psam Reinforcement.chr21.recombined.allsites.ID.psam --pvar Reinforcement.chr21.recombined.allsites.ID.pvar --pgen Reinforcement.chr21.recombined.allsites.ID.pgen --freq --out Reinforcement.chr21.recombined.allsites.ID")
system("mv Reinforcement.chr21.recombined.allsites.ID.log Reinforcement.chr21.recombined.allsites.ID.afreq.log")
system("plink2 --psam Reinforcement.chr21.recombined.allsites.ID.psam --pvar Reinforcement.chr21.recombined.allsites.ID.pvar --pgen Reinforcement.chr21.recombined.allsites.ID.pgen --read-freq Reinforcement.chr21.recombined.allsites.ID.afreq --pca --out Reinforcement.chr21.recombined.allsites.ID")

eigenValues <- read_delim("Reinforcement.chr21.recombined.allsites.ID.eigenval", delim = "\t")
eigenVectors <- read_delim("Reinforcement.chr21.recombined.allsites.ID.eigenvec", delim = "\t")

## Proportion of variation captured by each vector
eigen_percent <- round((eigenValues / (sum(eigenValues))*100), 2)

# PCA plot
ggplot(data = eigenVectors) +
  geom_point(mapping = aes(x = PC1, y = PC2, color = `#IID`, shape = SID), size = 3, show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "PCA of X-chromosomes",
       x = paste0("Principal component 1 (",eigen_percent[1,1]," %)"),
       y = paste0("Principal component 3 (",eigen_percent[3,1]," %)"),
       colour = "Samples", shape = "Populations") +
  theme_minimal()

