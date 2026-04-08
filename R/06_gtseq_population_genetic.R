rm(list = ls())

#########################
####### LIBRARY #########
#########################

library(tidyverse)
library(magrittr)
library(adegenet)
library(paletteer)
library(graph4lg)
source("99_genind2structure.R")
library(dartR)
library(readxl)

#########################
######### DATA ##########
#########################

#Read in genepop file including source individuals only
coar_all <- read.genepop("data/gtscore_analysis/polyGenResults_singleSNP_genepop_source.gen")

#Convert to ADMIXTURE format

#genind to genlight 
genlight_ref <- dartR::gi2gl(coar_all) 

#genlight to PLINK 
#assign SNP positions
pos_names <- genlight_ref@loc.names
pos_names <- str_remove(pos_names,"Chr\\d+_")
pos_names <- str_remove(pos_names,"_1")
chr_names <- str_extract(genlight_ref@loc.names,"Chr\\d+")
chr_names <- str_remove(chr_names,"Chr0")
chr_names <- str_remove(chr_names,"Chr")
genlight_ref@position <- as.factor(pos_names)
genlight_ref@chromosome <- as.factor(rep(1,410))

dartR :: gl2plink(genlight_ref, bed_file = T, plink_path = "~/Documents/10_Programmes/plink-1.07-mac-intel", outfile="reference_plink_gtseq",outpath=".")


#########################
######### PCA1 ##########
#########################


#In this part of the code I evaluate if there are weird samples
#in the dataset. Plot-twist, there are some I had to remove.

#Transform NA in the dataset
X <- tab(coar_all,NA.method = "mean")

#Make pca, with 20 axis
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE,nf=20)

#Evaluate % of variance explained for each axis
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

## 3 first axis explain more variation than expected by chance

#Transform pca in dataframe to plot with ggplot
pca <- as.data.frame(pca1$li)

#add column id to then join with the file which identifies population name 
#of each individual
pca$sample <- rownames(pca)

# Extract and reformat matches any letters ending in 's'
pca$sample <- stringr::str_replace(pca$sample, 
                       "([A-Za-z]*s)_(\\d+)_", 
                       function(x) {
                         prefix <- str_extract(x, "[A-Za-z]*s(?=_)")
                         num <- str_extract(x, "(?<=[A-Za-z]s_)\\d+")
                         paste0(prefix, "_", sprintf("%03d", as.numeric(num)), "_")
                       })

#import pop info
pop <- read.csv("data/JB_samples_metadata.csv",header = T)

pop %<>% mutate(sample_id = str_replace(sample_id,"-","_"))

pcapop <- pca %>% left_join(pop,by=c("sample"="sample_id"))


#import cluster info
clus <- read.table("data/baseline_assignment_gtseq.txt",header=T)

pcapop %<>% left_join(clus, by = c("sample" = "indiv"))

#Plot the PCA
ggplot(pcapop %>% filter(!is.na(collection)),aes(x = Axis1, y = Axis2, color = collection))+
  geom_point()+
  scale_color_manual(values = c("#8B8989","#E69F00","#0072B2","skyblue1"))+
  xlab("PC1 (19.67 %)")+
  ylab("PC2 (5.66 %)")+
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

ggplot(pcapop,aes(x = Axis3, y = Axis4,color = collection))+
  geom_point()+
  scale_color_manual(values = c("#8B8989","#E69F00","#0072B2","skyblue1"))+
  xlab("PC3 (1.43 %)")+
  ylab("PC4 (1.21 %)")+
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))


#---------------------------------------------------------------#
#Check allele frequency difference for a specific loci


genlight_ref

pop(genlight_ref) <- pcapop$pop

genlight_ref@pop <- as.factor(pcapop$pop)


freq_by_pop <- glMean(genlight_ref, alleleAsUnit = FALSE)

