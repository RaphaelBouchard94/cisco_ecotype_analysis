rm(list = ls())

#####################
#######Library#######
#####################

library(tidyverse)
library(magrittr)
library(paletteer)
library(plotly)
library(vegan)
library(ggridges)

#####################
#########DATA########
#####################

#Import covariance matrix from pcangsd

cov <- read.table("data/ALL_CHR_singletons.cov")

#How many PC axis to retain?

pca <- eigen(cov)

pca.mat<-as.matrix(pca$vectors %*% (diag(pca$values))^0.5)

#add column names
colnames(pca.mat) <- paste0("PC", 1:ncol(pca.mat))

#add rownames
id <- read.table("data/id_part3.tsv")
rownames(pca.mat)<-id$V1

#Make dataset to plot
pca_coar <- data.frame(pca.mat[,1:20])

#sample id column
pca_coar$id <- id$V1

#pop column
metadata <- read.csv("data/data_with_gillrakers.csv",header=T)

pca_coar %<>% left_join(metadata, by = c("id" = "FISHES_id"))

pca_coar %<>% mutate(pop = if_else(is.na(pop),"RUP",pop))

pca_coar$pop <-factor(pca_coar$pop, levels = c("RUP","NOT_fall","NOT_summer"))

pca_coar

#add length info for other samples
metadata_other <- read.csv("data/gtscore_analysis/data_2020_2021.csv")


pca_coar %<>% left_join(metadata_other %>% dplyr::select(FISHES_ID,length), by = c("id" = "FISHES_ID") )

gwas_pheno <- pca_coar %>% mutate(length_all = if_else(is.na(long_fourche),length, long_fourche)) %>% dplyr::select(id,length_all) %>% distinct()

gwas_pheno_ord <- id %>% left_join(gwas_pheno, by = c("V1" = "id"))

colnames(gwas_pheno_ord) <- c("FID","length")

gwas_pheno_ord$IID <- gwas_pheno_ord$FID

write.table(gwas_pheno_ord[,c(1,3,2)], "data/gwas_length/pheno.sample", quote = F, row.names = F)

