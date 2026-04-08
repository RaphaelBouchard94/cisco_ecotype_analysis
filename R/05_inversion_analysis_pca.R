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
library(apcluster)
library(plyr)

#####################
#########DATA########
#####################

###########################################
#ALL WINDOW
#Import covariance matrix from pcangsd

cov <- read.table("data/chr34/all_window.cov")

#How many PC axis to retain?

pca <- rda(cov)
screeplot(pca, bstick = TRUE, type = "lines")

pca <- eigen(cov)

pca.mat<-as.matrix(pca$vectors %*% (diag(pca$values))^0.5)

#add column names
colnames(pca.mat) <- paste0("PC", 1:ncol(pca.mat))

#add rownames
id <- read.table("data/id_part3.tsv")
rownames(pca.mat)<-id$V1

var_explained <- round(pca$values * 100 / sum(pca$values[pca$values >= 0]), 2)

# Then access as needed:
var1 <- var_explained[1]
var2 <- var_explained[2]

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

pca_coar %<>% mutate(true_eco = if_else(is.na(true_eco),"fall_run",true_eco))

##add info on pop

ggplot(pca_coar , aes(x = PC1,y = PC2,color = true_eco))+
  geom_point(alpha=0.8,size = 4)+
  xlab("PC1 (69.96%)")+
  ylab("PC2 (7.06%)")+
  scale_color_manual(values = c("skyblue1", "#E69F00"))+
  theme_bw(base_size = 16)+
  theme(panel.border = element_rect(colour = "black",fill=NA),
        panel.grid = element_blank(),
        legend.position = "none")


ggsave("figures/05_pca_chr36_outlier_window.png",
       width = 10, height = 10, unit = "cm", dpi = 900)

###############################
#APcluster to define genotypes
apclus_df <- pca_coar[,1]

s1 <- negDistMat(apclus_df , r=2)

apres1b <- apcluster(s1)

plot(apres1b, pca_coar[,1:2])

###############################

tidy_apclust_res <- reshape2::melt(apres1b@clusters)

tidy_apclust_res$value <- as.factor(tidy_apclust_res$value)

pca_coar$rowid <- as.factor(row.names(pca_coar))

plot_apclus <- left_join(pca_coar,tidy_apclust_res,by=c("rowid"="value"))

plot_apclus$L1 <- as.factor(plot_apclus$L1)

plot_apclus %<>% mutate(geno_clus = case_when(L1 == "1" ~ "ab",
                                             L1 == "2" ~ "aa",
                                             L1 == "3" ~ "bb"))

find_hull <- function(df)df[chull(df$PC1, df$PC2), ]
hulls <- ddply(plot_apclus, "geno_clus", find_hull)

ggplot()+
  geom_polygon(data = hulls,aes(x=PC1,y=PC2, fill = geno_clus),alpha = 0.6)+
  geom_point(data = plot_apclus, aes(x=PC1,y=PC2,colour = true_eco),alpha=0.8,size = 3)+
  scale_color_manual(values = c("skyblue1", "#E69F00"))+
  scale_fill_manual(values = c("grey90", "grey60", "grey10"))+
  xlab("PC1 (69.96%)")+
  ylab("PC2 (7.06%)")+
  labs(fill='Cluster',color = "Ecotype")+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank(),
        legend.position = "none")

ggsave("figures/05_local_pca_chr36_apclus.png",
       width = 10, height = 10, unit = "cm", dpi = 900)



write.table(aa_geno,"data/chr34/aa_geno_bam.id",quote = F, col.names = F, row.names = F)

ab_geno <- plot_apclus %>% filter(geno_clus == "ab") %>% pull(id)

write.table(ab_geno,"data/chr34/ab_geno_bam.id",quote = F, col.names = F, row.names = F)

bb_geno <- plot_apclus %>% filter(geno_clus == "bb") %>% pull(id)

write.table(bb_geno,"data/chr34/bb_geno_bam.id",quote = F, col.names = F, row.names = F)


karyo_file <- plot_apclus %>% dplyr::select(id,geno_clus)

write.table(karyo_file,"data/chr34/karyo_info.tsv",quote = F, row.names = F)

########################################
#Karyotype length/weight
data_20_21 <- read_csv("data/gtscore_analysis/data_2020_2021.csv")

data_20_21 %<>% dplyr::select(FISHES_ID,length,stock)

metadata <- data_20_21

apclus_meta <- plot_apclus %>% left_join(metadata, by =c("id" = "FISHES_ID"))

apclus_meta %<>% mutate(length_clean = if_else(is.na(long_fourche),length,long_fourche))

apclus_meta  %>% 
  ggplot(aes(x = fct_reorder(geno_clus,desc(geno_clus)),y = length_clean, fill = geno_clus))+
  geom_violin(width = 0.5,alpha = 0.6)+
  stat_summary(
    fun = mean,
    geom = "point",
    size = 1,
    color = "black"
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.1, 
    color = "black"
  ) +
  #geom_jitter(width = 0.1, alpha = 0.5)+
  scale_fill_manual(values = c("grey90", "grey60", "grey10"))+
  ylab("Fork length (mm)") +
  xlab("Putative karyotype") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = "none")


ggsave("figures/05_karyo_length_chr36_outlier.png",
       width = 10, height = 10, unit = "cm", dpi = 900)



###############################################################################
#First section

#Import covariance matrix from pcangsd

cov2 <- read.table("data/chr34/first_section_window.cov")

#How many PC axis to retain?

pca2 <- eigen(cov2)

pca.mat2<-as.matrix(pca2$vectors %*% (diag(pca2$values))^0.5)

#add column names
colnames(pca.mat2) <- paste0("PC", 1:ncol(pca.mat2))

#add rownames
id <- read.table("data/id_part3.tsv")
rownames(pca.mat2)<-id$V1

var_explained <- round(pca2$values * 100 / sum(pca2$values[pca2$values >= 0]), 2)

# Then access as needed:
var1.2 <- var_explained[1]
var2.2 <- var_explained[2]

#Make dataset to plot
pca_coar_first_section <- data.frame(pca.mat2[,1:20])

#sample id column
pca_coar_first_section$id <- id$V1

pca_coar_first_section %<>% left_join(metadata, by = c("id" = "FISHES_id"))

pca_coar_first_section %<>% mutate(pop = if_else(is.na(pop),"RUP",pop))

pca_coar_first_section$pop <-factor(pca_coar_first_section$pop, levels = c("RUP","NOT_fall","NOT_summer"))

pca_coar_first_section

##add info on pop

ggplot(pca_coar_first_section, aes(x = PC1,y = PC2,color = pop))+
  geom_point(alpha=0.8,size = 4)+
  xlab("PC1 (69.96%)")+
  ylab("PC2 (7.06%)")+
  scale_color_manual(values = c("#56B4E9", "#0072B2","#E69F00"))+
  theme_bw(base_size = 16)+
  theme(panel.border = element_rect(colour = "black",fill=NA),
        panel.grid = element_blank())


ggsave("figures/05_pca_chr34_first_section.png",
       width = 10, height = 10, unit = "cm", dpi = 900)


#################################################################################
#Last section

###############################################################################
#First section

#Import covariance matrix from pcangsd

cov3 <- read.table("data/chr34/outlier_last_section.cov")

#How many PC axis to retain?

pca3 <- eigen(cov3)

pca.mat3<-as.matrix(pca3$vectors %*% (diag(pca3$values))^0.5)

#add column names
colnames(pca.mat3) <- paste0("PC", 1:ncol(pca.mat3))

#add rownames
id <- read.table("data/id_part3.tsv")
rownames(pca.mat3)<-id$V1

var_explained3 <- round(pca3$values * 100 / sum(pca3$values[pca3$values >= 0]), 2)

# Then access as needed:
var1.3 <- var_explained3[1]
var2.3 <- var_explained3[2]

#Make dataset to plot
pca_coar_last_section <- data.frame(pca.mat3[,1:20])

#sample id column
pca_coar_last_section$id <- id$V1

pca_coar_last_section%<>% left_join(metadata, by = c("id" = "FISHES_id"))

pca_coar_last_section%<>% mutate(pop = if_else(is.na(pop),"RUP",pop))

pca_coar_last_section$pop <-factor(pca_coar_first_section$pop, levels = c("RUP","NOT_fall","NOT_summer"))


##add info on pop

ggplot(pca_coar_last_section, aes(x = PC1,y = PC2,color = pop))+
  geom_point(alpha=0.8,size = 4)+
  xlab("PC1 (69.96%)")+
  ylab("PC2 (7.06%)")+
  scale_color_manual(values = c("#56B4E9", "#0072B2","#E69F00"))+
  theme_bw(base_size = 16)+
  theme(panel.border = element_rect(colour = "black",fill=NA),
        panel.grid = element_blank(),
        legend.position = "none")


ggsave("figures/05_pca_chr34_second_section.png",
       width = 10, height = 10, unit = "cm", dpi = 900)

