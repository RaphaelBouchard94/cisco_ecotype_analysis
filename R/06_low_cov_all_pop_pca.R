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

cov <- read.table("data/pca/all_maf0.01pctind0.8_maxdepth10_ALL_CHR_singletons.canonical.pruned.cov")

#How many PC axis to retain?

pca <- rda(cov)
screeplot(pca, bstick = TRUE, type = "lines")

pca <- eigen(cov)

pca.mat<-as.matrix(pca$vectors %*% (diag(pca$values))^0.5)

#add column names
colnames(pca.mat) <- paste0("PC", 1:ncol(pca.mat))

var_explained <- round(pca$values * 100 / sum(pca$values[pca$values >= 0]), 2)

# Then access as needed:
var1 <- var_explained[1]
var2 <- var_explained[2]
var3 <- var_explained[3]
var4 <- var_explained[4]

#Make dataset to plot
pca_coar <- data.frame(pca.mat[,1:20])

#add id column
metadata <- read.table("data/id_lowcov_fulldata.tsv")

pca_coar$fishes_id <- metadata$V1

#add pop columns

metadata2 <- read.table("data/ngsadmix/admixture_assign_k4_allpop.txt",header=T)

pca_coar %<>% left_join(metadata2, by = c("fishes_id" = "individual_id"))

#sanity check

pca_coar %>% group_by(pop) %>% count()

#remove na

pca_coar %<>% filter(!is.na(pop))
##add info on pop
# pop_file <- read.table("data/pop_coordinates.txt",header=T)
# 
# pop_file %>% arrange(desc(lat)) 
# 
# pca_coar %<>% mutate(pop_tmp = if_else(pop %in% c("NOT_fall","NOT_summer"),"NOT",pop))
# 
# pca_coar %<>% left_join(pop_file,by = c("pop_tmp" = "pop_code"))
# 
# pca_coar %>% group_by(pop_tmp) %>% count()
# 
# pca_coar$pop.x <-factor(pca_coar$pop.x, 
#                         levels = c("LAG","SAB","EAS","RUP","NOT_fall","NOT_summer"))


ggplot(pca_coar, aes(x = PC1,y = PC2,color = genetic_cluster, alpha = assignment_probability))+
  geom_point(size = 4)+
  xlab("PC1 (4.44%)")+
  ylab("PC2 (2.58%)")+
  scale_color_manual(values = c("#8B8989","#E69F00","#0072B2","skyblue1"))+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = c(0.8,0.25),
        legend.title = element_blank())+
  guides(alpha = "none")


ggsave("figures/06_allpop_pc1_pc2_singletons_ld_pruned.jpeg",height=10,width =10, unit = "cm", dpi = 600)

ggplot(pca_coar, aes(x = PC3,y = PC4,color = genetic_cluster, alpha = assignment_probability))+
  geom_point(size = 4)+
  xlab("PC3 (2.16%)")+
  ylab("PC4 (1.08%)")+
  scale_color_manual(values = c("#8B8989","#E69F00","#0072B2","skyblue1"))+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "none")

ggsave("figures/06_allpop_pc3_pc4_singletons_ld_pruned.jpeg",height=10,width =10, unit = "cm", dpi = 600)

#######################################
### Check if pattern of missingness ###
#######################################

# missing <- read.table("00_data/pca/missing_per_ind.txt")
# 
# pca_coar$missing <- missing$V2
# 
# ggplot(pca_coar, aes(x = -PC1,y = PC2,color = missing))+
#   geom_point(alpha=0.4,size = 3)+
#   xlab("PC1 (2.96%)")+
#   ylab("PC2 (0.94%)")+
#   scale_color_viridis_c()+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.title.x = element_text(size = 18),
#         axis.title.y = element_text(size = 18),
#         axis.text.x = element_text(size = 18),
#         axis.text.y = element_text(size = 18),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 18))
# 
# ggplot(pca_coar, aes(x = PC3,y = PC4,color = missing))+
#   geom_point(alpha=0.8,size = 3)+
#   xlab("PC1 (2.96%)")+
#   ylab("PC2 (0.94%)")+
#   scale_color_viridis_c()+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.title.x = element_text(size = 18),
#         axis.title.y = element_text(size = 18),
#         axis.text.x = element_text(size = 18),
#         axis.text.y = element_text(size = 18),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 18))
# 
# 
# ggplot(pca_cocl, aes(x = fct_reorder(pop,desc(lat)), y = missing))+
#   geom_boxplot()+
#   xlab("Population (ordered by latitude)")+
#   ylab("Proportion of missing data")+
#   theme_bw()
# 
# ggsave("figure/missing_data_per_pop.jpeg")
# 
# mean(pca_coar$missing)
# mean(pca_coar$missing) + 2*sd(pca_coar$missing)
# 
# pca_coar[which(pca_coar$missing > mean(pca_coar$missing) + 2*sd(pca_coar$missing)),]
# 
# 
