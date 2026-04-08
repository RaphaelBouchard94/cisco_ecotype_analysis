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
var3 <- var_explained[3]
var4 <- var_explained[4]

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

##add info on pop

ggplot(pca_coar, aes(x = PC1,y = PC2,color = pop ))+
  geom_point(alpha=0.8,size = 4)+
  xlab("PC1 (8.39%)")+
  ylab("PC2 (1.15%)")+
  scale_color_manual(values = c("skyblue3", "skyblue1","#E69F00"))+
  theme_bw(base_size = 14)+
  theme(panel.border = element_rect(colour = "black",fill=NA),
        panel.grid = element_blank(),
        legend.position = c(0.7,0.25))


ggsave("figures/04_pca.png",
       width = 10, height = 10, unit = "cm", dpi = 900)



pca_coar %>% filter(PC1 > 0 & pop == "NOT_fall")



