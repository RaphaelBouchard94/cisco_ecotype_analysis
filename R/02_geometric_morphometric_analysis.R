rm(list = ls())

#################
#### Library ####
#################

library(tidyverse)
library(magrittr)
library(geomorph)
library(abind)
library(Morpho)
library(MASS)
library(ggpubr)
library(plyr)
library(cowplot)
library(ggExtra)
library(png)
library(patchwork)
library(car)
library(multcomp)
library(RRPP)
library(caret)
library(klaR)
library(performance)
library(factoextra)
library(mclust)

##########################################
#MANOVA

final_data_gpa_projected <- readRDS("tmp/final_data_gpa_projected.Rds")

data <- read.csv("data/data_with_gillrakers.csv", header = T,sep = ",")

data %<>% mutate(river = if_else(str_detect(sample_id,'NOT'),"NOT","RUP"))

gdf <- geomorph.data.frame(shape = final_data_gpa_projected$coords,
                           pop = data$river,
                           sex = data$sex,
                           eco = data$true_eco,
                           cs = final_data_gpa_projected$Csize)

##############################################
###### PCA with corrected dataset ##########
##############################################


pca_gm <- gm.prcomp(final_data_gpa_projected$coords)

summary(pca_gm)

#Check proportion of variance explained by each PC
screeplot(pca_gm, type = "lines")



#Data exploration of PC axis

pca_ggplot <- data.frame(ecotype = data$true_eco,
                         pop = data$pop,
                         sex = data$sex,
                         eco = data$true_eco,
                         cs = data$long_fourche,
                         PC1 = pca_gm$x[,"Comp1"],
                         PC2 = pca_gm$x[,"Comp2"],
                         PC3 = pca_gm$x[,"Comp3"],
                         PC4 = pca_gm$x[,"Comp4"])


geo_pc12 <- pca_ggplot %>% filter(!sex == "R") %>% ggplot(aes(x = PC1 , y = PC2 , color = eco, shape = sex))+
  geom_point(size = 2)+
  xlab("PC1 (19.89%)")+
  ylab("PC2 (16.01%)")+
  stat_ellipse(aes(linetype = sex))+
  scale_color_manual(name = "Ecotype", labels = c("summer_run","fall-run"), values = c("#E69F00","skyblue1"))+
  scale_fill_manual(name = "Ecotype", labels = c("summer_run","fall-run"), values = c("#E69F00","skyblue1"))+
  scale_shape_manual(name = "Sex", labels = c("Female","Male"), values = c(1,16))+
  scale_linetype_manual(name = "Sex", labels = c("Female","Male"), values = c("solid","dotted"))+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "none")

geo_pc12

ggsave("figures/02_geomorph_pc12.png", unit = "cm", width = 10, height = 10, dpi = 600)

pca_ggplot %>% filter(!sex == "R") %>% ggplot(aes(x = PC1 , y = PC2 , color = eco, shape = sex))+
  geom_point(size = 2)+
  xlab("PC1 (19.89%)")+
  ylab("PC2 (16.01%)")+
  stat_ellipse(aes(linetype = sex))+
  scale_color_manual(name = "Ecotype", labels = c("summer_run","fall-run"), values = c("#E69F00","skyblue1"))+
  scale_fill_manual(name = "Ecotype", labels = c("summer_run","fall-run"), values = c("#E69F00","skyblue1"))+
  scale_shape_manual(name = "Sex", labels = c("Female","Male"), values = c(1,16))+
  scale_linetype_manual(name = "Sex", labels = c("Female","Male"), values = c("solid","dotted"))+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA))




geo_pc34 <- pca_ggplot %>% filter(!sex == "R") %>% ggplot(aes(x = PC3 , y = PC4 , color = eco, shape = sex))+
  geom_point(size = 2)+
  xlab("PC3 (12.16%)")+
  ylab("PC4 (7.37%)")+
  stat_ellipse(aes(linetype = sex))+
  scale_color_manual(name = "Ecotype", labels = c("Kâcikâsikumekw","Nûtamesânîw"), values = c("#E69F00","skyblue1"))+
  scale_fill_manual(name = "Ecotype", labels = c("Kâcikâsikumekw","Nûtamesânîw"), values = c("#E69F00","skyblue1"))+
  scale_shape_manual(name = "Sex", labels = c("Female","Male"), values = c(1,16))+
  scale_linetype_manual(name = "Sex", labels = c("Female","Male"), values = c("solid","dotted"))+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "none")

geo_pc34

ggsave("figures/02_geomorph_pc34.png", unit = "cm", width = 10, height = 10, dpi = 600)


geo_pc12 + geo_pc34

ggsave("figures/02_geomorph_pc1234.png", unit = "cm", width = 20, height = 10, dpi = 600)


#############################

# Describe visually morphological variation associated with each PC 

data_gpa_proj <- readRDS("tmp/data_gpa_proj.Rds")

PCA_proj <- readRDS("tmp/PCA_proj.Rds")

ref <- mshape(data_gpa_proj)

x<-c(2,1,4,5,6,7,8,9,10,11,12,13,14,15,16)
y<-c(1,4,5,6,7,8,9,10,11,12,13,14,15,16,2)

links <- data.frame(x,y)
links <- as.matrix(links)


##PC1

png("figures/02_PC1_morph_change.png",
    width = 3000, height = 3000,  
    res = 300,                     
    bg = "transparent")

plotRefToTarget(PCA_proj$shapes$shapes.comp1$min, 
                PCA_proj$shapes$shapes.comp1$max,
                method="points",
                mag=2,
                links = links,
                gridPars=gridPar(pt.bg = "grey40", 
                                 tar.pt.bg = "black",
                                 pt.size = 3,link.lwd = 3,
                                 link.col = "grey40",
                                 tar.link.col = "black", tar.link.lwd = 3,
                                 tar.pt.size = 3))

dev.off()

##PC2

png("figures/02_PC2_morph_change.png",
    width = 3000, height = 3000,  
    res = 300,                     
    bg = "transparent")

plotRefToTarget(PCA_proj$shapes$shapes.comp2$min, 
                PCA_proj$shapes$shapes.comp2$max, 
                method="points",mag=2, links = links,
                gridPars=gridPar(pt.bg = "grey40", 
                                 tar.pt.bg = "black",
                                 pt.size = 3,link.lwd = 3,
                                 link.col = "grey40",
                                 tar.link.col = "black", tar.link.lwd = 3,
                                 tar.pt.size = 3))

dev.off()

##PC3

png("figures/02_PC3_morph_change.png",
    width = 3000, height = 3000,  
    res = 300,                     
    bg = "transparent")

plotRefToTarget(PCA_proj$shapes$shapes.comp3$min,
                PCA_proj$shapes$shapes.comp3$max,
                method="points",mag=2, links = links,
                gridPars=gridPar(pt.bg = "grey40", 
                                 tar.pt.bg = "black",
                                 pt.size = 3,link.lwd = 3,
                                 link.col = "grey40",
                                 tar.link.col = "black", tar.link.lwd = 3,
                                 tar.pt.size = 3))

dev.off()

# #PC4

png("figures/02_PC4_morph_change.png",
    width = 3000, height = 3000,  
    res = 300,                     
    bg = "transparent")

plotRefToTarget(PCA_proj$shapes$shapes.comp4$min,
                PCA_proj$shapes$shapes.comp4$max,
                method="points",mag=2, links = links,
                gridPars=gridPar(pt.bg = "grey40", 
                                 tar.pt.bg = "black",
                                 pt.size = 3,link.lwd = 3,
                                 link.col = "grey40",
                                 tar.link.col = "black", tar.link.lwd = 3,
                                 tar.pt.size = 3))

dev.off()
##############################
#Cross-validated LDA on PCA scores


pcs <- pca_gm$x[,1:4]

lda_cv <- lda(
  x = pcs,
  grouping = gdf$eco,
  CV = TRUE
)

table(True = gdf$eco, Predicted = lda_cv$class)
mean(lda_cv$class == gdf$eco)

conf <- table(True = gdf$eco, Predicted = lda_cv$class)

sens <- diag(prop.table(conf, 1))

balanced_acc <- mean(sens)

