rm(list = ls())

########################
# 1 - Library  #
########################

library(tidyverse)
library(magrittr)
library(data.table)
library(ggdist)
library(windowscanr)

#######################
######## Data #########
#######################
hobs_aa <- fread("data/chr34/aa_geno_CM078251.1.hwe")
hobs_ab <- fread("data/chr34/ab_geno_CM078251.1.hwe")
hobs_bb <- fread("data/chr34/bb_geno_CM078251.1.hwe")

hobs <- rbind(hobs_aa, hobs_ab, hobs_bb)

hobs$karyo <- c(rep("aa",nrow(hobs_aa)),rep("ab",nrow(hobs_ab)),rep("bb",nrow(hobs_bb)))

hobs$Hexp<-2*(hobs$hweFreq)*(1-hobs$hweFreq) #Hexp = 2*(hweFreq)(1-hweFreq)
hobs$Hobs<-hobs$Hexp-(hobs$F*hobs$Hexp) #Hobs= Hexp - F* Hexp


#Regarder distribution des Hobs pour chaque dataset

targ_win <- c(34800000,35619999)
targ_win_plus <- c(34800000-100000,35619999+100000)

hobs %>% filter(karyo == "aa")

hobsaa <- winScan(x = hobs %>% filter(karyo == "aa"), groups = "Chromo", position = "Position",values = "Hobs",win_size = 20000,win_step = 1000,funs = c("mean"))

hobsab <- winScan(x = hobs %>% filter(karyo == "ab"), groups = "Chromo", position = "Position",values = "Hobs",win_size = 20000,win_step = 1000,funs = c("mean"))

hobsbb <- winScan(x = hobs %>% filter(karyo == "bb"), groups = "Chromo", position = "Position",values = "Hobs",win_size = 20000,win_step = 1000,funs = c("mean"))


winhobs <- rbind(hobsaa,hobsab, hobsbb)

winhobs$karyo <- c(rep("aa",nrow(hobsaa)),rep("ab",nrow(hobsab)),rep("bb",nrow(hobsbb)))


ggplot(winhobs %>% filter(win_mid > targ_win_plus[1] & win_mid < targ_win_plus[2]),
       aes(x = fct_reorder(karyo,desc(karyo)), y = Hobs_mean, fill = karyo)) +
  geom_boxplot(width = 0.5, alpha = 0.8)+
  scale_fill_manual(values = c("grey90", "grey60", "grey10"))+
  ylab("Hobs") +
  xlab("Putative karyotype") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = "none")

ggsave("figures/05_het_chr36_outlier.png",
       width = 10, height = 10, unit = "cm", dpi = 900)


#######################################











