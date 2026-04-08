rm(list = ls())

########################
# 1 - Library  #
########################

library(tidyverse)
library(magrittr)
library(data.table)
library(readxl)

#######################
######## Data #########
#######################

#######################
#Metadata

data_20_21 <- read_csv("data/gtscore_analysis/data_2020_2021.csv")
data_22 <- read_csv("data/gtscore_analysis/data_2022.csv")

data_20_21 %<>% dplyr::select(FISHES_ID,length,stock)
data_22 %<>% dplyr::select(FISHES_ID,length,stock)

metadata <- rbind(data_20_21,data_22)

metadata_filt <- metadata %>% filter(length_corr > 100, stock == "Source")

ggplot(metadata_filt , aes(x = length_corr))+
  geom_histogram()

##################################################################
#Cluster info

source <- read.table("data/ngsadmix/admixture_assign_k3_allpop.txt",header=T)

source %<>% dplyr::select(individual_id,genetic_cluster)

colnames(source) <- c("indiv","pop_clus")

source %<>% mutate(indiv = str_replace(indiv,"-","_"))


#######################


tmp <-  read.table("data/gtscore_analysis/polyGenResults_singleSNP_rubias.txt",header=T, 
                   stringsAsFactors=FALSE, 
                   colClasses = c("character"))


tmp %>% filter(is.na(sample_type)) %>% select(indiv)

tmp$indiv <- str_replace(tmp$indiv, 
                         "([A-Za-z]*s)_(\\d+)_", 
                         function(x) {
                           prefix <- str_extract(x, "[A-Za-z]*s(?=_)")
                           num <- str_extract(x, "(?<=[A-Za-z]s_)\\d+")
                           paste0(prefix, "_", sprintf("%03d", as.numeric(num)), "_")
                         })

tmp %<>% mutate(sample_type = if_else(is.na(sample_type),"reference",sample_type))


tmp %<>% 
  left_join(source, by = c("indiv"))

#Retain only the information on causal SNP

geno_all <- tmp %>% dplyr::select(indiv,starts_with("Chr34_34060819"),pop_clus)

geno_all %<>% mutate(geno = paste(Chr34_34060819_1,Chr34_34060819_1.1, sep = ""))

geno_pheno <- geno_all %>% 
  left_join(metadata_filt, by = c("indiv" = "id_corr")) %>%
  filter(!is.na(length),geno != "NANA", !is.na(pop_clus))


ggplot(geno_pheno, aes(x = geno, y = length_corr, color = pop_clus))+
  geom_boxplot(outliers = F)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
              alpha = 0.5)+
  labs(x = "MC4R Genotype",
       y = "Fork Length (mm)")+
  scale_color_manual(values = c("#0072B2","skyblue1","#E69F00"))+
  theme_bw(base_size = 16)+
  theme(legend.position = "none",
        panel.grid = element_blank())


aov_mc4r <- lm(length_corr ~ geno + pop_clus, data = geno_pheno)

summary(aov_mc4r)





test <- tmp %>% dplyr::select(indiv,starts_with("Chr34_34060819"))












