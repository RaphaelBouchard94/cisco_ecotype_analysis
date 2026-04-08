setwd("~/2024-11-01_cisco_james_bay/01_scripts/")

#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(magrittr)
library(data.table)
library(patchwork)

#------------#
#--- DATA ---#
#------------#

#By pos

NOT_fall_sum_pos <- fread("../12_fst/NOT_fall_NOT_summer.bypos.sfs")

colnames(NOT_fall_sum_pos) <- c("chrom","pos","V1","V2")

NOT_fall_sum_pos %<>% mutate(fst = round(V1/V2,digits = 5)) %>% select(-V1,-V2)

fwrite(NOT_fall_sum_pos, 
       file = "../12.2_outlier_detection/not_fall_sum_fst.tsv",
       quote = F, col.names = T, row.names = F,sep = "\t")

NOT_sum_LAG_pos <- fread("../12_fst/LAG_NOT_summer.bypos.sfs")

colnames(NOT_sum_LAG_pos) <- c("chrom","pos","V1","V2")

NOT_sum_LAG_pos %<>% mutate(fst = round(V1/V2,digits = 5)) %>% select(-V1,-V2)

fwrite(NOT_sum_LAG_pos, 
       file = "../12.2_outlier_detection/not_sum_lag_fst.tsv",
       quote = F, col.names = T, row.names = F,sep = "\t")




