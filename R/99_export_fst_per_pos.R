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

NOT_fall_sum_pos <- fread("data/genome_scan/notsummer_rup.bypos.sfs")

colnames(NOT_fall_sum_pos) <- c("chrom","pos","V1","V2")

NOT_fall_sum_pos %<>% mutate(fst = round(V1/V2,digits = 5)) %>% select(-V1,-V2)

fwrite(NOT_fall_sum_pos, 
       file = "data/genome_scan/notsummer_rup_fst.tsv",
       quote = F, col.names = T, row.names = F,sep = "\t")


