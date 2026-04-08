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

win10kb <- fread("../12.2_outlier_detection/not_fall_sum.window_10kb.tsv")
win20kb <- fread("../12.2_outlier_detection/not_fall_sum.window_20kb.tsv")
win50kb <- fread("../12.2_outlier_detection/not_fall_sum.window_50kb.tsv")

nsnp_win10kb <- win10kb %>% group_by(chrom,window) %>% count() %>% ungroup()

nsnp_win10kb %>% filter(n < 10) %>% count()

16230/217942

hist(nsnp_win10kb$n,breaks = 100)
abline(v = 10)

nsnp_win20kb <- win20kb %>% group_by(chrom,window) %>% count() %>% ungroup()

nsnp_win20kb %>% filter(n < 10) %>% count()

4915/112205

hist(nsnp_win20kb$n,breaks = 1000)
abline(v = 10)

nsnp_win50kb <- win50kb %>% group_by(chrom,window) %>% count() %>% ungroup()

nsnp_win50kb %>% filter(n < 10) %>% count()

1472/46059

hist(nsnp_win50kb$n,breaks = 100)
abline(v = 10)

#Windows of 20kb seems to be the most legit with less than 5% of windows
#with less than 10 SNP








