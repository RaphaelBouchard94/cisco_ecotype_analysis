

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

win05kb <- fread("data/genome_scan/notsummer_rup_5kb.outlier.fst")
win10kb <- fread("data/genome_scan/notsummer_rup_10kb.outlier.fst")
win15kb <- fread("data/genome_scan/notsummer_rup_15kb.outlier.fst")
win20kb <- fread("data/genome_scan/notsummer_rup_20kb.outlier.fst")
win25kb <- fread("data/genome_scan/notsummer_rup_25kb.outlier.fst")


nsnp_win05kb <- win05kb %>% group_by(chrom,window) %>% count() %>% ungroup()

less10snp_05kb <- nsnp_win05kb %>% filter(n < 10) %>% count()

less10snp_05kb/nrow(nsnp_win05kb)

#32% window with less than 10 SNP

hist(nsnp_win10kb$n,breaks = 100)
abline(v = 10)


#10kb windows

nsnp_win10kb <- win10kb %>% group_by(chrom,window) %>% count() %>% ungroup()

less10snp_10kb <- nsnp_win10kb %>% filter(n < 10) %>% count()

less10snp_10kb/nrow(nsnp_win10kb)

#12% windows with less than 10 SNP

hist(nsnp_win20kb$n,breaks = 1000)
abline(v = 10)

nsnp_win15kb <- win15kb %>% group_by(chrom,window) %>% count() %>% ungroup()

less10snp_15kb <- nsnp_win15kb %>% filter(n < 10) %>% count()

less10snp_15kb/nrow(nsnp_win15kb)

#7% windows with less than 10 SNPs

hist(nsnp_win10kb$n,breaks = 100)
abline(v = 10, col = "red")

nsnp_win20kb <- win20kb %>% group_by(chrom,window) %>% count() %>% ungroup()

less10snp_20kb <- nsnp_win20kb %>% filter(n < 10) %>% count()

less10snp_20kb/nrow(nsnp_win20kb)

#6% windows with less than 10 SNPs

nsnp_win25kb <- win25kb %>% group_by(chrom,window) %>% count() %>% ungroup()

less10snp_25kb <- nsnp_win25kb %>% filter(n < 10) %>% count()

less10snp_25kb/nrow(nsnp_win25kb)

#5% windows with less than 10 SNPs


#Let's keep 20kb windows.



