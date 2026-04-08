setwd("~/2024-11-01_cisco_james_bay/01_scripts/")

#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(magrittr)
library(data.table)
library(patchwork)

# Function to subset SNPs within outlier windows
subset_snps <- function(snps, windows) {
  # Initialize an empty data frame for results
  snps_within_windows <- data.frame()
  
  # Loop over each outlier window
  for (i in 1:nrow(windows)) {
    window_chr <- windows$chrom[i]
    window_start <- windows$window_start[i]
    window_end <- windows$window_end[i]
    
    # Find SNPs that match the chromosome and fall within the window
    subset <- snps %>%
      filter(chrom == window_chr & pos >= window_start & pos <= window_end)
    
    # Append to the results
    snps_within_windows <- rbind(snps_within_windows, subset)
  }
  
  return(snps_within_windows)
}

#------------#
#--- DATA ---#
#------------#

outlier_windows <- fread("../12.2_outlier_detection/common_outlier_window_merged.tsv")

outlier_notsum_notfall <- fread("../12.2_outlier_detection/notsum_fall_wind20kb_top_candidate.tsv")
outlier_notsum_lag <- fread("../12.2_outlier_detection/notsum_lag_wind20kb_top_candidate.tsv")

#Subset outlier positions within outlier overlapping windows

outlier_notsum_notfall %<>% filter(status == "outlier", top_candidate == "top_candidate") %>% select("chrom","pos")
outlier_notsum_lag %<>% filter(status == "outlier", top_candidate == "top_candidate") %>% select("chrom","pos")

all_outliers_snp <- rbind(outlier_notsum_notfall,outlier_notsum_lag) %>% unique()

filtered_snps <- subset_snps(all_outliers_snp, outlier_windows)

maf_notsum <- fread("../11_saf_maf_by_pop/NOT_summer.mafs") 
maf_notfall <- fread("../11_saf_maf_by_pop/NOT_fall.mafs")
maf_lag <- fread("../11_saf_maf_by_pop/LAG.mafs")

maf_notsum %<>% select(chromo, position, knownEM) %>% rename(NOTsummer = knownEM)
maf_notfall %<>% select(chromo, position, knownEM) %>% rename(NOTfall = knownEM)
maf_lag %<>% select(chromo, position, knownEM) %>% rename(LAG = knownEM)

all_freq <- left_join(left_join(maf_notsum, maf_notfall, by = c("chromo","position")), maf_lag, by = c("chromo","position"))

all_freq_filt <- left_join(filtered_snps, all_freq, by = c("chrom" = "chromo","pos" = "position"))

all_freq_filt %<>% pivot_longer(NOTsummer:LAG, names_to = "pop") %>% rename(all_freq = value)

#Plot

axis_man <- all_freq_filt %>% 
  group_by(chrom) %>%
  summarize(center = mean(pos))

ggplot(all_freq_filt, aes(x = pos , y = all_freq, color = as_factor(pop))) +
  geom_point(alpha = 0.2, size=0.5) +
  facet_wrap(~chrom, nrow = 8, ncol = 5, scales = "free_x") +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.05,1.05))+
  labs(
    x = "Position",
    y = "Allele frequencies") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank())

#CM078251


ggplot(all_freq_filt %>% filter(chrom == "CM078251.1", pos > 34000000), aes(x = pos , y = all_freq, color = as_factor(pop))) +
  geom_point(alpha = 0.2, size=0.5) +
  facet_wrap(~chrom, nrow = 8, ncol = 5, scales = "free_x") +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.05,1.05))+
  labs(
    x = "Position",
    y = "Allele frequencies") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank())


