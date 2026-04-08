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

win20kb <- fread("../12.2_outlier_detection/not_sum_lag.window_20kb.tsv")

nsnp_win20kb <- win20kb %>% group_by(chrom,window) %>% count() %>%  ungroup()

window2remove <- nsnp_win20kb %>% filter(n <= 10) %>% select(chrom,window)

win20kb_filt <- anti_join(win20kb,window2remove,by = c("chrom","window"))

thresh <- quantile(win20kb_filt$fst, probs = 0.99, na.rm = T)

win20kb_filt %<>% mutate(status = if_else(fst >= thresh, "outlier","background"))

win20kb_outlier_per_win <- win20kb_filt %>% group_by(chrom, window) %>% count(status)

win20kb_outlier_per_win %<>% 
  filter(!is.na(status)) %>% 
  pivot_wider(names_from = status, values_from = n) %>% 
  replace_na(list(outlier = 0)) %>% 
  mutate(nsnp = background + outlier)

# Calculate the overall proportion of outliers
p <- sum(win20kb_outlier_per_win$outlier) / sum(win20kb_outlier_per_win$nsnp)

# Generate binomial distribution of outliers
binom_dist_outliers <- sapply(1:max(win20kb_outlier_per_win$nsnp), function(j) {
  qbinom(0.9999, j, p)
})

binom_dist <- data.frame(outlier_null_dist = binom_dist_outliers, nsnp = 1:max(win20kb_outlier_per_win$nsnp))

# Identify top-candidate genes
top_candidate <- win20kb_outlier_per_win %>%
  left_join(binom_dist, by = "nsnp") %>%
  dplyr::mutate(top_candidate = if_else(outlier <= outlier_null_dist, "non_top_candidate", "top_candidate"))

top_candidate %>% ungroup() %>% filter(top_candidate == "top_candidate") %>% count()

#4038 windows enriched in outliers

wind20kb_top_candidate <- win20kb_filt %>% left_join(top_candidate %>% select(-background,-outlier,-nsnp, -outlier_null_dist),
                                                     by = c("chrom","window"))

fwrite(wind20kb_top_candidate,"../12.2_outlier_detection/notsum_lag_wind20kb_top_candidate.tsv",quote = F, sep = "\t")

top_windows_plot <- 
  wind20kb_top_candidate %>% 
  filter(top_candidate == "top_candidate") %>% 
  select(-pos,-fst,-window,-status) %>% 
  group_by(chrom, window_start, window_end) %>% 
  slice(1)


top_windows_plot %>% group_by(chrom) %>% count() %>% View()
