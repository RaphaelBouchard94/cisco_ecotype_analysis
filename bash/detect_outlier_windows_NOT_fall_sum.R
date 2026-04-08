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

win20kb <- fread("../12.2_outlier_detection/not_fall_sum.window_20kb.tsv")

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

top_candidate %>% ungroup() %>% filter(top_candidate == "top_candidate") %>% View()

#4023 windows enriched in outliers

wind20kb_top_candidate <- win20kb_filt %>% left_join(top_candidate %>% select(-background,-outlier,-nsnp, -outlier_null_dist),
                                      by = c("chrom","window"))

fwrite(wind20kb_top_candidate,"../12.2_outlier_detection/notsum_fall_wind20kb_top_candidate.tsv",quote = F, sep = "\t")


top_windows_plot <- 
  wind20kb_top_candidate %>% 
  filter(top_candidate == "top_candidate") %>% 
  select(-pos,-fst,-window,-status) %>% 
  group_by(chrom, window_window_start, window_window_end) %>% 
  slice(1)


top_windows_plot %>% group_by(chrom) %>% count()

merge_windows <- function(df, max_gap) {
  df <- df %>%
    arrange(chrom, window_start) # Ensure data is sorted
  
  # Initialize output
  merged <- data.frame()
  
  # Variables to keep track of current region
  current_chrom <- df$chrom[1]
  current_window_start <- df$window_start[1]
  current_window_end <- df$window_end[1]
  current_status <- df$top_candidate[1]
  
  for (i in 2:nrow(df)) {
    if (df$chrom[i] == current_chrom && df$window_start[i] - current_window_end <= max_gap) {
      # Extwindow_end the current region
      current_window_end <- max(current_window_end, df$window_end[i])
    } else {
      # Save the previous region and window_start a new one
      merged <- rbind(merged, data.frame(chrom = current_chrom, window_start = current_window_start, window_end = current_window_end, status = current_status))
      current_chrom <- df$chrom[i]
      current_window_start <- df$window_start[i]
      current_window_end <- df$window_end[i]
      current_status <- df$top_candidate[i]
    }
  }
  
  # Save the last region
  merged <- rbind(merged, data.frame(chrom = current_chrom, window_start = current_window_start, window_end = current_window_end, status = current_status))
  
  return(merged)
}

# Apply the function
merged_data <- merge_windows(top_windows_plot, max_gap = 50000)

merged_data %>% group_by(chrom) %>% count() %>% View()


#Plot manhattan with outlier regions

NOT_fall_sum <- fread("../12_fst/NOT_fall_NOT_summer.slidingwindow")[,-1]

colnames(NOT_fall_sum) <- c("chrom","pos","nsites","fst")

axis_man <- NOT_fall_sum %>% 
  group_by(chrom) %>%
  summarize(center = mean(pos))

ggplot(NOT_fall_sum, aes(x = pos , y = fst, color = as_factor(chrom))) +
  geom_point(alpha = 0.2, size=0.5) +
  geom_smooth(color = "black",linewidth = .5)+
  geom_rect(data=merged_data, 
            inherit.aes=FALSE, 
            aes(xmin=window_start, xmax=window_end,
                ymin=0,ymax=1, 
                group=chrom), 
            color="transparent", 
            fill="darkred", 
            alpha=0.3)+
  scale_x_continuous(
    label = axis_man$chrom,
    breaks = axis_man$center,
    expand = c(0.05,0)) +
  scale_color_manual(values = rep(
    c("grey70", "grey90"),
    unique(length(axis_man$chrom)))) +
  facet_grid(~chrom, scales="free_x") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  labs(
    x = "NOT_summer - NOT_fall",
    y = "Fst") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0,'lines'),
        axis.title.y =  element_text(angle = 90, size = 12, vjust = 0.5),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90))

ggsave("../12.1_fst_manhattan/NOT_sum_fall_outlier.jpeg")


ggplot(wind20kb_top_candidate %>% filter(chrom == "CM078251.1"), 
       aes(x = pos , y = fst, color = as_factor(chrom))) +
  geom_point(alpha = 0.2, size=0.5) +
  geom_smooth(color = "black",linewidth = .5)+
  geom_rect(data=merged_data %>% filter(chrom == "CM078251.1"), 
            inherit.aes=FALSE, 
            aes(xmin=window_start, xmax=window_end,
                ymin=0,ymax=1, 
                group=chrom), 
            color="transparent", 
            fill="indianred", 
            alpha=0.3)+
  scale_x_continuous(
    label = axis_man$chrom,
    breaks = axis_man$center,
    expand = c(0.05,0)) +
  scale_color_manual(values = rep(
    c("grey70", "grey90"),
    unique(length(axis_man$chrom)))) +
  facet_grid(~chrom, scales="free_x") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  labs(
    x = "NOT_summer - NOT_fall",
    y = "Fst") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0,'lines'),
        axis.title.y =  element_text(angle = 90, size = 12, vjust = 0.5),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90))









