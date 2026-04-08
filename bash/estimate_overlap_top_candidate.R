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

tc_notsum_fall <- fread("../12.2_outlier_detection/notsum_fall_wind20kb_top_candidate.tsv")
tc_notsum_lag <- fread("../12.2_outlier_detection/notsum_lag_wind20kb_top_candidate.tsv")


top_windows_notsum_fall <- 
  tc_notsum_fall %>% 
  filter(top_candidate == "top_candidate") %>% 
  select(-pos,-fst,-window,-status) %>% 
  group_by(chrom, window_start, window_end) %>% 
  slice(1)

top_windows_notsum_lag <- 
  tc_notsum_lag %>% 
  filter(top_candidate == "top_candidate") %>% 
  select(-pos,-fst,-window,-status) %>% 
  group_by(chrom, window_start, window_end) %>% 
  slice(1)


common_top_windows <- 
  left_join(top_windows_notsum_fall,top_windows_notsum_lag,by = c("chrom","window_start","window_end")) %>% 
  filter(!is.na(top_candidate.y))

common_top_windows %>% group_by(chrom) %>% count() %>% View()

merge_windows <- function(df, max_gap) {
  df <- df %>%
    arrange(chrom, window_start) # Ensure data is sorted
  
  # Initialize output
  merged <- data.frame()
  
  # Variables to keep track of current region
  current_chrom <- df$chrom[1]
  current_window_start <- df$window_start[1]
  current_window_end <- df$window_end[1]
  
  for (i in 2:nrow(df)) {
    if (df$chrom[i] == current_chrom && df$window_start[i] - current_window_end <= max_gap) {
      # Extwindow_end the current region
      current_window_end <- max(current_window_end, df$window_end[i])
    } else {
      # Save the previous region and window_start a new one
      merged <- rbind(merged, data.frame(chrom = current_chrom, window_start = current_window_start, window_end = current_window_end))
      current_chrom <- df$chrom[i]
      current_window_start <- df$window_start[i]
      current_window_end <- df$window_end[i]
    }
  }
  
  # Save the last region
  merged <- rbind(merged, data.frame(chrom = current_chrom, window_start = current_window_start, window_end = current_window_end))
  
  return(merged)
}

# Apply the function
merged_data <- merge_windows(common_top_windows, max_gap = 1000000)

fwrite(merged_data,"../12.2_outlier_detection/common_outlier_window_merged.tsv",quote = F, sep = "\t")

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
            alpha=0.5)+
  scale_x_continuous(
    label = axis_man$chrom,
    breaks = axis_man$center,
    expand = c(0.05,0)) +
  scale_color_manual(values = rep(
    c("grey70", "grey85"),
    unique(length(axis_man$chrom)))) +
  facet_wrap(~chrom, nrow = 8, ncol = 5,scales="free_x") +
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

ggsave("../12.1_fst_manhattan/NOT_sum_fall_outlier.jpeg",width = 10, height = 3)


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




















