
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

win20kb <- fread("data/genome_scan/notsummer_rup_20kb.outlier.fst")

nsnp_win20kb <- win20kb %>% group_by(chrom,window) %>% count() %>%  ungroup()

window2remove <- nsnp_win20kb %>% filter(n <= 10) %>% dplyr::select(chrom,window)

win20kb_filt <- anti_join(win20kb,window2remove,by = c("chrom","window"))

thresh <- quantile(win20kb_filt$fst, probs = 0.999, na.rm = T)

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
  qbinom(0.999, j, p)
})

binom_dist <- data.frame(outlier_null_dist = binom_dist_outliers, nsnp = 1:max(win20kb_outlier_per_win$nsnp))

# Identify top-candidate windows
top_candidate <- win20kb_outlier_per_win %>%
  left_join(binom_dist, by = "nsnp") %>%
  dplyr::mutate(top_candidate = if_else(outlier <= outlier_null_dist, "non_top_candidate", "top_candidate"))

top_candidate %>% ungroup() %>% filter(top_candidate == "top_candidate") %>% count()

#500 windows enriched in outliers

wind20kb_top_candidate <- win20kb_filt %>% left_join(top_candidate %>% dplyr::select(-background,-outlier,-nsnp, -outlier_null_dist),
                                      by = c("chrom","window"))

fwrite(wind20kb_top_candidate,"data/genome_scan/notsum_rup_wind20kb_top_candidate.tsv",quote = F, sep = "\t")

top_windows_plot <- 
  wind20kb_top_candidate %>% 
  filter(top_candidate == "top_candidate") %>% 
  dplyr::select(-pos,-fst,-window,-status) %>% 
  group_by(chrom, window_start, window_end) %>% 
  dplyr::slice(1)

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
  current_count <- 1  # Initialize window count
  
  for (i in 2:nrow(df)) {
    if (df$chrom[i] == current_chrom && df$window_start[i] - current_window_end <= max_gap) {
      # Extend the current region
      current_window_end <- max(current_window_end, df$window_end[i])
      current_count <- current_count + 1  # Increment count
    } else {
      # Save the previous region and start a new one
      merged <- rbind(merged, data.frame(
        chrom = current_chrom, 
        window_start = current_window_start, 
        window_end = current_window_end, 
        status = current_status,
        n_windows_merged = current_count  # Add count
      ))
      current_chrom <- df$chrom[i]
      current_window_start <- df$window_start[i]
      current_window_end <- df$window_end[i]
      current_status <- df$top_candidate[i]
      current_count <- 1  # Reset count
    }
  }
  
  # Save the last region
  merged <- rbind(merged, data.frame(
    chrom = current_chrom, 
    window_start = current_window_start, 
    window_end = current_window_end, 
    status = current_status,
    n_windows_merged = current_count  # Add count
  ))
  
  return(merged)
}

# Apply the function
merged_data <- merge_windows(top_windows_plot, max_gap = 100000)

merged_data %<>% mutate(window_size = window_end - window_start)


#Stats on merge windows
max(merged_data$window_size)
median(merged_data$window_size)

#Stats on largest outlier window

merged_data %>% filter(window_size == max(merged_data$window_size))

win20kb %>% filter(chrom == "CM078251.1" & pos > 34800000 & pos < 35619999) %>% summarize(fst_mean = mean(fst))

n_snp_chr34 <- win20kb %>% filter(chrom == "CM078251.1" & pos > 34800000 & pos < 35619999) %>% count() %>% pull()

n_snp_chr34_0.9fst <- win20kb %>% filter(chrom == "CM078251.1" & pos > 34800000 & pos < 35619999 & fst > 0.9) %>% count() %>% pull()

n_snp_chr34_0.9fst/n_snp_chr34

fwrite(merged_data,"data/genome_scan/notsum_rup_merged_top_candidate.tsv",quote = F, sep = "\t")


#Save in Supplementary table:

wb <- loadWorkbook("res/supplementary_table.xlsx")

addWorksheet(wb, "Fst_outlier_window")
writeData(wb, sheet = "Fst_outlier_window", merged_data)

saveWorkbook(wb, "res/supplementary_table.xlsx", overwrite = TRUE)


#############################

#Plot manhattan with outlier regions

notsum_rup <- fread("data/genome_scan/notsummer_rup_genome_wide.25kb_windows.fst")[,-1]

colnames(notsum_rup) <- c("chrom","pos","nsites","fst")

axis_man <- notsum_rup %>% 
  group_by(chrom) %>%
  summarize(center = mean(pos))

# First, create a lookup table that maps your current chromosome names to numbers
# Make sure chromosomes are ordered as they appear in your data
chrom_order <- unique(notsum_rup$chrom)  # or use sort() if you want alphabetical order

chrom_lookup <- setNames(
  sprintf("%02d", 1:length(chrom_order)),  # Creates 01, 02, 03, ..., 38
  chrom_order
)

# Apply this mapping to your data
notsum_rup <- notsum_rup %>%
  mutate(chrom_num = chrom_lookup[chrom])

# Also apply to merged_data
merged_data <- merged_data %>%
  mutate(chrom_num = chrom_lookup[chrom])

# Now create the plot using chrom_num instead of chrom
fst_outlier_plot <- ggplot(notsum_rup, aes(x = pos , y = fst, color = as_factor(chrom_num))) +
  geom_point(alpha = 0.2, size=0.5) +
  geom_smooth(color = "black",linewidth = .5)+
  geom_rect(data=merged_data, 
            inherit.aes=FALSE, 
            aes(xmin=window_start, xmax=window_end,
                ymin=0,ymax=1, 
                group=chrom_num), 
            fill="darkred")+
  scale_color_manual(values = rep(
    c("grey40", "grey60"),
    19)) + 
  facet_grid(~chrom_num, scales="free_x", switch = "x") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  labs(y = "Fst",
       x = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.spacing = unit(0,'lines'),
        panel.border = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.title.y =  element_text(angle = 90, vjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(angle = 0,size = 8, vjust = 1, margin = margin(t = 2, b = 0)),
        strip.placement = "outside",
        strip.background = element_blank())


ggsave("figures/04_genome_scan_fst_outlier_window.png",
       fst_outlier_plot,height =5,width = 18, unit = "cm", dpi = 1000)
