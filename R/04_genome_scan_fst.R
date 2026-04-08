#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(magrittr)
library(data.table)
library(patchwork)
library(GenomicRanges)

#------------#
#--- DATA ---#
#------------#

#Sliding window

RUP_NOT_summer <- fread("data/genome_scan/notsummer_rup_genome_wide.25kb_windows.fst")[,-1]

colnames(RUP_NOT_summer) <- c("Chromosome","Position","nsites","fst")


chr_offsets <- RUP_NOT_summer %>%
  group_by(Chromosome) %>%
  summarize(chr_len = max(Position)) %>%
  mutate(
    pos_add = lag(cumsum(as.numeric(chr_len)), default = 0)
  ) %>%
  dplyr::select(Chromosome, pos_add)

RUP_NOT_summer <- RUP_NOT_summer %>%
  left_join(chr_offsets, by = "Chromosome") %>%
  mutate(pos_cum = Position + pos_add)


# Function to add cumulative positions
add_cumulative_pos <- function(df) {
  df %>%
    mutate(Chromosome = factor(Chromosome, levels = unique(Chromosome))) %>%
    group_by(Chromosome) %>%
    summarize(max_pos = max(Position)) %>%
    mutate(pos_add = lag(cumsum(as.numeric(max_pos)), default = 0)) %>%
    select(Chromosome, pos_add) %>%
    left_join(df, ., by = "Chromosome") %>%
    mutate(pos_cum = Position + pos_add)
}


#----------------#
#--- Analysis ---#
#----------------#

# Calculate axis positions for continuous plot
axis_man <- RUP_NOT_summer %>%
  group_by(Chromosome) %>%
  summarize(center = mean(pos_cum))

##########################
#Chromosome offset
#add outlier window
outlier_windows <- fread("data/genome_scan/notsum_rup_merged_top_candidate.tsv")

outlier_windows <- outlier_windows %>%
  left_join(chr_offsets, by = c("chrom" = "Chromosome")) %>%
  mutate(
    start_cum = as.numeric(window_start) + as.numeric(pos_add),
    end_cum   = as.numeric(window_end)   + as.numeric(pos_add)
  )

#--------------------------------------------------------#
#--- METHOD 3: Equal-Spaced Continuous Plot ------------#
#--------------------------------------------------------#

# Give each chromosome equal space in continuous plot
# regardless of actual length

# Calculate equal spacing
n_chroms <- length(unique(RUP_NOT_summer$Chromosome))
chrom_spacing <- 1e7  # 10 Mb spacing between chromosomes

RUP_NOT_summer_equal <- RUP_NOT_summer %>%
  mutate(Chromosome = factor(Chromosome, levels = unique(Chromosome))) %>%
  group_by(Chromosome) %>%
  mutate(
    # Normalize within chromosome
    pos_norm = (Position - min(Position)) / (max(Position) - min(Position)),
    chrom_min = min(Position),
    chrom_max = max(Position),
    chrom_range = chrom_max - chrom_min
  ) %>%
  ungroup() %>%
  mutate(
    # Add equal spacing
    chrom_num = as.numeric(Chromosome),
    pos_equal = pos_norm * chrom_spacing + (chrom_num - 1) * chrom_spacing * 1.2
  )

chrom_info <- RUP_NOT_summer_equal %>%
  group_by(Chromosome) %>%
  summarize(
    chrom_min = dplyr::first(chrom_min),
    chrom_max = dplyr::first(chrom_max),
    chrom_range = dplyr::first(chrom_range),
    chrom_num = dplyr::first(chrom_num)
  )


outlier_windows_equal <- outlier_windows %>%
  dplyr::rename(Chromosome = chrom,
                start = window_start,
                end = window_end) %>%
  left_join(chrom_info, by = "Chromosome") %>%
  mutate(
    # Normalize positions within chromosome
    start_norm = (start - chrom_min) / chrom_range,
    end_norm = (end - chrom_min) / chrom_range,
    # Transform to equal spacing
    start_equal = start_norm * chrom_spacing + (chrom_num - 1) * chrom_spacing * 1.2,
    end_equal = end_norm * chrom_spacing + (chrom_num - 1) * chrom_spacing * 1.2,
    mid_equal = (start_equal + end_equal) / 2
  )

# Calculate axis positions
axis_equal <- RUP_NOT_summer_equal %>%
  group_by(Chromosome) %>%
  summarize(center = mean(pos_equal))

p_equal_continuous <- ggplot(RUP_NOT_summer_equal, 
                             aes(x = pos_equal, y = fst, color = Chromosome)) +
  geom_point(alpha = 0.75, size = 0.8) +
  geom_segment(data = outlier_windows_equal,
               aes(x = start_equal, xend = end_equal, 
                   y = -0.03, yend = -0.03),
               color = "red", linewidth = 6, lineend = "butt") +
  # Scales
  scale_x_continuous(
    label = seq(1,38),
    breaks = axis_equal$center,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(-0.05, max(RUP_NOT_summer_equal$fst, na.rm = TRUE) * 1.05)
  ) +
  scale_color_manual(
    values = rep(c("grey50", "grey70"), 
                 length.out = length(unique(RUP_NOT_summer_normalized$Chromosome)))
  ) +
  # Labels
  labs(
    x = "Chromosome",
    y = expression(F[ST]),
  ) +
  
  # Theme
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid= element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

ggsave("figures/04_genome_fstscan_seg_notsum_rup.png", p_equal_continuous,width = 30, height = 10, unit = "cm", dpi = 900)

##############################################

outlier_windows_all <- fread("data/genome_scan/notsum_rup_wind20kb_top_candidate.tsv")

outlier_windows_all %<>% filter(top_candidate == "top_candidate") %>% 
  group_by(chrom, window_start, window_end) %>% 
  dplyr::slice(1)

# Count outlier windows per chromosome
outlier_counts <- outlier_windows_all %>%
  group_by(chrom) %>%
  summarize(outlier_windows = n())

outlier_counts %<>% mutate(windows_prop = outlier_windows/sum(outlier_windows))


# Create barplot
p_barplot <- ggplot(outlier_counts, aes(x = chrom, y = windows_prop)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_x_discrete(
    expand = c(0.01, 0)
  ) +
  
  # Labels
  labs(
    x = NULL,  # Remove x-axis label (will be in bottom plot)
    y = "Outlier Windows (%)",
  ) +
  
  # Theme
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid = element_blank()
  )

print(p_barplot)




