#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(data.table)
library(GenomicRanges)

#------------#
#--- DATA ---#
#------------#

# Sliding windows
RUP_NOT_summer <- fread("data/genome_scan/notsummer_rup_genome_wide.25kb_windows.fst")[, -1]
colnames(RUP_NOT_summer) <- c("Chromosome", "Position", "nsites", "fst")

# Outlier windows
outlier_windows <- fread("data/genome_scan/notsum_rup_merged_top_candidate.tsv")

#------------------------------#
#--- Identify outlier SNPs ----#
#------------------------------#

# Sliding windows as GRanges
gr_windows <- GRanges(
  seqnames = RUP_NOT_summer$Chromosome,
  ranges   = IRanges(start = RUP_NOT_summer$Position,
                     end   = RUP_NOT_summer$Position)
)

# Outlier windows as GRanges
gr_outliers <- GRanges(
  seqnames = outlier_windows$chrom,
  ranges   = IRanges(start = outlier_windows$window_start,
                     end   = outlier_windows$window_end)
)

# Overlaps
hits <- findOverlaps(gr_windows, gr_outliers)

RUP_NOT_summer$outlier <- FALSE
RUP_NOT_summer$outlier[queryHits(hits)] <- TRUE

#------------------------------------------#
#--- Equal-spaced chromosome coordinates --#
#------------------------------------------#

chrom_spacing <- 1e7

RUP_NOT_summer_equal <- RUP_NOT_summer %>%
  mutate(Chromosome = factor(Chromosome, levels = unique(Chromosome))) %>%
  group_by(Chromosome) %>%
  mutate(
    pos_norm = (Position - min(Position)) / (max(Position) - min(Position)),
    chrom_num = as.numeric(Chromosome)
  ) %>%
  ungroup() %>%
  mutate(
    pos_equal = pos_norm * chrom_spacing +
      (chrom_num - 1) * chrom_spacing * 1.2
  )

axis_equal <- RUP_NOT_summer_equal %>%
  group_by(Chromosome) %>%
  summarize(center = mean(pos_equal))

#--------------------------------#
#--- Prepare outlier segments ---#
#--------------------------------#

chrom_info <- RUP_NOT_summer_equal %>%
  group_by(Chromosome) %>%
  summarize(
    chrom_min = min(Position),
    chrom_max = max(Position),
    chrom_range = chrom_max - chrom_min,
    chrom_num = dplyr::first(chrom_num)
  )

outlier_windows_equal <- outlier_windows %>%
  dplyr::rename(Chromosome = chrom,
         start = window_start,
         end   = window_end) %>%
  left_join(chrom_info, by = "Chromosome") %>%
  mutate(
    start_norm = (start - chrom_min) / chrom_range,
    end_norm   = (end   - chrom_min) / chrom_range,
    start_equal = start_norm * chrom_spacing +
      (chrom_num - 1) * chrom_spacing * 1.2,
    end_equal   = end_norm   * chrom_spacing +
      (chrom_num - 1) * chrom_spacing * 1.2
  )

#----------------#
#--- Plot ------#
#----------------#

bg <- RUP_NOT_summer_equal %>% filter(!outlier)
ol <- RUP_NOT_summer_equal %>% filter(outlier)

p <- ggplot() +
  geom_point(
    data = bg,
    aes(pos_equal, fst, color = Chromosome),
    size = 0.4,
    alpha = 0.6
  ) +
  geom_point(
    data = ol,
    aes(pos_equal, fst),
    color = "red4",
    size = 0.4
  ) +
  scale_x_continuous(
    breaks = axis_equal$center,
    labels = seq_len(n_distinct(RUP_NOT_summer_equal$Chromosome)),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    limits = c(0, max(RUP_NOT_summer_equal$fst, na.rm = TRUE) * 1.05),
    expand = c(0, 0)
  ) +
  scale_color_manual(
    values = rep(c("grey50", "grey70"),
                 length.out = n_distinct(bg$Chromosome))
  ) +
  labs(x = "Chromosome", y = expression(F[ST])) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 10),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )


##############################################

all_chromosomes <- sort(unique(RUP_NOT_summer$Chromosome))

chrom_info <- RUP_NOT_summer %>%
  group_by(Chromosome) %>%
  summarize(
    chrom_length = max(Position),
    n_windows = n(),
    .groups = "drop"
  ) %>%
  dplyr::rename(chrom = Chromosome)


outlier_windows_all <- fread("data/genome_scan/notsum_rup_wind20kb_top_candidate.tsv")

outlier_counts <- outlier_windows_all %>%
  filter(top_candidate == "top_candidate") %>%
  group_by(chrom) %>%
  summarize(outlier_windows = n(), .groups = "drop") %>%
  right_join(
    tibble(chrom = all_chromosomes),
    by = "chrom"
  ) %>%
  mutate(outlier_windows = replace_na(outlier_windows, 0)) %>%
  left_join(chrom_info, by = "chrom") %>%
  mutate(
    outlier_rate = outlier_windows / (chrom_length/1e6),
  )


p_barplot <- ggplot(outlier_counts, aes(x = chrom, y = outlier_rate)) +
  geom_bar(stat = "identity", width = 0.8, fill = "grey40") +
  scale_x_discrete(
    labels = seq(1,38),
    expand = c(0.01, 0)
  ) +
  labs(
    x = NULL,
    y = "Outlier windows\n per Mb"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

pplot <- p_barplot + p + plot_layout(ncol = 1, nrow = 2, height = c(0.3,0.7))

ggsave(
  "figures/04_genome_fstscan_outlier_snps_barplot.png",
  pplot,
  width = 30,
  height = 10,
  units = "cm",
  dpi = 900
)
