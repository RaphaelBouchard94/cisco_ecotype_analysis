library(tidyverse)
library(ggrepel)

window_df <- tibble(
  chrom = "chr36",
  start = 34800000,
  end   = 35619999,
)

gene_df <- tibble(
  gene  = c("mc4r", "cad-20"),
  start = c(35060636,35589506),
  end   = c(35059641, 35618151)
)


gene_df <- gene_df %>%
  mutate(strand = c("+", "+"))

ggplot() +
  # outlier window
  geom_segment(
    data = window_df,
    aes(x = start, xend = end, y = 0, yend = 0),
    linewidth = 8,
    color = "indianred", alpha = 0.4
  ) +
  
  # genes
  geom_segment(
    data = gene_df,
    aes(x = start, xend = end, y = 0, yend = 0, color = gene),
    linewidth = 8
  ) +
  scale_color_manual(values = c("skyblue4","skyblue4"))+
  scale_x_continuous(
    name = "Position (Mb)",
    labels = function(x) x / 1e6
  ) +
  scale_y_continuous(NULL, breaks = NULL) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )
  
ggsave("figures/05_gene_position_outlier.png",height = 2, width = 8, unit = "cm", dpi = 900)
