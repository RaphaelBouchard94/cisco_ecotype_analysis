## Load libraries
library(data.table)
library(tidyverse)
library(ggpubr)



ld_chr36 <- fread("data/chr34/all_ind_chr34_outlier.ld.gz")

colnames(ld_chr36) <- c("site1","site2","dist","r2_ExpG","D","Dp","r2")

q_all <- quantile(ld_chr36$r2,0.99,na.rm=T)

# Keep valid R2
ld_chr36 <- ld_chr36[is.finite(r2) & r2 <= 1]

# Keep only high LD pairs
LD_all_hi <- ld_chr36[r2 > q_all, .(site1, site2, r2)]

# Count n_before (pairs ending at pos) and n_after (pairs starting at pos)
if (nrow(LD_all_hi) > 0) {
  before_all <- LD_all_hi[site1 < site2, .N, by = .(site2)]
  setnames(before_all, c("site2", "n_before"))
  after_all  <- LD_all_hi[site1 < site2, .N, by = .(site1)]
  setnames(after_all, c("site1", "n_after"))
} else {
  before_all <- data.table(site2 = integer(), n_before = integer())
  after_all  <- data.table(site1 = integer(), n_after = integer())
}

# Full set of positions
pos_all <- sort(unique(c(ld_chr36$site1, ld_chr36$site2)))
dd2_all <- data.table(pos = pos_all)

# Merge counts, replace NA by 0
dd2_all <- merge(dd2_all,
                 before_all,
                 by.x = "pos",
                 by.y = "site2",
                 all.x = TRUE)
dd2_all <- merge(dd2_all,
                 after_all,
                 by.x = "pos",
                 by.y = "site1",
                 all.x = TRUE)
dd2_all[is.na(n_before), n_before := 0]
dd2_all[is.na(n_after), n_after  := 0]

#outlier window
targ_win <- c(34800000,35619999)

dd2_all %<>% separate(col = pos,sep = ":",into = c("chr","pos"))

# Plot
ggplot(dd2_all, aes(x = as.numeric(pos))) +
  geom_point(
    aes(y = n_before),
    col = "grey70",
    alpha = 0.5,
    size = 0.8
  ) +
  geom_point(aes(y = n_after),
             col = "black",
             alpha = 0.5,
             size = 0.8) +
  annotate("rect",
           xmin = targ_win[1], xmax = targ_win[2],
           ymin = 0, ymax = 600,
           fill = "indianred",
           alpha = 0.2) +
  theme_bw(base_size = 16) +
  ylab("# linked SNPs") +
  xlab("Position (bp)") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_blank()
  )

ggsave("figures/05_ldxplot_outlier_window_chr36.png",
       width = 10, height = 10, unit = "cm",dpi = 1200)

