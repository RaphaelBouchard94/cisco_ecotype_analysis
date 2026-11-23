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


notsum_div <- fread("data/genome_scan/notsummer.thetas")
notfall_div <- fread("data/genome_scan/notfall.thetas")
rup_div <- fread("data/genome_scan/rup.thetas")


pi <- bind_rows(
  notsummer = notsum_div,
  rup       = rup_div,
  notfall = notfall_div,
  .id = "pop"
) %>%
  rename(chr = `#Chromo`) %>%
  mutate(pi = exp(Pairwise))

win_size <- 50000

pi_win <- pi %>%
  mutate(win = Pos %/% win_size * win_size + win_size / 2) %>%
  group_by(pop, chr, win) %>%
  summarise(pi = mean(pi), n_sites = n(), .groups = "drop") 

ggplot(pi_win, aes(pop, pi, fill = pop)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  labs(x = NULL, y = "Nucleotide diversity (π)") +
  theme_bw() +
  theme(legend.position = "none")

ggplot(pi_win, aes(win / 1e6, pi, colour = pop)) +
  geom_line(alpha = 0.7) +
  facet_grid(~ chr, scales = "free_x", space = "free_x") +
  labs(x = "Position (Mb)", y = "Nucleotide diversity (π)", colour = NULL) +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"))

pi_win %>% group_by(pop) %>% summarize(mean_pi = mean(pi, na.rm = T))



block_size <- 1e6

blocks <- pi %>%
  mutate(block = paste(chr, Pos %/% block_size)) %>%
  group_by(pop, block) %>%
  summarise(sum_pi = sum(pi), n = n(), .groups = "drop")

boot <- replicate(1000, {
  b <- tibble(block = sample(unique(blocks$block), replace = TRUE))
  blocks %>%
    inner_join(b, by = "block", relationship = "many-to-many") %>%
    group_by(pop) %>%
    summarise(pi = sum(sum_pi) / sum(n)) %>%
    deframe()
})

apply(boot, 1, quantile, c(0.025, 0.975))                     # CI per pop
quantile(boot["notsummer", ] / boot["notfall", ], c(0.025, 0.975))
                                                      