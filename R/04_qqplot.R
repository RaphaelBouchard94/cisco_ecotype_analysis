#!/usr/bin/env Rscript
# ============================================================
# QQ plot of window median p-values (summer and fall runs)
#
# After Montejo-Kovacevich et al. 2021, Fig. S9, with one change:
# the expected distribution comes from the 200 permutations,
# not from Uniform(0, 1). The median of 20 (or 50) p-values is
# not uniform under the null, so the uniform line would be wrong.
#
# For each run:
#   observed  = sorted -log10(window median p) of the real scan
#   expected  = for each rank, the median across permutations of
#               the window median p at that rank
#   band      = 2.5-97.5% across permutations at that rank
#   lambda    = SNP-level genomic inflation factor (as estlambda,
#               method "median") and a window-level version
#               measured against the permutation null
#
# Input per run:
#   winmed: EITHER the folder holding winmed_<CHR>.rds from
#           perm_windows_chr.R, OR one genome-wide .rds made from them
#        (rows = windows; col 1 = observed, cols 2.. = permutations)
#   <lrt0_file>  observed ANGSD -doAsso output, all chromosomes
#
# Output:
#   res/qq_window_median_p.pdf and .png
#   res/qq_window_median_p_lambda.tsv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
  library(ggplot2)
})

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------
runs <- data.frame(
  panel      = c("a", "b"),
  label      = c("Summer run", "Fall run"),
  winmed     = c("data/gwas_length/summer_genome_winmed.rds",            # CHECK these paths
                 "data/gwas_length/fall_genome_winmed.rds"),             # folder or genome .rds
  lrt0_file  = c("data/gwas_length/summergwas_ALL_CHR.lrt0.gz",
                 "data/gwas_length/fallgwas_ALL_CHR.lrt0.gz")
)
n_top      <- 2000             # smallest-p windows: every one is plotted
n_bulk     <- 3000             # points kept for the rest of the distribution
band_probs <- c(0.025, 0.975)  # permutation band
out_prefix <- "res/qq_window_median_p"
dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Functions
# ------------------------------------------------------------

# p-value -> 1-df chi-square statistic (what estlambda does internally)
to_chisq <- function(p) qchisq(p, df = 1, lower.tail = FALSE)

# Genome-wide matrix: read one combined .rds, or stack the per-chromosome ones
read_winmed <- function(path) {
  if (grepl("\\.rds$", path)) {
    M <- readRDS(path)
  } else {
    files <- list.files(path, pattern = "^winmed_.*\\.rds$", full.names = TRUE)
    if (length(files) == 0) stop("No winmed_<CHR>.rds files in ", path)
    M <- do.call(rbind, lapply(files, readRDS))
  }
  bad <- rowAnyNAs(M)                      # windows with no tested SNP in some run
  if (any(bad)) message(sprintf("  dropping %d windows with NA medians", sum(bad)))
  M[!bad, , drop = FALSE]
}

# Which ranks to plot: all of the low-p tail, then log-spaced ranks in the bulk
rank_grid <- function(N, n_top, n_bulk) {
  top <- seq_len(min(N, n_top))
  if (N <= n_top) return(top)
  bulk <- round(exp(seq(log(n_top + 1), log(N), length.out = n_bulk)))
  unique(c(top, bulk))
}

# Observed and expected order statistics
qq_table <- function(M, grid) {
  obs <- sort(M[, 1])[grid]
  # one column per permutation: its window medians sorted, kept at the grid ranks
  null_sorted <- vapply(2:ncol(M), function(k) sort(M[, k])[grid],
                        numeric(length(grid)))
  nl <- -log10(null_sorted)
  data.table(rank     = grid,
             observed = -log10(obs),
             expected = rowMedians(nl),
             lo       = rowQuantiles(nl, probs = band_probs[1]),
             hi       = rowQuantiles(nl, probs = band_probs[2]))
}

# Window-level lambda against the permutation null
lambda_window <- function(M) {
  perm_meds <- colMedians(M[, -1, drop = FALSE])      # genome-wide median of each permutation
  ref       <- to_chisq(median(perm_meds))            # typical null value
  lam_perm  <- to_chisq(perm_meds) / ref
  c(lambda_window = to_chisq(median(M[, 1])) / ref,
    perm_lo = unname(quantile(lam_perm, band_probs[1])),
    perm_hi = unname(quantile(lam_perm, band_probs[2])))
}

# SNP-level lambda from the observed ANGSD output
lambda_snp <- function(file) {
  if (is.na(file) || !file.exists(file)) {
    message("  no lrt0 file found; SNP-level lambda skipped")
    return(c(lambda_snp = NA_real_, n_snps = NA_real_))
  }
  cmd <- paste("zcat", shQuote(file))
  col <- intersect(c("LRT", "LRTscore"), names(fread(cmd = cmd, nrows = 0)))[1]
  x   <- fread(cmd = cmd, select = col)[[1]]
  x   <- x[x != -999 & is.finite(x)]
  c(lambda_snp = median(x) / qchisq(0.5, df = 1), n_snps = length(x))
}

# ------------------------------------------------------------
# Compute
# ------------------------------------------------------------
qq_list  <- list()
lam_list <- list()

for (i in seq_len(nrow(runs))) {
  message(runs$label[i])
  M    <- read_winmed(runs$winmed[i])
  grid <- rank_grid(nrow(M), n_top, n_bulk)
  
  qq_list[[i]] <- qq_table(M, grid)[, `:=`(panel = runs$panel[i], label = runs$label[i])]
  
  lam_list[[i]] <- data.table(panel = runs$panel[i], label = runs$label[i],
                              n_windows = nrow(M), n_perm = ncol(M) - 1L,
                              t(lambda_snp(runs$lrt0_file[i])),
                              t(lambda_window(M)))
  rm(M); invisible(gc())
}

qq  <- rbindlist(qq_list)
lam <- rbindlist(lam_list)
# panel order follows `runs`, not the alphabet
qq[,  label := factor(label, levels = runs$label)]
lam[, label := factor(label, levels = runs$label)]
fwrite(lam, paste0(out_prefix, "_lambda.tsv"), sep = "\t")
print(lam, digits = 3)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
# lambda labels as plotmath strings (drawn with parse = TRUE)
lam[, txt_snp := ifelse(is.na(lambda_snp), NA_character_,
                        sprintf('lambda[SNP] == "%.3f"', lambda_snp))]
lam[, txt_win := sprintf('lambda[window] == "%.3f"~~"(null: %.3f-%.3f)"',
                         lambda_window, perm_lo, perm_hi)]
lam_long <- melt(lam[, .(label, txt_snp, txt_win)], id.vars = "label",
                 value.name = "txt", na.rm = TRUE)
lam_long[, line := ifelse(variable == "txt_snp", 2.6, 1.2)]   # vertical position in lines

# same range on both axes within each panel, so the y = x line is a true diagonal
lims <- qq[, .(m = max(observed, hi, expected)), by = label]

p <- ggplot(qq, aes(expected, observed)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey85") +
  geom_abline(intercept = 0, slope = 1, colour = "red3", linewidth = 0.4) +
  geom_point(size = 0.5, colour = "black") +
  geom_blank(data = lims, aes(x = m, y = m)) +
  geom_blank(data = lims, aes(x = 0, y = 0)) +
  geom_text(data = unique(qq[, .(label, panel)]), aes(label = panel),
            x = -Inf, y = Inf, hjust = -0.6, vjust = 1.5,
            fontface = "bold", size = 7 / .pt * 1.2, inherit.aes = FALSE) +
  geom_text(data = lam_long, aes(label = txt, vjust = -line + 0.1),
            x = Inf, y = -Inf, hjust = 1.05, parse = TRUE,
            size = 6 / .pt, inherit.aes = FALSE) +
  facet_wrap(~ label, scales = "free") +
  labs(x = expression(Expected~-log[10](italic(P))~"(permutations)"),
       y = expression(Observed~-log[10](italic(P)))) +
  theme_classic(base_size = 7) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(size = 7, hjust = 0),
        aspect.ratio     = 1)


ggsave("figures/04_gwasqqplot.png", p, width = 183, height = 95, units = "mm", dpi = 300)



