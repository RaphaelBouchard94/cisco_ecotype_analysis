# ============================================================
# Length-at-maturity GWAS: overlap between summer-run and fall-run
# significant windows, and with summer-fall Fst outlier regions
#
# 1. Pairwise overlaps (bp and windows):
#      fall x summer, fall x Fst, summer x Fst, (fall & summer) x Fst
#    each tested against a null built by circularly shifting the second
#    set of regions along each chromosome (keeps their number, sizes
#    and clustering)
# 2. Table of regions significant in both runs, with peak statistics,
#    Fst outlier status and genes
# 3. Direction of effect: is the length-increasing allele the same in
#    both runs? (one SNP per shared region, sign test)
# 4. Miami plot: summer up, fall down, Fst outlier regions shaded
#
# Inputs are the outputs of the two 04_gwas_*_length_perm_windows.R
# scripts. Windows are compared by genomic coordinates, since the
# 50-SNP windows of the two runs are built from different SNP sets.
# ============================================================

library(tidyverse)
library(data.table)
library(GenomicRanges)

set.seed(2026)

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------
summer_gwas_file  <- "data/gwas_length/summergwas_ALL_CHR.lrt0.gz"
fall_gwas_file    <- "data/gwas_length/fallgwas_ALL_CHR.lrt0.gz"
summer_win_file   <- "res/gwas_summer_length_sig_windows.tsv"
fall_win_file     <- "res/gwas_fall_length_sig_windows.tsv"
summer_gene_file  <- "res/gwas_summer_length_sig_window_genes.tsv"
fall_gene_file    <- "res/gwas_fall_length_sig_window_genes.tsv"
fst_out_file      <- "data/genome_scan/notsum_rup_merged_top_candidate.tsv"
fst_scan_file     <- "data/genome_scan/notsummer_rup_genome_wide.25kb_windows.fst"
n_shift           <- 1000

dir.create("res", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ============================================================
# Load data
# ============================================================

read_gwas <- function(f) {
  fread(f) %>%
    mutate(LRT = suppressWarnings(as.numeric(LRT))) %>%
    filter(!is.na(LRT), LRT >= 0) %>%
    mutate(P = pchisq(LRT, df = 1, lower.tail = FALSE),
           neglog10p = -log10(P),
           Chromosome = as.character(Chromosome))
}

gwas_s <- read_gwas(summer_gwas_file)
gwas_f <- read_gwas(fall_gwas_file)

win_s <- fread(summer_win_file) %>% mutate(Chromosome = as.character(Chromosome))
win_f <- fread(fall_win_file)   %>% mutate(Chromosome = as.character(Chromosome))

fst_out <- fread(fst_out_file) %>%
  transmute(Chromosome = as.character(chrom), start = window_start, end = window_end)

fst_scan <- fread(fst_scan_file)[, -1]
colnames(fst_scan) <- c("Chromosome", "Position", "nsites", "fst")
fst_scan$Chromosome <- as.character(fst_scan$Chromosome)

# Chromosomes present in both GWAS and the Fst scan
chroms <- Reduce(intersect, list(unique(gwas_s$Chromosome),
                                 unique(gwas_f$Chromosome),
                                 unique(fst_scan$Chromosome)))
if (length(chroms) == 0) stop("No chromosome names shared by the GWAS and Fst files")

# Chromosome length = furthest position seen in any input
chr_len_v <- bind_rows(
  gwas_s %>% dplyr::select(Chromosome, Position),
  gwas_f %>% dplyr::select(Chromosome, Position),
  fst_scan %>% dplyr::select(Chromosome, Position),
  fst_out %>% transmute(Chromosome, Position = end),
  win_s %>% transmute(Chromosome, Position = end_pos),
  win_f %>% transmute(Chromosome, Position = end_pos)
) %>%
  filter(Chromosome %in% chroms) %>%
  group_by(Chromosome) %>%
  summarize(len = max(Position), .groups = "drop") %>%
  { setNames(.$len, .$Chromosome) }

to_gr <- function(chr, s, e, ...) {
  gr <- GRanges(chr, IRanges(s, e), ...)
  gr[as.character(seqnames(gr)) %in% chroms]
}

win_s_gr <- to_gr(win_s$Chromosome, win_s$start_pos, win_s$end_pos, window_id = win_s$window_id)
win_f_gr <- to_gr(win_f$Chromosome, win_f$start_pos, win_f$end_pos, window_id = win_f$window_id)

reg_s   <- reduce(win_s_gr)                 # bp covered by significant summer windows
reg_f   <- reduce(win_f_gr)                 # bp covered by significant fall windows
reg_fst <- reduce(to_gr(fst_out$Chromosome, fst_out$start, fst_out$end))
reg_sf  <- GenomicRanges::intersect(reg_s, reg_f, ignore.strand = TRUE)  # significant in both

bp <- function(gr) sum(as.numeric(width(gr)))
genome_bp <- sum(as.numeric(chr_len_v))

cat(sprintf("Genome analysed: %.1f Mb on %d chromosomes\n", genome_bp / 1e6, length(chroms)))
cat(sprintf("Summer: %d significant windows covering %.2f Mb (%.2f%%)\n",
            length(win_s_gr), bp(reg_s) / 1e6, 100 * bp(reg_s) / genome_bp))
cat(sprintf("Fall:   %d significant windows covering %.2f Mb (%.2f%%)\n",
            length(win_f_gr), bp(reg_f) / 1e6, 100 * bp(reg_f) / genome_bp))
cat(sprintf("Fst outliers: %d regions covering %.2f Mb (%.2f%%)\n\n",
            length(reg_fst), bp(reg_fst) / 1e6, 100 * bp(reg_fst) / genome_bp))



GenomicRanges::intersect(reg_fst, reg_s, ignore.strand = T)
GenomicRanges::intersect(reg_fst, reg_f, ignore.strand = T)
GenomicRanges::intersect(reg_s, reg_f, ignore.strand = T)

# ============================================================
# 1. Overlap tests (circular shift null)
# ============================================================

circ_shift <- function(gr, len) {
  out <- lapply(split(gr, seqnames(gr)), function(x) {
    if (length(x) == 0) return(NULL)
    L  <- len[[as.character(seqnames(x))[1]]]
    sh <- sample.int(as.integer(L), 1)
    s  <- (start(x) - 1 + sh) %% L + 1
    e  <- s + width(x) - 1
    wrap <- e > L                                   # region crosses the end: split it
    pieces <- GRanges(seqnames(x), IRanges(s, pmin(e, L)))
    if (any(wrap))
      pieces <- c(pieces, GRanges(seqnames(x)[wrap], IRanges(1, e[wrap] - L)))
    pieces
  })
  out <- out[!vapply(out, is.null, logical(1))]
  reduce(do.call(c, unname(out)))
}

shared_bp <- function(a, b) bp(GenomicRanges::intersect(a, b, ignore.strand = TRUE))

# a = fixed set, b = set that is shifted; a_win = windows of a (for window counts)
overlap_test <- function(name, a, b, a_win = NULL) {
  obs_bp  <- shared_bp(a, b)
  null_bp <- replicate(n_shift, shared_bp(a, circ_shift(b, chr_len_v)))
  res <- tibble(
    comparison  = name,
    a_Mb        = bp(a) / 1e6,
    b_Mb        = bp(b) / 1e6,
    shared_kb   = obs_bp / 1e3,
    expected_kb = mean(null_bp) / 1e3,
    fold        = obs_bp / mean(null_bp),
    emp_p       = (1 + sum(null_bp >= obs_bp)) / (n_shift + 1),
    n_windows_a = if (is.null(a_win)) NA_integer_ else length(a_win),
    n_windows_a_overlapping_b = if (is.null(a_win)) NA_integer_ else sum(overlapsAny(a_win, b))
  )
  attr(res, "null") <- null_bp
  res
}

tests <- list(
  overlap_test("fall x summer",          reg_f,  reg_s,   win_f_gr),
  overlap_test("fall x Fst",             reg_f,  reg_fst, win_f_gr),
  overlap_test("summer x Fst",           reg_s,  reg_fst, win_s_gr),
  overlap_test("(fall & summer) x Fst",  reg_sf, reg_fst)
)

test_tab <- bind_rows(tests)

print(test_tab %>% mutate(across(where(is.double), ~ signif(.x, 3))), width = Inf)

fwrite(test_tab, "res/gwas_summer_fall_fst_overlap_tests.tsv", sep = "\t")

# ============================================================
# 2. Regions significant in both runs
# ============================================================

peak_in <- function(gwas, regions) {
  snp_gr <- GRanges(gwas$Chromosome, IRanges(gwas$Position, width = 1))
  h <- findOverlaps(regions, snp_gr)
  tibble(region = queryHits(h),
         Position = gwas$Position[subjectHits(h)],
         neglog10p = gwas$neglog10p[subjectHits(h)]) %>%
    group_by(region) %>%
    summarize(peak_pos = Position[which.max(neglog10p)],
              max_neglog10p = max(neglog10p), .groups = "drop")
}

genes_in <- function(gene_file, regions) {
  if (!file.exists(gene_file) || length(regions) == 0) return(tibble(region = integer(), gene = character()))
  g <- fread(gene_file) %>%
    distinct(Chromosome, gene_start, gene_end, gene_name, .keep_all = TRUE)
  lab <- if ("external_gene_name" %in% names(g))
    ifelse(is.na(g$external_gene_name) | g$external_gene_name == "", g$gene_name, g$external_gene_name)
  else g$gene_name
  g_gr <- GRanges(as.character(g$Chromosome), IRanges(g$gene_start, g$gene_end))
  h <- findOverlaps(regions, g_gr)
  tibble(region = queryHits(h), gene = lab[subjectHits(h)])
}

if (length(reg_sf) > 0) {
  pk_s <- peak_in(gwas_s, reg_sf) %>% dplyr::rename(summer_peak_pos = peak_pos, summer_max_neglog10p = max_neglog10p)
  pk_f <- peak_in(gwas_f, reg_sf) %>% dplyr::rename(fall_peak_pos = peak_pos, fall_max_neglog10p = max_neglog10p)
  
  genes <- bind_rows(genes_in(summer_gene_file, reg_sf), genes_in(fall_gene_file, reg_sf)) %>%
    distinct() %>%
    group_by(region) %>%
    summarize(genes = paste(sort(unique(gene)), collapse = ", "), .groups = "drop")
  
  fst_hit <- findOverlaps(reg_sf, reg_fst)
  
  shared_regions <- tibble(
    region     = seq_along(reg_sf),
    Chromosome = as.character(seqnames(reg_sf)),
    start      = start(reg_sf),
    end        = end(reg_sf),
    width_kb   = width(reg_sf) / 1e3,
    n_summer_windows = countOverlaps(reg_sf, win_s_gr),
    n_fall_windows   = countOverlaps(reg_sf, win_f_gr),
    fst_outlier      = overlapsAny(reg_sf, reg_fst)
  ) %>%
    left_join(pk_s, by = "region") %>%
    left_join(pk_f, by = "region") %>%
    left_join(genes, by = "region") %>%
    arrange(desc(fst_outlier), desc(pmin(summer_max_neglog10p, fall_max_neglog10p)))
  
  fwrite(shared_regions, "res/gwas_summer_fall_shared_regions.tsv", sep = "\t")
  cat(sprintf("\nRegions significant in both runs: %d (%.1f kb), %d of them in Fst outlier regions\n",
              nrow(shared_regions), sum(shared_regions$width_kb), sum(shared_regions$fst_outlier)))
} else {
  shared_regions <- tibble()
  cat("\nNo region is significant in both runs\n")
}


# ============================================================
# 3. Miami plot
# ============================================================

chrom_spacing <- 1e7
chrom_gap     <- 2e6

chr_order <- unique(c(gwas_s$Chromosome, gwas_f$Chromosome))
chr_order <- chr_order[chr_order %in% chroms]
chr_map <- tibble(Chromosome = chr_order,
                  chr_index  = seq_along(chr_order) - 1,
                  max_pos    = chr_len_v[chr_order])

to_equal <- function(chr, pos) {
  m <- chr_map[match(chr, chr_map$Chromosome), ]
  m$chr_index * (chrom_spacing + chrom_gap) + pmin(pos / m$max_pos, 1) * chrom_spacing
}

# Regions significant in both runs AND inside an Fst outlier region
reg_sf_fst <- GenomicRanges::intersect(reg_sf, reg_fst, ignore.strand = TRUE)

# Status of each SNP, from most to least specific. case_when() keeps the
# first condition that is TRUE, so a SNP in reg_sf_fst is "shared_fst"
# even though it is also inside reg_sf, own_fst and own.
#   own     = significant windows of this run
#   own_fst = significant windows of this run that lie in Fst outlier regions
flag <- function(gwas, own) {
  own_fst <- GenomicRanges::intersect(own, reg_fst, ignore.strand = TRUE)
  gwas <- gwas %>% filter(Chromosome %in% chroms)
  snp_gr <- GRanges(gwas$Chromosome, IRanges(gwas$Position, width = 1))
  gwas %>%
    mutate(pos_equal = to_equal(Chromosome, Position),
           chr_shade = factor(chr_map$chr_index[match(Chromosome, chr_map$Chromosome)] %% 2),
           status = case_when(overlapsAny(snp_gr, reg_sf_fst) ~ "shared_fst",
                              overlapsAny(snp_gr, reg_sf)     ~ "shared",
                              overlapsAny(snp_gr, own_fst)    ~ "run_fst",
                              overlapsAny(snp_gr, own)        ~ "run",
                              TRUE                            ~ "none"))
}

plot_df <- bind_rows(
  flag(gwas_s, reg_s) %>% mutate(run = "summer", y = neglog10p),
  flag(gwas_f, reg_f) %>% mutate(run = "fall",   y = -neglog10p)
)

print(count(plot_df, run, status))

fst_rect <- tibble(Chromosome = as.character(seqnames(reg_fst)),
                   start = start(reg_fst), end = end(reg_fst)) %>%
  mutate(xmin = to_equal(Chromosome, start),
         xmax = pmax(to_equal(Chromosome, end), xmin + 3e4))   # keep narrow regions visible

axis_df <- chr_map %>%
  mutate(center = chr_index * (chrom_spacing + chrom_gap) + chrom_spacing / 2)

ymax <- max(abs(plot_df$y)) * 1.05

col_run    <- "#ABC7D4FF"   # significant in one run
col_run_fst <- "#FF9E9DFF"  # significant in one run and inside an Fst outlier region
col_shared <- "#175149FF"   # significant in both runs
col_shared_fst <- "darkred"

# Legend labels (also used as the keys that link each layer to its colour)
lab_run     <- "Significant in one ecotype"
lab_shared  <- "Significant in both ecotypes"
lab_run_fst <- "Significant in one ecotype within Fst outlier window"
lab_shared_fst <- "Significant in both ecotypes within Fst outlier window"

miami <- ggplot() +
  geom_point(data = filter(plot_df, status == "none"),
             aes(pos_equal, y, color = chr_shade), size = 0.1, alpha = 0.3) +
  geom_point(data = filter(plot_df, status == "run"),
             aes(pos_equal, y, fill = lab_run),
             shape = 21, color = col_run, size = 0.1, alpha = 0.3) +
  geom_point(data = filter(plot_df, status == "shared"),
             aes(pos_equal, y, fill = lab_shared),
             shape = 21, color = col_shared, size = 0.1, alpha = 0.3) +
  geom_point(data = filter(plot_df, status == "run_fst"),
             aes(pos_equal, y, fill = lab_run_fst),
             shape = 21, color = col_run_fst, size = 0.1, alpha = 0.3) +
  geom_point(data = filter(plot_df, status == "shared_fst"),
             aes(pos_equal, y, fill = lab_shared_fst),
             shape = 21, color = col_shared_fst, size = 0.1, alpha = 0.3) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
  annotate("text", x = -Inf, y =  ymax * 0.92, label = "Summer run", hjust = -0.1, size = 4.5) +
  annotate("text", x = -Inf, y = -ymax * 0.92, label = "Fall run",   hjust = -0.1, size = 4.5) +
  scale_color_manual(values = c("0" = "grey50", "1" = "grey70"), guide = "none") +
  scale_fill_manual(name   = NULL,
                    values = setNames(c(col_run, col_run_fst, col_shared, col_shared_fst),
                                      c(lab_run, lab_run_fst, lab_shared, lab_shared_fst)),
                    breaks = c(lab_run, lab_run_fst, lab_shared, lab_shared_fst)) +
  guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1, stroke = 0))) +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr_index + 1,
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(-20, 20), labels = function(b) abs(b),
                     expand = c(0, 0)) +
  labs(x = "Chromosome", y = expression(-log[10](P))) +
  theme_bw(base_size = 16) +
  theme(legend.position        = "inside",
        legend.position.inside = c(0.995, 0.98),
        legend.justification   = c(1, 1),
        legend.background      = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.key             = element_blank(),
        legend.text            = element_text(size = 10),
        legend.key.size        = unit(0.4, "cm"),
        panel.grid             = element_blank(),
        axis.text.x            = element_text(size = 9))

ggsave("figures/04_gwas_summer_fall_miami_fst.png", miami,
       width = 30, height = 12, units = "cm", dpi = 900)

# Null distributions of the four overlap tests
null_df <- map_dfr(tests, ~ tibble(comparison = .x$comparison,
                                   shared_kb = attr(.x, "null") / 1e3))
p_null <- ggplot(null_df, aes(shared_kb)) +
  geom_histogram(bins = 40, fill = "grey60", color = "white") +
  geom_vline(data = test_tab, aes(xintercept = shared_kb), color = col_run, linewidth = 1) +
  facet_wrap(~ comparison, scales = "free") +
  labs(x = "Shared length (kb)", y = "Random shifts") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

ggsave("figures/06_gwas_summer_fall_fst_overlap_null.png", p_null,
       width = 18, height = 14, units = "cm", dpi = 600)