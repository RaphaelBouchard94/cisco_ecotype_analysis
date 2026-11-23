# ============================================================
# GWAS of length at maturity within one ecotype
#   - Manhattan plot, SNPs in permutation-significant windows in col_sig
#   - significant windows, with their peak SNP
#   - genes overlapping each significant window, with zebrafish symbols
#
# Windows are used as-is (no merging). Windows longer than a length cap
# are not used. The cap is applied here, ONCE, and the final flag is
# saved for every window in res/gwas_<run>_length_all_windows.tsv, so
# later scripts (overlap/candidate regions, GO, report) read the same
# window set instead of re-deriving it from the raw window file.
#
# The summer and fall scripts are identical except for the `run` line.
# ============================================================

library(tidyverse)
library(data.table)
library(GenomicRanges)
library(stringdist)

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------
run <- "fall"     # "summer" or "fall": the only line that differs between the two scripts

# sprintf() replaces each %s with the next argument, in order:
# sprintf("data/%s_windows.tsv", "summer") -> "data/summer_windows.tsv"
gwas_file    <- sprintf("data/gwas_length/%sgwas_ALL_CHR.lrt0.gz", run)
windows_file <- sprintf("data/gwas_length/%s_windows_ALL_CHR_10snp.tsv", run)  # output of perm_windows_combine.R
annot_file   <- "data/genome_annotation_table.tsv"
zf_file      <- "data/Danio_rerio.gene_info.gz"                                # NCBI Gene, zebrafish

out_allwin   <- sprintf("res/gwas_%s_length_all_windows.tsv", run)    # every window + final outlier99 flag
out_sigwin   <- sprintf("res/gwas_%s_length_sig_windows.tsv", run)    # significant windows only
out_genes    <- sprintf("res/gwas_%s_length_sig_window_genes.tsv", run)
out_plot     <- sprintf("figures/04_gwas_%s_length_perm_windows.png", run)

col_sig      <- "darkred"
p_snp        <- 0.05    # SNPs drawn in col_sig: inside a significant window AND P below this

# Window length cap: a 10-SNP window that spans a very long distance
# straddles a gap or a SNP-poor region, so its SNPs are not one local
# LD unit. Such windows are not used.
max_len_q    <- 0.95    # cap = this quantile of window length (computed per run)
max_len_bp   <- NA      # or a fixed cap in bp (e.g. from LD decay); replaces max_len_q when set

# Window design, used only for the LD diagnostic below
w_snps       <- 10      # SNPs per window
step_snps    <- 2       # SNPs between the starts of consecutive windows
ld_dist      <- 30000   # bp - placeholder, replace with your ngsLD decay distance

dir.create("res", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Large numbers with thousands separators: 1234567 -> "1,234,567"
fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

# ============================================================
# Load & prepare GWAS data
# ============================================================

gwas <- fread(gwas_file) %>%
  mutate(LRT = suppressWarnings(as.numeric(LRT))) %>%      # "-nan" -> NA
  filter(!is.na(LRT), LRT >= 0) %>%                         # also drops -999
  mutate(
    P          = pchisq(LRT, df = 1, lower.tail = FALSE),
    neglog10p  = -log10(P),
    Chromosome = as.character(Chromosome)
  )

chr_levels <- unique(gwas$Chromosome)

# Equal-width chromosome coordinates: every chromosome gets the same
# width on the x-axis (chrom_spacing), separated by chrom_gap, whatever
# its length in bp. A SNP's x position = start of its chromosome's slot
# + its relative position along the chromosome (0 to 1) * slot width.
chrom_spacing <- 1e7
chrom_gap     <- 2e6

chr_map <- tibble(Chromosome = chr_levels,
                  chr_index  = seq_along(chr_levels) - 1)

gwas <- gwas %>%
  group_by(Chromosome) %>%
  mutate(pos_scaled = Position / max(Position)) %>%
  ungroup() %>%
  left_join(chr_map, by = "Chromosome") %>%
  mutate(pos_equal = chr_index * (chrom_spacing + chrom_gap) +
           pos_scaled * chrom_spacing)

axis_df <- gwas %>%
  group_by(Chromosome, chr_index) %>%
  summarize(center = mean(pos_equal), .groups = "drop")

# ============================================================
# Permutation-significant windows and the length cap
# ============================================================

wins <- fread(windows_file)
stopifnot(all(c("Chromosome", "start_pos", "end_pos", "outlier99") %in% names(wins)))

wins <- wins %>%
  mutate(Chromosome = as.character(Chromosome),
         outlier99  = as.logical(outlier99),
         win_len    = end_pos - start_pos + 1)

if (!any(unique(wins$Chromosome) %in% chr_levels))
  stop("Chromosome names in windows file do not match the GWAS file")

# Length distribution of all windows vs outlier windows
wins %>%
  group_by(outlier99) %>%
  summarize(n = n(), mean_len = mean(win_len), median_len = median(win_len)) %>%
  print()
print(quantile(wins$win_len, c(0.90, 0.95, 0.99, 0.999)))

# The cap is computed over every row of the window file (all
# chromosomes, tested or not), as in the original analysis.
# names = FALSE makes quantile() return a plain number rather than one
# named "95%".
max_win_len <- if (is.na(max_len_bp)) {
  quantile(wins$win_len, max_len_q, names = FALSE)
} else {
  max_len_bp
}

# outlier99_raw = flag from the permutation step, kept for reference
# outlier99     = the flag every later step (and later script) uses
# R's & has three-valued logic: NA & FALSE is FALSE, NA & TRUE is NA.
# So an untested window (NA) stays NA, unless it is long, in which case
# it becomes FALSE. Later scripts drop long windows anyway.
wins <- wins %>%
  mutate(long_window   = win_len > max_win_len,
         outlier99_raw = outlier99,
         outlier99     = outlier99 & !long_window)

# Counts for the Results, on a consistent denominator: the percentage of
# significant windows is computed among usable windows (tested, not long)
tested <- !is.na(wins$outlier99_raw)
usable <- tested & !wins$long_window
n_raw  <- sum(wins$outlier99_raw, na.rm = TRUE)
n_used <- sum(wins$outlier99,     na.rm = TRUE)

cat(sprintf("\n%s run: %s windows tested; length cap %.1f kb removes %s (%s usable)\n",
            run, fmt_n(sum(tested)), max_win_len / 1e3,
            fmt_n(sum(tested & wins$long_window)), fmt_n(sum(usable))))
cat(sprintf("Outlier windows: %s from the permutations, %s removed as long, %s used = %.2f%% of usable windows (1%% expected by chance)\n",
            fmt_n(n_raw), fmt_n(n_raw - n_used), fmt_n(n_used), 100 * n_used / sum(usable)))

# LD diagnostic: how many contiguous windows does one LD block span?
# A block of ld_dist bp holds about ld_dist / bp_per_snp SNPs, and a
# window of w_snps SNPs sliding by step_snps fits
# (n_snps - w_snps) / step_snps + 1 times in n_snps SNPs.
# Use it to choose min_adjacent in the combined script (per run: SNP
# density, hence the answer, differs between summer and fall).
bp_per_snp <- median(wins$win_len) / (w_snps - 1)     # 10 SNPs = 9 gaps
k_ld <- max(1, ceiling((ld_dist / bp_per_snp - w_snps) / step_snps + 1))
cat(sprintf("Median spacing %.0f bp per SNP: an LD block of %.0f kb spans ~%d contiguous windows\n\n",
            bp_per_snp, ld_dist / 1e3, k_ld))

# Every window, with the final flag: read this file downstream, never
# the raw window file
fwrite(wins, out_allwin, sep = "\t")

sig_wins <- wins %>%
  filter(outlier99) %>%                 # filter() drops NA as well as FALSE
  mutate(window_id = row_number())

sig_wins_gr <- GRanges(
  seqnames  = sig_wins$Chromosome,
  ranges    = IRanges(sig_wins$start_pos, sig_wins$end_pos),
  window_id = sig_wins$window_id
)

# ============================================================
# Flag GWAS SNPs falling in any significant window
# ============================================================

snps_gr <- GRanges(gwas$Chromosome, IRanges(gwas$Position, width = 1))
hits    <- findOverlaps(snps_gr, sig_wins_gr)   # a SNP can sit in several (overlapping) windows

gwas$in_sig    <- overlapsAny(snps_gr, sig_wins_gr)
gwas$highlight <- gwas$in_sig & gwas$P < p_snp

# Peak SNP of each significant window.
# hits has one row per (SNP, window) pair: queryHits() gives the SNP
# (row of gwas), subjectHits() the window (row of sig_wins_gr).
window_stats <- tibble(
  snp_idx   = queryHits(hits),
  window_id = sig_wins_gr$window_id[subjectHits(hits)]
) %>%
  mutate(Position  = gwas$Position[snp_idx],
         neglog10p = gwas$neglog10p[snp_idx],
         pos_equal = gwas$pos_equal[snp_idx]) %>%
  group_by(window_id) %>%
  summarize(
    n_snps_tested  = n(),
    max_neglog10p  = max(neglog10p),
    peak_position  = Position[which.max(neglog10p)],
    peak_pos_equal = pos_equal[which.max(neglog10p)],
    .groups = "drop"
  )

sig_wins <- sig_wins %>% left_join(window_stats, by = "window_id")

# Length of the significant windows
sig_wins %>%
  summarize(n         = n(),
            median_kb = median(win_len) / 1e3,
            mean_kb   = mean(win_len) / 1e3,
            sd_kb     = sd(win_len) / 1e3) %>%
  print()

cat(sprintf("Significant windows: %s | SNPs in significant windows: %s of %s | highlighted (P < %g): %s\n",
            fmt_n(nrow(sig_wins)), fmt_n(sum(gwas$in_sig)), fmt_n(nrow(gwas)),
            p_snp, fmt_n(sum(gwas$highlight))))

# ============================================================
# Gene annotation & overlap with significant windows
# ============================================================

gene_annotation <- fread(annot_file) %>%
  dplyr::select(ScaffoldName, FromPosition, ToPosition, GeneName)

genes_gr <- GRanges(
  seqnames  = gene_annotation$ScaffoldName,
  ranges    = IRanges(gene_annotation$FromPosition, gene_annotation$ToPosition),
  gene_name = gene_annotation$GeneName
)

# Genes overlapping a significant window, even partly (type = "any").
# For genes lying entirely inside a window, use type = "within"
# (query = genes, subject = windows: the gene lies within the window).
ov <- findOverlaps(genes_gr, sig_wins_gr, type = "any")

overlap_df <- tibble(
  window_id = sig_wins_gr$window_id[subjectHits(ov)],
  gene_idx  = queryHits(ov)
) %>%
  left_join(sig_wins, by = "window_id") %>%
  mutate(
    gene_name    = genes_gr$gene_name[gene_idx],
    gene_start   = start(genes_gr)[gene_idx],
    gene_end     = end(genes_gr)[gene_idx],
    peak_in_gene = peak_position >= gene_start & peak_position <= gene_end
  ) %>%
  filter(gene_name != "-", gene_name != "") %>%
  dplyr::select(-gene_idx) %>%
  arrange(desc(max_neglog10p))

# ============================================================
# Fuzzy gene name mapping (NCBI Gene, zebrafish)
# ============================================================

clean_string <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all("\\(.*?\\)", "") %>%    # drop text in brackets
    str_replace_all("[^a-z0-9 ]", " ") %>%  # punctuation -> space
    str_squish()                            # collapse repeated spaces
}

if (!file.exists(zf_file)) {
  options(timeout = 600)
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Non-mammalian_vertebrates/Danio_rerio.gene_info.gz",
    zf_file, mode = "wb")
}

zf_genes <- fread(zf_file, sep = "\t", quote = "") %>%
  filter(type_of_gene == "protein-coding") %>%
  transmute(
    ensembl_gene_id    = str_extract(dbXrefs, "ENSDARG\\d+"),
    external_gene_name = Symbol,
    description        = description
  ) %>%
  mutate(target_desc = clean_string(description)) %>%
  filter(target_desc != "",
         !str_detect(target_desc, "prothymosin|karyopherin"))

# amatch() returns, for each query, the row of zf_genes$target_desc with
# the smallest Jaro-Winkler distance (the first such row when tied).
# Same answer as computing every distance, sorting and keeping the first
# row (the previous version), but in one vectorised call instead of
# building and sorting the whole zebrafish table once per gene.
# maxDist = Inf: always return a match. match_dist is saved so weak
# matches (large values) can be checked by hand.
gene_map <- overlap_df %>%
  distinct(gene_name) %>%
  mutate(query = clean_string(gene_name)) %>%
  filter(query != "") %>%
  mutate(hit                = amatch(query, zf_genes$target_desc, method = "jw", maxDist = Inf),
         external_gene_name = zf_genes$external_gene_name[hit],
         ensembl_gene_id    = zf_genes$ensembl_gene_id[hit],
         zf_description     = zf_genes$description[hit],
         match_dist         = stringdist(query, zf_genes$target_desc[hit], method = "jw")) %>%
  dplyr::select(-query, -hit)

# ============================================================
# Plot
# ============================================================

gwas <- gwas %>% mutate(chr_shade = factor(chr_index %% 2))   # alternate grey shades

gwas_plot <- ggplot() +
  geom_point(data = filter(gwas, !highlight),
             aes(pos_equal, neglog10p, color = chr_shade),
             size = 0.2, alpha = 0.6) +
  geom_point(data = filter(gwas, highlight),          # drawn last = on top
             aes(pos_equal, neglog10p), size = 0.2, alpha = 0.6,
             color = col_sig) +
  scale_color_manual(values = c("0" = "grey50", "1" = "grey70")) +
  scale_x_continuous(breaks = axis_df$center,
                     labels = axis_df$chr_index + 1,
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, max(gwas$neglog10p) * 1.05)) +
  labs(x = "Chromosome", y = expression(-log[10](P))) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        panel.grid = element_blank())

ggsave(out_plot, gwas_plot, width = 30, height = 7, units = "cm", dpi = 900)

# ============================================================
# Tables
# ============================================================

fwrite(sig_wins, out_sigwin, sep = "\t")
fwrite(overlap_df %>% left_join(gene_map, by = "gene_name"), out_genes, sep = "\t")

cat(sprintf("Written: %s, %s, %s, %s\n", out_allwin, out_sigwin, out_genes, out_plot))
