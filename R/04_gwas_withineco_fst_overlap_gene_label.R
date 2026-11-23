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
# 3. Candidate genes, following Montejo-Kovacevich et al. 2021
#    (Mol Ecol 30: 6387-6402), separately for each run:
#      - candidate region = >= min_adjacent contiguous outlier windows,
#        at least one of them with median p < max_median_p
#      - core interval    = between the strongest and second strongest
#        SNP of the region
#      - candidate genes  = the gene containing the core interval, or
#        the two genes closest to it when it is intergenic or spans
#        more than one gene
# 4. Miami plot: summer up, fall down, exported at its printed size in
#    Fig. 3e (candidate gene labels optional, off by default)
# 5. Numbers for the Results section, in one report
#
# Inputs, all written by the two 04_gwas_<run>_perm.R scripts:
#   res/gwas_<run>_length_sig_windows.tsv       significant windows
#   res/gwas_<run>_length_all_windows.tsv       every window, with the
#       final outlier99 flag (windows longer than the length cap already
#       set to FALSE), outlier99_raw, long_window and win_len
#   res/gwas_<run>_length_sig_window_genes.tsv  genes per significant window
# The outlier flag is never re-derived here, so every section uses the
# same set of significant windows. Windows are compared by genomic
# coordinates, since the windows of the two runs are built from
# different SNP sets.
# ============================================================

library(tidyverse)
library(data.table)
library(GenomicRanges)
library(stringdist)
library(ggrepel)

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------
summer_gwas_file   <- "data/gwas_length/summergwas_ALL_CHR.lrt0.gz"
fall_gwas_file     <- "data/gwas_length/fallgwas_ALL_CHR.lrt0.gz"
summer_win_file    <- "res/gwas_summer_length_sig_windows.tsv"
fall_win_file      <- "res/gwas_fall_length_sig_windows.tsv"
summer_allwin_file <- "res/gwas_summer_length_all_windows.tsv"
fall_allwin_file   <- "res/gwas_fall_length_all_windows.tsv"
summer_gene_file   <- "res/gwas_summer_length_sig_window_genes.tsv"
fall_gene_file     <- "res/gwas_fall_length_sig_window_genes.tsv"
fst_out_file       <- "data/genome_scan/notsum_rup_merged_top_candidate.tsv"
fst_scan_file      <- "data/genome_scan/notsummer_rup_genome_wide.25kb_windows.fst"
annot_file         <- "data/genome_annotation_table.tsv"
zf_file            <- "data/Danio_rerio.gene_info.gz"   # NCBI Gene, zebrafish (downloaded by the per-run scripts)

# Overlap tests
n_shift <- 1000      # circular shifts per test
seed    <- 2026      # random seed for the shifts

# Candidate regions (Montejo-Kovacevich et al. 2021). 10 is their value
# for 50-SNP windows at ~33 bp/SNP; derive yours from LD decay (the
# per-run scripts print how many windows one LD block spans).
min_adjacent <- 5      # min. number of contiguous outlier windows in a region
max_median_p <- 0.01   # at least one window of the region must have median p below this

# Miami plot
label_genes      <- FALSE   # TRUE: write the candidate genes on the plot (needs a taller plot than Fig. 3e)
max_genes_label  <- 3       # genes written per label on the plot; the rest become "+n"
label_classes    <- c("shared_fst", "shared", "run_fst", "run")   # region classes labelled on the plot
drop_label_regex <- "^LOC"  # gene names matching this are not shown on the plot (uncharacterised genes)
y_lim            <- NA      # y-axis runs from -y_lim to y_lim; NA = just above the strongest SNP of either run
label_nudge      <- 3       # how far labels are pushed away from their peak (in -log10 P units)
shade_fst        <- FALSE   # TRUE: shade Fst outlier regions behind the points

# Size and text of the exported Miami plot.
# The PNG is saved at the size it is printed: the slot of panel e in
# Fig. 3 (Figure_031.png is 183 x 170 mm; panel e is its bottom 49 mm).
# Placed in the figure at 100%, nothing is shrunk, so the font sizes
# below are the sizes in the printed paper.
miami_width_cm  <- 18.3
miami_height_cm <- 4.9
miami_dpi       <- 867      # resolution of Figure_031.png, so the pixels line up 1:1
font_pt         <- 8        # axis numbers and legend text, in points (other Fig. 3 panels: ~8 pt)
title_pt        <- 9        # axis titles and run names, in points
font_family     <- ""       # "" = R's default sans-serif; or a font installed on your computer, e.g. "Open Sans"

# Results report (section 5)
# Phenotype files given to ANGSD (-yQuant): one value per fish, -999 =
# missing. Used only to count the fish in each GWAS.
pheno_files <- c(summer = "data/gwas_length/summer_pheno.yquant",
                 fall   = "data/gwas_length/fall_pheno.yquant")
lambda_file <- "res/qq_window_median_p_lambda.tsv"   # written by qq_windows.R
go_file     <- "res/gwas_length_go_enrichment.tsv"   # written by 05_gwas_length_go.R
go_fdr      <- 0.1      # GO terms listed if their FDR is below this
inv_label   <- 36       # chromosome number of the inversion, as labelled on the plots (Fig. 4g)
inv_range   <- c(34.8e6, 35.6e6)
n_top       <- 10       # strongest candidate regions listed per run
weak_peak   <- 2        # a significant window is "weak" if its best SNP has -log10 P below this
report_file <- "res/gwas_results_numbers.txt"

dir.create("res", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ============================================================
# Load data
# ============================================================

read_gwas <- function(f) {
  fread(f) %>%
    mutate(LRT = suppressWarnings(as.numeric(LRT))) %>%   # "-nan" -> NA
    filter(!is.na(LRT), LRT >= 0) %>%                      # also drops -999
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

fst_scan <- fread(fst_scan_file)[, -1]           # drop the first (index) column
colnames(fst_scan) <- c("Chromosome", "Position", "nsites", "fst")
fst_scan$Chromosome <- as.character(fst_scan$Chromosome)

# Chromosomes present in both GWAS and the Fst scan
chroms <- Reduce(intersect, list(unique(gwas_s$Chromosome),
                                 unique(gwas_f$Chromosome),
                                 unique(fst_scan$Chromosome)))
if (length(chroms) == 0) stop("No chromosome names shared by the GWAS and Fst files")

# Every window of one run, with the flags written by the per-run script.
# fread(select = ...) only WARNS when a column is missing and carries on,
# so the explicit check makes an old window file (e.g. the raw
# *_windows_ALL_CHR_10snp.tsv, which has no long_window column) stop the
# script instead of silently bringing back the unfiltered flag.
read_allwin <- function(f) {
  need <- c("Chromosome", "start_pos", "end_pos", "win_len",
            "outlier99", "outlier99_raw", "long_window")
  w <- suppressWarnings(fread(f, select = need))
  if (!all(need %in% names(w)))
    stop(f, " lacks ", paste(setdiff(need, names(w)), collapse = ", "),
         ": use res/gwas_<run>_length_all_windows.tsv written by 04_gwas_<run>_perm.R")
  w %>%
    mutate(Chromosome = as.character(Chromosome),
           across(c(outlier99, outlier99_raw, long_window), as.logical)) %>%
    filter(Chromosome %in% chroms)
}

allwin_s <- read_allwin(summer_allwin_file)
allwin_f <- read_allwin(fall_allwin_file)

# Chromosome length = furthest position seen in any input.
# The braces { } after the last pipe make `.` stand for the whole table,
# so setNames(.$len, .$Chromosome) turns it into a named vector:
# c(chrA = 81234567, chrB = ...), looked up later as chr_len_v[["chrA"]].
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

# The significant-window file and the all-window file come from the same
# run of the per-run script, so their outlier counts must agree. If they
# do not, one of the two files is stale.
check_sync <- function(allw, win_gr, run) {
  n_all <- sum(allw$outlier99, na.rm = TRUE)
  if (n_all != length(win_gr))
    stop(sprintf("%s: %d outlier windows in the all-window file but %d in the significant-window file; rerun 04_gwas_%s_perm.R",
                 run, n_all, length(win_gr), run))
}
check_sync(allwin_s, win_s_gr, "summer")
check_sync(allwin_f, win_f_gr, "fall")

reg_s      <- reduce(win_s_gr)                 # bp covered by significant summer windows
reg_f      <- reduce(win_f_gr)                 # bp covered by significant fall windows
reg_fst    <- reduce(to_gr(fst_out$Chromosome, fst_out$start, fst_out$end))
reg_sf     <- GenomicRanges::intersect(reg_s, reg_f, ignore.strand = TRUE)       # significant in both
reg_sf_fst <- GenomicRanges::intersect(reg_sf, reg_fst, ignore.strand = TRUE)    # ... and in an Fst outlier region

bp <- function(gr) sum(as.numeric(width(gr)))
genome_bp <- sum(as.numeric(chr_len_v))

cat(sprintf("Genome analysed: %.1f Mb on %d chromosomes\n", genome_bp / 1e6, length(chroms)))
cat(sprintf("Summer: %d significant windows covering %.2f Mb (%.2f%%)\n",
            length(win_s_gr), bp(reg_s) / 1e6, 100 * bp(reg_s) / genome_bp))
cat(sprintf("Fall:   %d significant windows covering %.2f Mb (%.2f%%)\n",
            length(win_f_gr), bp(reg_f) / 1e6, 100 * bp(reg_f) / genome_bp))
cat(sprintf("Fst outliers: %d regions covering %.2f Mb (%.2f%%)\n\n",
            length(reg_fst), bp(reg_fst) / 1e6, 100 * bp(reg_fst) / genome_bp))

# ============================================================
# 1. Overlap tests (circular shift null)
# ============================================================

# Each region set is shifted by a random offset along its chromosome,
# wrapping round the end: a region pushed past the end is split, its
# overhang placed at the start. %% is the remainder after division, so
# (position - 1 + shift) %% L + 1 always lands between 1 and L.
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

# The seed is set right before the only random step of the script, so
# the shifts depend on the seed alone, not on code that runs earlier.
set.seed(seed)
tests <- list(
  overlap_test("fall x summer",          reg_f,  reg_s,   win_f_gr),
  overlap_test("fall x Fst",             reg_f,  reg_fst, win_f_gr),
  overlap_test("summer x Fst",           reg_s,  reg_fst, win_s_gr),
  overlap_test("(fall & summer) x Fst",  reg_sf, reg_fst)
)

test_tab <- bind_rows(tests)
fwrite(test_tab, "res/gwas_summer_fall_fst_overlap_tests.tsv", sep = "\t")
print(test_tab %>% mutate(across(where(is.double), ~ signif(.x, 3))), width = Inf)

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

# Genes of the per-run gene tables overlapping each region, named by
# their zebrafish symbol when one was matched
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
  
  shared_regions <- tibble(
    region           = seq_along(reg_sf),
    Chromosome       = as.character(seqnames(reg_sf)),
    start            = start(reg_sf),
    end              = end(reg_sf),
    width_kb         = width(reg_sf) / 1e3,
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
# 3. Candidate genes (Montejo-Kovacevich et al. 2021)
# ============================================================
#
# Done separately for each run, in three steps:
#   a. Candidate regions: stretches of >= min_adjacent contiguous outlier
#      windows (outlier99) in which at least one window has a median
#      p-value < max_median_p
#   b. Core interval: from the SNP with the strongest association to the
#      SNP with the second strongest association in the region
#   c. Candidate genes for the core interval:
#        within_gene   it lies inside a gene              -> that gene
#        otherwise                                        -> the two genes closest to it
#                                                            (a gene overlapping it is at distance 0)
#          intergenic    it overlaps no gene
#          partly_genic  it overlaps one gene but is not inside it
#          across_genes  it overlaps two or more genes (the two closest
#                        to the strongest SNP are kept)
#
# Each region also gets a class, compared with the significant WINDOWS
# of sections 1-2, so it matches the Miami plot colours:
#   shared_fst  overlaps significant windows of the other run, in an Fst outlier region
#   shared      overlaps significant windows of the other run
#   run_fst     significant in this run only, overlaps an Fst outlier region
#   run         significant in this run only

# Genes without a name ("-") are dropped first, so a candidate gene is
# never an unnamed one
gene_annotation <- fread(annot_file) %>%
  transmute(Chromosome = as.character(ScaffoldName),
            gene_start = FromPosition,
            gene_end   = ToPosition,
            gene_name  = GeneName) %>%
  filter(gene_name != "-", gene_name != "", Chromosome %in% chroms)

# --- a + b. Candidate regions and core intervals of one run --
# allw: every window of this run (read_allwin); gwas: this run's SNPs;
# other: significant windows of the other run (reduced)
find_regions <- function(allw, gwas, other, run) {
  # rleid() starts a new id every time Chromosome or outlier99 changes
  # from one row to the next:  outlier99 = F F T T T F T  ->  1 1 2 2 2 3 4
  # Once the windows are sorted along each chromosome, one id is one
  # stretch of contiguous windows with the same status, and n() counts
  # its windows. Any window that is not a usable outlier breaks a
  # stretch instead of joining two:
  #   - an untested window (NA), set to FALSE here
  #   - a window longer than the length cap, already FALSE in the file.
  #     It spans a gap in the SNP data, so its neighbours are not
  #     treated as contiguous across it.
  wins <- allw %>%
    mutate(outlier99 = coalesce(outlier99, FALSE)) %>%
    arrange(Chromosome, start_pos) %>%
    mutate(stretch = rleid(Chromosome, outlier99)) %>%
    group_by(stretch) %>%
    mutate(n_adjacent = n()) %>%
    ungroup() %>%
    filter(outlier99, n_adjacent >= min_adjacent)
  
  if (nrow(wins) == 0) return(NULL)
  snp_gr <- GRanges(gwas$Chromosome, IRanges(gwas$Position, width = 1))
  
  # Median p-value of each window, from the GWAS SNPs it contains
  win_gr <- GRanges(wins$Chromosome, IRanges(wins$start_pos, wins$end_pos))
  h <- findOverlaps(win_gr, snp_gr)
  win_p <- tibble(w = queryHits(h), P = gwas$P[subjectHits(h)]) %>%
    group_by(w) %>%
    summarize(median_p = median(P), .groups = "drop")
  
  # One row per stretch; keep stretches with a window below max_median_p
  regions <- wins %>%
    mutate(w = row_number()) %>%
    left_join(win_p, by = "w") %>%
    group_by(Chromosome, stretch) %>%
    summarize(start     = min(start_pos),
              end       = max(end_pos),
              n_windows = n(),
              min_window_median_p = min(median_p, na.rm = TRUE),
              .groups = "drop") %>%
    filter(min_window_median_p < max_median_p)
  
  cat(sprintf("%s: %d stretches of >= %d contiguous outlier windows, %d with a window median p < %g\n",
              run, n_distinct(wins$stretch), min_adjacent, nrow(regions), max_median_p))
  if (nrow(regions) == 0) return(NULL)
  
  # Strongest and second strongest SNP of each region.
  # SNPs are ranked by LRT rather than P: for very strong associations
  # pchisq() returns P = 0, and SNPs tied at 0 could not be ranked.
  # arrange() sorts the SNPs within each region, and group_by() keeps
  # that order, so in summarize() [1] is the top SNP and [2] the second.
  reg_gr <- GRanges(regions$Chromosome, IRanges(regions$start, regions$end))
  h <- findOverlaps(reg_gr, snp_gr)
  top2 <- tibble(r         = queryHits(h),
                 Position  = gwas$Position[subjectHits(h)],
                 LRT       = gwas$LRT[subjectHits(h)],
                 neglog10p = gwas$neglog10p[subjectHits(h)]) %>%
    arrange(r, desc(LRT)) %>%
    group_by(r) %>%
    summarize(n_snps           = n(),
              top_pos          = Position[1],
              max_neglog10p    = neglog10p[1],
              second_pos       = Position[2],
              second_neglog10p = neglog10p[2],
              .groups = "drop")
  
  # left_join() keeps the row order of `regions`, which is also the
  # order of reg_gr, so the overlapsAny() results line up with the rows
  regions %>%
    mutate(r = row_number()) %>%
    left_join(top2, by = "r") %>%
    mutate(run         = run,
           region_id   = sprintf("%s_%03d", run, r),
           width_kb    = (end - start + 1) / 1e3,
           core_start  = pmin(top_pos, second_pos, na.rm = TRUE),
           core_end    = pmax(top_pos, second_pos, na.rm = TRUE),
           core_kb     = (core_end - core_start + 1) / 1e3,
           shared      = overlapsAny(reg_gr, other),
           fst_outlier = overlapsAny(reg_gr, reg_fst),
           class       = case_when(overlapsAny(reg_gr, reg_sf_fst) ~ "shared_fst",
                                   shared                          ~ "shared",
                                   fst_outlier                     ~ "run_fst",
                                   TRUE                            ~ "run")) %>%
    dplyr::select(run, region_id, Chromosome, start, end, width_kb, n_windows, n_snps,
                  min_window_median_p, top_pos, max_neglog10p, second_pos, second_neglog10p,
                  core_start, core_end, core_kb, shared, fst_outlier, class)
}

# --- c. Candidate genes of one core interval -----------------
genes_for_region <- function(chr, core_start, core_end, top_pos) {
  # Distance between two intervals: if the gene is to the right of the
  # core, gene_start - core_end is positive; if it is to the left,
  # core_start - gene_end is positive; if they overlap, both are <= 0.
  # pmax(0, a, b) therefore gives the gap in bp, or 0 for an overlap.
  # Genes are ranked by distance to the core, then by distance to the top
  # SNP. Exact ties are then broken by fixed rules, so the result never
  # depends on the row order of the annotation file (the GO script uses
  # the same order and must pick the same genes):
  #   !inside_gene  a copy containing the core comes first (FALSE sorts
  #                 before TRUE, so !inside_gene puts inside copies first)
  #   gene_start, gene_end, gene_name   then position, then name
  # distinct() then keeps the first copy of each gene name, so two
  # annotations of the same gene cannot fill both places.
  g <- gene_annotation %>%
    filter(Chromosome == chr) %>%
    mutate(inside_gene  = gene_start <= core_start & gene_end >= core_end,
           distance_bp  = pmax(0, gene_start - core_end, core_start - gene_end),
           dist_top_snp = pmax(0, gene_start - top_pos, top_pos - gene_end)) %>%
    arrange(distance_bp, dist_top_snp, !inside_gene, gene_start, gene_end, gene_name) %>%
    distinct(gene_name, .keep_all = TRUE)
  
  if (nrow(g) == 0) return(g)            # no named gene on this chromosome
  
  n_overlap <- sum(g$distance_bp == 0)
  if (any(g$inside_gene)) {
    # Usually one gene; two only if annotated genes overlap each other
    g %>% filter(inside_gene) %>% mutate(case = "within_gene")
  } else {
    g %>%
      slice_head(n = 2) %>%
      mutate(case = case_when(n_overlap == 0 ~ "intergenic",
                              n_overlap == 1 ~ "partly_genic",
                              TRUE           ~ "across_genes"))
  }
}

cand_regions <- bind_rows(
  find_regions(allwin_s, gwas_s, reg_f, "summer"),
  find_regions(allwin_f, gwas_f, reg_s, "fall")
)
if (nrow(cand_regions) == 0)
  stop("No candidate region in either run; check min_adjacent / max_median_p")

# pmap_dfr() calls the function once per row of the table it receives,
# passing each column as the argument of the same name, then stacks the
# results into one table (like a for loop over rows + bind_rows()).
cand_genes <- cand_regions %>%
  dplyr::select(run, region_id, Chromosome, core_start, core_end, top_pos) %>%
  pmap_dfr(function(run, region_id, Chromosome, core_start, core_end, top_pos) {
    genes_for_region(Chromosome, core_start, core_end, top_pos) %>%
      mutate(run = run, region_id = region_id, .before = 1)
  }) %>%
  dplyr::select(run, region_id, case, gene_name, gene_start, gene_end, distance_bp)

# --- Zebrafish names (same fuzzy matching as the per-run scripts) ---
clean_string <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all("\\(.*?\\)", "") %>%    # drop text in brackets
    str_replace_all("[^a-z0-9 ]", " ") %>%  # punctuation -> space
    str_squish()                            # collapse repeated spaces
}

zf_genes <- fread(zf_file, sep = "\t", quote = "") %>%
  filter(type_of_gene == "protein-coding") %>%
  transmute(ensembl_gene_id = str_extract(dbXrefs, "ENSDARG\\d+"),
            zf_symbol       = Symbol,
            zf_description  = description,
            target_desc     = clean_string(description)) %>%
  filter(target_desc != "",
         !str_detect(target_desc, "prothymosin|karyopherin"))

# amatch() gives, for each query, the row of zf_genes with the smallest
# Jaro-Winkler distance (the first one if tied). With maxDist = Inf it
# always returns a match, however poor, so match_dist is kept to check
# weak matches (large values) by hand.
gene_map <- tibble(gene_name = unique(na.omit(cand_genes$gene_name))) %>%
  mutate(query = clean_string(gene_name)) %>%
  filter(query != "") %>%
  mutate(hit             = amatch(query, zf_genes$target_desc, method = "jw", maxDist = Inf),
         zf_symbol       = zf_genes$zf_symbol[hit],
         ensembl_gene_id = zf_genes$ensembl_gene_id[hit],
         zf_description  = zf_genes$zf_description[hit],
         match_dist      = stringdist(query, zf_genes$target_desc[hit], method = "jw")) %>%
  dplyr::select(-query, -hit)

cand_genes <- cand_genes %>%
  left_join(gene_map, by = "gene_name") %>%
  # name shown on the plot: zebrafish symbol, or the (shortened) original
  # annotation when no match was found
  mutate(label = coalesce(zf_symbol, str_trunc(gene_name, 25))) %>%
  arrange(run, region_id, gene_start)                 # genes in genomic order

# --- One row per region, with its genes and plot label ------
# Plot label: names matching drop_label_regex are left out, at most
# max_genes_label names, and intergenic regions in brackets: "(left, right)"
make_label <- function(lab, case) {
  lab <- unique(lab[!is.na(lab) & !str_detect(lab, drop_label_regex)])
  if (length(lab) == 0) return(NA_character_)
  txt <- paste(head(lab, max_genes_label), collapse = ", ")
  if (length(lab) > max_genes_label) txt <- paste0(txt, ", +", length(lab) - max_genes_label)
  if (case == "intergenic") txt <- paste0("(", txt, ")")
  txt
}

region_genes <- cand_genes %>%
  group_by(run, region_id) %>%
  summarize(case       = case[1],
            n_genes    = n(),
            genes      = paste(label, collapse = ", "),
            plot_label = make_label(label, case[1]),
            .groups = "drop")

# Regions without any named gene on their chromosome keep NA here
cand_regions <- cand_regions %>%
  left_join(region_genes, by = c("run", "region_id")) %>%
  mutate(class = factor(class, levels = c("shared_fst", "shared", "run_fst", "run"))) %>%
  arrange(run, class, desc(max_neglog10p))

cat("\nCandidate regions per run and class:\n")
print(count(cand_regions, run, class))
cat("\nCore interval relative to genes:\n")
print(count(cand_regions, run, case))

fwrite(cand_regions, "res/gwas_length_candidate_regions.tsv", sep = "\t")
fwrite(cand_genes,   "res/gwas_length_candidate_region_genes.tsv", sep = "\t")

# ============================================================
# 4. Miami plot
# ============================================================

chrom_spacing <- 1e7
chrom_gap     <- 2e6

chr_order <- unique(c(gwas_s$Chromosome, gwas_f$Chromosome))
chr_order <- chr_order[chr_order %in% chroms]
chr_map <- tibble(Chromosome = chr_order,
                  chr_index  = seq_along(chr_order) - 1,
                  max_pos    = chr_len_v[chr_order])

# Genomic position -> x position with equal-width chromosomes
to_equal <- function(chr, pos) {
  m <- chr_map[match(chr, chr_map$Chromosome), ]
  m$chr_index * (chrom_spacing + chrom_gap) + pmin(pos / m$max_pos, 1) * chrom_spacing
}

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

# y-axis limit. With y_lim = NA it is set just above the strongest SNP of
# either run, so the data fill both halves instead of leaving the empty
# headroom the gene labels needed. is.finite() leaves out SNPs with
# P = 0 (-log10 P = Inf), which would otherwise make the maximum Inf.
# Everything below (labels, run names) reads y_lim, so it is set once here.
if (is.na(y_lim)) {
  y_lim <- ceiling(1.05 * max(plot_df$neglog10p[is.finite(plot_df$neglog10p)]))
}
cat(sprintf("Miami plot y-axis: -%g to %g\n", y_lim, y_lim))

axis_df <- chr_map %>%
  mutate(center = chr_index * (chrom_spacing + chrom_gap) + chrom_spacing / 2)

# Optional shading of the Fst outlier regions (shade_fst = TRUE).
# A layer set to NULL adds nothing to the plot: `p + NULL` is just `p`.
fst_rect <- tibble(Chromosome = as.character(seqnames(reg_fst)),
                   start = start(reg_fst), end = end(reg_fst)) %>%
  mutate(xmin = to_equal(Chromosome, start),
         xmax = pmax(to_equal(Chromosome, end), xmin + 3e4))   # keep narrow regions visible
fst_layer <- if (shade_fst) {
  geom_rect(data = fst_rect, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey88", colour = NA)
} else NULL

# Gene labels: one per candidate region, anchored at the region's
# strongest SNP. Summer labels point up, fall labels point down.
#   y     anchors are capped at y_lim so a peak above the axis limit
#         still gets its label
#   nudge the label's starting point is moved by `nudge` BEFORE the axis
#         limits are applied, and ggplot deletes any label that lands
#         outside them. So the nudge shrinks for peaks close to the limit
#         (0 for a peak at the limit) and ggrepel then moves the label
#         back inside the plot itself.
lab_df <- cand_regions %>%
  filter(class %in% label_classes, !is.na(plot_label)) %>%
  mutate(sign  = if_else(run == "summer", 1, -1),
         x     = to_equal(Chromosome, top_pos),
         y     = sign * pmin(max_neglog10p, y_lim),
         nudge = sign * pmax(0, pmin(label_nudge, 0.95 * y_lim - abs(y))))

# A function that returns a ggplot layer, so the summer and fall labels
# share all their settings and differ only in direction (+1 up, -1 down).
#   nudge_y       one value per label, taken from the data (see above)
#   ylim          keeps each run's labels on its own half of the plot
#   max.overlaps  Inf: never silently drop a label (default drops crowded ones)
#   seed          same label placement every time the script runs
gene_label_layer <- function(df, direction) {
  if (nrow(df) == 0) return(NULL)   # no labels on this side; `+ NULL` adds nothing
  geom_text_repel(data = df, aes(x, y, label = plot_label),
                  nudge_y = df$nudge,
                  ylim = if (direction > 0) c(0, y_lim) else c(-y_lim, 0),
                  size = 2.6, fontface = "italic", colour = "grey10",
                  segment.colour = "grey40", min.segment.length = 0,
                  box.padding = 0.3, max.overlaps = Inf, seed = 1)
}

col_run        <- "#ABC7D4FF"   # significant in one run
col_run_fst    <- "#FF9E9DFF"   # significant in one run and inside an Fst outlier region
col_shared     <- "#175149FF"   # significant in both runs
col_shared_fst <- "darkred"     # significant in both runs and inside an Fst outlier region

# Keys that link each layer to its colour in scale_fill_manual()
lab_run        <- "Significant in one ecotype"
lab_shared     <- "Significant in both ecotypes"
lab_run_fst    <- "Significant in one ecotype within Fst outlier window"
lab_shared_fst <- "Significant in both ecotypes within Fst outlier window"

# Text shown in the legend, in the same order as `breaks` below. The keys
# above stay as they are; the legend only shows these shorter texts under
# the title "Significant in:", so it fits on one line under the plot.
# expression() is R's math notation (see ?plotmath): F[ST] writes ST as a
# subscript, and * joins the pieces without adding a space.
legend_text <- c(expression("one ecotype"),
                 expression("one ecotype, " * F[ST] * " outlier window"),
                 expression("both ecotypes"),
                 expression("both ecotypes, " * F[ST] * " outlier window"))

miami <- ggplot() +
  fst_layer +
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
  # `if` without `else` gives NULL when label_genes is FALSE, and
  # `+ NULL` adds nothing to the plot (same trick as fst_layer above)
  { if (label_genes) gene_label_layer(filter(lab_df, run == "summer"),  1) } +
  { if (label_genes) gene_label_layer(filter(lab_df, run == "fall"),   -1) } +
  scale_color_manual(values = c("0" = "grey50", "1" = "grey70"), guide = "none") +
  scale_fill_manual(name   = "Significant in:",
                    values = setNames(c(col_run, col_run_fst, col_shared, col_shared_fst),
                                      c(lab_run, lab_run_fst, lab_shared, lab_shared_fst)),
                    breaks = c(lab_run, lab_run_fst, lab_shared, lab_shared_fst),
                    labels = legend_text) +
  # nrow = 1: all four keys side by side in a single row
  guides(fill = guide_legend(nrow = 1,
                             override.aes = list(size = 2, alpha = 1, stroke = 0))) +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr_index + 1,
                     expand = c(0.01, 0.01)) +
  # The run names sit on a second y-axis on the right, outside the
  # data, so no peak can hide them. sec_axis(~ .) copies the left axis
  # (~ . means "same values"); its only two "ticks" are the middle of
  # each half (+y_lim/2 and -y_lim/2), labelled with the run names.
  scale_y_continuous(limits = c(-y_lim, y_lim), labels = function(b) abs(b),
                     expand = c(0, 0),
                     sec.axis = sec_axis(~ ., breaks = c(y_lim, -y_lim) / 2,
                                         labels = c("Summer run", "Fall run"))) +
  labs(x = "Chromosome", y = expression(-log[10](P))) +
  # Text sizes in theme() are in points. Because the plot is saved at its
  # printed size, 8 pt here is 8 pt in the paper. (Sizes inside geoms,
  # e.g. geom_text(size = ), are in mm instead: size = pt / .pt.)
  theme_bw(base_size = title_pt, base_family = font_family) +
  theme(legend.position    = "bottom",
        legend.title       = element_text(size = font_pt),
        legend.text        = element_text(size = font_pt),
        legend.key         = element_blank(),
        legend.key.size    = unit(3, "mm"),
        legend.margin      = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "mm"),          # gap between the x-axis title and the legend
        panel.grid         = element_blank(),
        axis.title         = element_text(size = title_pt),
        axis.text          = element_text(size = font_pt),
        # Run names: rotated to read downwards, centred on their half
        axis.text.y.right  = element_text(size = title_pt, angle = -90, hjust = 0.5,
                                          colour = "black", margin = margin(l = 1, unit = "mm")),
        axis.ticks.y.right = element_blank(),
        plot.margin        = margin(1, 1, 0.5, 1, unit = "mm"))

ggsave("figures/04_gwas_summer_fall_miami_fst.png", miami,
       width = miami_width_cm, height = miami_height_cm, units = "cm", dpi = miami_dpi)

# Null distributions of the four overlap tests, with the observed value
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

# ============================================================
# 5. Numbers for the Results section
# ============================================================
#
# Collects every number the Results text needs into one report,
# printed to the console and saved to:
#   res/gwas_results_numbers.txt              the report
#   res/gwas_results_run_stats.tsv            one row per run
#   res/gwas_results_inversion.tsv            chr36 inversion, one row per run
#   res/supp_table_gwas_candidate_regions.tsv SI table of candidate regions

# The inversion's chromosome: the one labelled inv_label on the plots.
# The report prints its name so you can check it is the one in Fig. 4.
inv_chr <- chr_map$Chromosome[chr_map$chr_index + 1 == inv_label]
if (length(inv_chr) != 1)
  stop("Could not find chromosome ", inv_label, " in chr_map; set inv_chr by hand")

# ------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------
# Chromosome name -> number used on the plots (1, 2, ... 38)
chr_num <- function(chr) chr_map$chr_index[match(chr, chr_map$Chromosome)] + 1

# Large numbers with thousands separators: 1234567 -> "1,234,567"
fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

# Empirical P values from n_rand randomisations cannot go below
# 1 / (n_rand + 1). A P value at that floor is written "P <= floor",
# because it is a bound, not an estimate.
fmt_p <- function(p, n_rand) {
  floor_p <- 1 / (n_rand + 1)
  ifelse(p <= floor_p * (1 + 1e-9),
         sprintf("P <= %s", signif(floor_p, 1)),
         sprintf("P = %s", signif(p, 2)))
}

# min() / max() / median() of an empty vector give -Inf/Inf/NA plus a
# warning. safe() returns NA instead, so a run with no candidate region
# does not stop the report.
safe <- function(f, x) if (length(x) > 0) f(x) else NA_real_

# Number of fish with a phenotype in an ANGSD -yQuant file
count_fish <- function(f) {
  if (is.na(f) || !file.exists(f)) return(NA_integer_)
  y <- scan(f, quiet = TRUE)        # reads every number in the file into a vector
  sum(y != -999)
}

# Strongest SNP (max -log10 P) inside each window of win_gr.
# The tibble has one row per (window, SNP) pair; summarize() keeps the
# max per window. right_join() onto 1..length(win_gr) puts the result
# back in window order and gives NA to a window that holds no SNP.
win_peak <- function(win_gr, gwas) {
  snp_gr <- GRanges(gwas$Chromosome, IRanges(gwas$Position, width = 1))
  h <- findOverlaps(win_gr, snp_gr)
  tibble(w = queryHits(h), neglog10p = gwas$neglog10p[subjectHits(h)]) %>%
    group_by(w) %>%
    summarize(peak = max(neglog10p), .groups = "drop") %>%
    right_join(tibble(w = seq_along(win_gr)), by = "w") %>%
    arrange(w) %>%
    pull(peak)
}

# ============================================================
# 5a. Per-run summary: fish, SNPs, lambda, windows, regions, genes
# ============================================================
# Window counts use the flags of the all-window file:
#   tested  outlier99_raw is TRUE or FALSE (NA = window not tested)
#   long    longer than the length cap of the per-run script
#   usable  tested and not long: the denominator of the percentage of
#           significant windows, so numerator and denominator come from
#           the same set of windows
# cap_kb is the longest window kept, so every excluded window is longer.
#
# The argument is called this_run, not run, on purpose. Inside filter(),
# a bare name is looked up among the columns first, so
# filter(cand_regions, run == run) would compare the column with itself
# and keep every row. A different name avoids the clash.
run_summary <- function(this_run, gwas, win_gr, reg, allw) {
  peaks  <- win_peak(win_gr, gwas)
  tested <- !is.na(allw$outlier99_raw)
  usable <- tested & !allw$long_window
  cr <- filter(cand_regions, run == this_run)
  cg <- filter(cand_genes,   run == this_run)
  tibble(
    run               = this_run,
    n_fish            = count_fish(pheno_files[this_run]),
    n_snps            = nrow(gwas),                           # sites with LRT >= 0 (ANGSD -999 excluded)
    lambda_snp        = median(gwas$LRT) / qchisq(0.5, df = 1),
    n_tested          = sum(tested),
    n_long            = sum(tested & allw$long_window),
    cap_kb            = max(allw$win_len[!allw$long_window]) / 1e3,
    n_usable          = sum(usable),
    n_sig_windows     = length(win_gr),
    pct_sig_windows   = 100 * length(win_gr) / sum(usable),   # compare with 1% expected by chance
    sig_Mb            = bp(reg) / 1e6,
    pct_genome        = 100 * bp(reg) / genome_bp,
    n_weak_windows    = sum(peaks < weak_peak, na.rm = TRUE),
    n_cand_regions    = nrow(cr),
    cand_Mb           = sum(cr$width_kb) / 1e3,
    median_region_kb  = safe(median, cr$width_kb),
    max_region_kb     = safe(max, cr$width_kb),
    n_cand_genes      = n_distinct(cg$gene_name),
    max_intergenic_kb = safe(max, cg$distance_bp[cg$case == "intergenic"]) / 1e3
  )
}

run_stats <- bind_rows(
  run_summary("summer", gwas_s, win_s_gr, reg_s, allwin_s),
  run_summary("fall",   gwas_f, win_f_gr, reg_f, allwin_f)
)
fwrite(run_stats, "res/gwas_results_run_stats.tsv", sep = "\t")

# Window-level lambda and its permutation range, from qq_windows.R
lambda_tab <- if (file.exists(lambda_file)) fread(lambda_file) else NULL

# ============================================================
# 5b. Candidate regions: counts, strongest regions, SI table
# ============================================================

# class is a factor (set in section 3), so .drop = FALSE also lists
# classes with zero regions instead of leaving them out
class_counts <- count(cand_regions, run, class, .drop = FALSE)
case_counts  <- count(cand_regions, run, case)

# slice_max() keeps the n_top rows with the largest max_neglog10p in each
# group; with_ties = FALSE keeps exactly n_top even if two peaks are equal
top_regions <- cand_regions %>%
  group_by(run) %>%
  slice_max(max_neglog10p, n = n_top, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(run, chr = chr_num(Chromosome),
            start_Mb = round(start / 1e6, 3), end_Mb = round(end / 1e6, 3),
            width_kb = round(width_kb, 1),
            peak_neglog10p = round(max_neglog10p, 1),
            class, case, genes)

supp_table <- cand_regions %>%
  transmute(Ecotype = if_else(run == "summer", "Summer run", "Fall run"),
            Chromosome = chr_num(Chromosome),
            Start = start, End = end,
            `Width (kb)` = round(width_kb, 1),
            `Significant windows` = n_windows,
            SNPs = n_snps,
            `Lowest window median P` = signif(min_window_median_p, 2),
            `Strongest SNP` = top_pos,
            `Max -log10 P` = round(max_neglog10p, 2),
            `Second SNP` = second_pos,
            `Shared with other ecotype` = shared,
            `In FST outlier region` = fst_outlier,
            `Core interval` = case,
            `Candidate genes` = genes) %>%
  arrange(Ecotype, Chromosome, Start)
fwrite(supp_table, "res/supp_table_gwas_candidate_regions.tsv", sep = "\t")

# ============================================================
# 5c. Chromosome 36 inversion seen from the within-ecotype GWAS
# ============================================================
# Two intervals [a1, a2] and [b1, b2] overlap when a1 <= b2 and a2 >= b1.
# The gap between them is max(0, b1 - a2, a1 - b2): 0 when they overlap,
# otherwise the distance in bp (same formula as in genes_for_region()).

inv_summary <- function(this_run, gwas, win) {
  in_inv <- gwas$Chromosome == inv_chr & between(gwas$Position, inv_range[1], inv_range[2])
  win_in <- win$Chromosome == inv_chr & win$end_pos >= inv_range[1] & win$start_pos <= inv_range[2]
  cr <- cand_regions %>%
    filter(run == this_run, Chromosome == inv_chr) %>%
    mutate(gap_bp = pmax(0, start - inv_range[2], inv_range[1] - end))
  nearest <- slice_min(cr, gap_bp, n = 1, with_ties = FALSE)
  tibble(
    run                      = this_run,
    snps_in_inv              = sum(in_inv),
    snps_per_10kb_inv        = sum(in_inv) / (diff(inv_range) / 1e4),
    snps_per_10kb_genome     = nrow(gwas) / (genome_bp / 1e4),
    median_p_inv             = safe(median, gwas$P[in_inv]),     # ~0.5 under the null
    min_p_inv                = safe(min, gwas$P[in_inv]),
    expected_min_p           = 1 / (sum(in_inv) + 1),           # expected lowest P among that many null tests
    sig_windows_in_inv       = sum(win_in),
    cand_regions_in_inv      = sum(cr$gap_bp == 0),
    nearest_region_gap_kb    = safe(min, cr$gap_bp) / 1e3,
    nearest_region_Mb        = safe(identity, (nearest$start + nearest$end) / 2e6),
    nearest_region_neglog10p = safe(identity, nearest$max_neglog10p),
    nearest_region_genes     = if (nrow(nearest) > 0) nearest$genes else NA_character_
  )
}

inv_stats <- bind_rows(inv_summary("summer", gwas_s, win_s),
                       inv_summary("fall",   gwas_f, win_f))
fwrite(inv_stats, "res/gwas_results_inversion.tsv", sep = "\t")

# ============================================================
# 5d. GO enrichment
# ============================================================
# as_tibble(): fread() returns a data.table, whose print() ignores n = Inf
go_tab <- if (file.exists(go_file)) as_tibble(fread(go_file)) else NULL

# ============================================================
# 5e. Write the report
# ============================================================
# sink(file, split = TRUE) sends everything printed to the console to
# the file as well. on.exit(sink()) runs when the function ends, even if
# it stops on an error, so the console is never left redirected.
#
# options(width = 200) lets wide tables print on one line instead of
# wrapping at 80 characters. options() returns the old value, which
# on.exit() restores at the end, so the change stays inside this function.
write_report <- function(file) {
  sink(file, split = TRUE)
  on.exit(sink(), add = TRUE)
  old_opt <- options(width = 200)
  on.exit(options(old_opt), add = TRUE)
  
  hdr <- function(x) cat("\n", x, "\n", strrep("-", nchar(x)), "\n", sep = "")
  
  cat("Length-at-maturity GWAS: numbers for the Results\n")
  cat("Generated", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
  
  hdr("1. Samples, SNPs and calibration")
  for (i in seq_len(nrow(run_stats))) {
    r <- run_stats[i, ]
    cat(sprintf("%s run: n = %s fish; %s SNPs tested; lambda_SNP = %.3f\n",
                r$run, ifelse(is.na(r$n_fish), "? (set pheno_files)", fmt_n(r$n_fish)),
                fmt_n(r$n_snps), r$lambda_snp))
  }
  if (!is.null(lambda_tab)) {
    cat(sprintf("\nWindow-level lambda (from %s):\n", lambda_file))
    print(lambda_tab)
  } else {
    cat(sprintf("\n(no %s - run qq_windows.R for window-level lambda)\n", lambda_file))
  }
  
  hdr("2. Significant windows (empirical P <= 0.01)")
  for (i in seq_len(nrow(run_stats))) {
    r <- run_stats[i, ]
    cat(sprintf("%s run: %s of %s usable windows significant (%.2f%%; 1%% expected by chance), covering %.2f Mb (%.2f%% of the genome analysed)\n",
                r$run, fmt_n(r$n_sig_windows), fmt_n(r$n_usable), r$pct_sig_windows,
                r$sig_Mb, r$pct_genome))
    cat(sprintf("   usable = tested and no longer than %.1f kb; %s of %s tested windows were longer and excluded\n",
                r$cap_kb, fmt_n(r$n_long), fmt_n(r$n_tested)))
    cat(sprintf("   significant windows whose strongest SNP has -log10 P < %g: %s (%.1f%%)\n",
                weak_peak, fmt_n(r$n_weak_windows), 100 * r$n_weak_windows / r$n_sig_windows))
  }
  
  hdr(sprintf("3. Overlap tests (circular shifts, n = %s)", fmt_n(n_shift)))
  for (i in seq_len(nrow(test_tab))) {
    t <- test_tab[i, ]
    cat(sprintf("%s: %.1f kb shared vs %.1f kb expected (%.2f-fold, %s)",
                t$comparison, t$shared_kb, t$expected_kb, t$fold, fmt_p(t$emp_p, n_shift)))
    if (!is.na(t$n_windows_a))
      cat(sprintf("; %s of %s windows of the first set overlap the second",
                  fmt_n(t$n_windows_a_overlapping_b), fmt_n(t$n_windows_a)))
    cat("\n")
  }
  
  hdr(sprintf("4. Candidate regions (>= %d contiguous significant windows, one with median P < %g)",
              min_adjacent, max_median_p))
  for (i in seq_len(nrow(run_stats))) {
    r <- run_stats[i, ]
    cat(sprintf("%s run: %d regions, %.2f Mb in total (median %.1f kb, max %.1f kb); %d candidate genes\n",
                r$run, r$n_cand_regions, r$cand_Mb, r$median_region_kb, r$max_region_kb, r$n_cand_genes))
    if (!is.na(r$max_intergenic_kb))
      cat(sprintf("   farthest candidate gene of an intergenic region: %.1f kb from the core interval\n",
                  r$max_intergenic_kb))
  }
  cat("\nRegions per class (shared = overlaps significant windows of the other run):\n")
  print(class_counts, n = Inf)
  cat("\nCore interval relative to genes:\n")
  print(case_counts, n = Inf)
  cat(sprintf("\nStrongest %d regions per run:\n", n_top))
  print(top_regions, n = Inf, width = Inf)
  
  hdr(sprintf("5. Chromosome %d inversion (%s:%s-%s Mb)", inv_label, inv_chr,
              inv_range[1] / 1e6, inv_range[2] / 1e6))
  for (i in seq_len(nrow(inv_stats))) {
    v <- inv_stats[i, ]
    cat(sprintf("%s run: %s SNPs tested in the inversion (%.1f per 10 kb vs %.1f genome-wide)\n",
                v$run, fmt_n(v$snps_in_inv), v$snps_per_10kb_inv, v$snps_per_10kb_genome))
    cat(sprintf("   median P = %.2f; lowest P = %.1e (expected lowest P under the null ~ %.1e)\n",
                v$median_p_inv, v$min_p_inv, v$expected_min_p))
    cat(sprintf("   significant windows in the inversion: %d; candidate regions overlapping it: %d\n",
                v$sig_windows_in_inv, v$cand_regions_in_inv))
    if (!is.na(v$nearest_region_gap_kb))
      cat(sprintf("   nearest candidate region on that chromosome: %.0f kb away (centre %.2f Mb, max -log10 P = %.1f; genes: %s)\n",
                  v$nearest_region_gap_kb, v$nearest_region_Mb, v$nearest_region_neglog10p,
                  v$nearest_region_genes))
  }
  
  hdr("6. GO enrichment")
  if (is.null(go_tab)) {
    cat(sprintf("(no %s - run 04_gwas_hit_GO_analaysi.R)\n", go_file))
  } else {
    cat(sprintf("%d terms tested; %d with FDR < %g\n",
                nrow(go_tab), sum(go_tab$fdr < go_fdr), go_fdr))
    go_tab %>%
      filter(fdr < go_fdr) %>%
      dplyr::select(any_of(c("go_id", "term", "n_regions", "expected", "fold", "emp_p", "fdr"))) %>%
      mutate(across(any_of(c("expected", "fold")), ~ round(.x, 2)),
             across(any_of(c("emp_p", "fdr")), ~ signif(.x, 2))) %>%
      print(n = Inf, width = Inf)
  }
  invisible(NULL)
}

write_report(report_file)
cat("\nReport written to", report_file, "\n")