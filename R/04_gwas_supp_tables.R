# ============================================================
# Supplementary tables for the within-ecotype GWAS of length at
# reproduction
#
#   Supplementary Table 17  summary per ecotype: fish, SNPs, lambda,
#                           windows, candidate regions, chr36 inversion
#   Supplementary Table 18  overlap tests (circular-shift null)
#   Supplementary Table 19  candidate regions and candidate genes
#
# Inputs are the files written by 04_gwas_summer_perm.R,
# 04_gwas_fall_perm.R and 04_gwas_withineco_fst_overlap_gene_label.R
# (run those first). Each table replaces the sheet of the same name in
# the supplementary workbook, at the same position, and is also saved
# as a TSV in res/.
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(openxlsx)

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------
supp_file <- "res/RB_Chapter_3_Supplementary_Tables.xlsx"

run_stats_file   <- "res/gwas_results_run_stats.tsv"
inversion_file   <- "res/gwas_results_inversion.tsv"
overlap_file     <- "res/gwas_summer_fall_fst_overlap_tests.tsv"
regions_file     <- "res/gwas_length_candidate_regions.tsv"
genes_file       <- "res/gwas_length_candidate_region_genes.tsv"
chr_number_file  <- "res/supp_table_gwas_candidate_regions.tsv"   # gives the chromosome numbers used on the figures

sheet_summary <- "Supplementary Table 17"
sheet_overlap <- "Supplementary Table 18"
sheet_regions <- "Supplementary Table 19"

ecotype_label <- c(summer = "Summer run (kâcikâsikumekw)",
                   fall   = "Fall run (nûtamesânîw)")

# Same values as in 04_gwas_withineco_fst_overlap_gene_label.R (labels only)
inv_label    <- 36
inv_range_mb <- c(34.8, 35.6)
n_shift      <- 1000

comparison_label <- c(
  "fall x summer"         = "Fall-run vs summer-run significant windows",
  "fall x Fst"            = "Fall-run significant windows vs FST outlier regions",
  "summer x Fst"          = "Summer-run significant windows vs FST outlier regions",
  "(fall & summer) x Fst" = "Windows significant in both ecotypes vs FST outlier regions"
)

case_label <- c(within_gene  = "Within gene",
                partly_genic = "Partly genic",
                across_genes = "Across genes",
                intergenic   = "Intergenic")

out_summary <- "res/supp_gwas_summary.tsv"
out_overlap <- "res/supp_gwas_overlap_tests.tsv"
out_regions <- "res/supp_gwas_candidate_regions_genes.tsv"

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
read_tsv <- function(f) {
  if (!file.exists(f)) stop(f, " not found: run the GWAS scripts first")
  as_tibble(fread(f))
}

yes_no <- function(x) if_else(as.logical(x), "Yes", "No")

drop_evidence_tags <- function(x) str_squish(str_remove_all(x, "\\{ECO:[^}]*\\}"))

# Replace a sheet by a new one with the same name, keeping its position
replace_sheet <- function(wb, sheet, data) {
  if (!sheet %in% names(wb)) stop("Sheet '", sheet, "' not found in ", supp_file)
  order_names <- names(wb)[worksheetOrder(wb)]
  removeWorksheet(wb, sheet)
  addWorksheet(wb, sheet)
  writeData(wb, sheet, data)
  worksheetOrder(wb) <- match(order_names, names(wb))
}

# ============================================================
# Supplementary Table 17: summary per ecotype
# ============================================================
# One row per statistic; `digits` gives decimal places, `sig_digits`
# significant digits, NA = written as is.

inv_section <- sprintf("Chromosome %d inversion (%.1f–%.1f Mb)",
                       inv_label, inv_range_mb[1], inv_range_mb[2])

summary_spec <- tribble(
  ~Section,                     ~Statistic,                                           ~column,               ~digits, ~sig_digits,
  "Samples and SNPs",           "Fish with a length phenotype",                       "n_fish",              NA,     NA,
  "Samples and SNPs",           "SNPs tested",                                        "n_snps",              NA,     NA,
  "Samples and SNPs",           "Genomic inflation factor, λSNP",                     "lambda_snp",          3,      NA,
  "Windows (10 SNPs, step 2)",  "Windows tested",                                     "n_tested",            NA,     NA,
  "Windows (10 SNPs, step 2)",  "Length cap (kb)",                                    "cap_kb",              1,      NA,
  "Windows (10 SNPs, step 2)",  "Windows longer than the cap (excluded)",             "n_long",              NA,     NA,
  "Windows (10 SNPs, step 2)",  "Windows retained",                                   "n_usable",            NA,     NA,
  "Windows (10 SNPs, step 2)",  "Significant windows (empirical P ≤ 0.01)",           "n_sig_windows",       NA,     NA,
  "Windows (10 SNPs, step 2)",  "Significant windows, % of retained (1% expected)",   "pct_sig_windows",     2,      NA,
  "Windows (10 SNPs, step 2)",  "Genome covered by significant windows (Mb)",         "sig_Mb",              2,      NA,
  "Windows (10 SNPs, step 2)",  "Genome covered by significant windows (%)",          "pct_genome",          2,      NA,
  "Candidate regions",          "Candidate regions",                                  "n_cand_regions",      NA,     NA,
  "Candidate regions",          "Total length (Mb)",                                  "cand_Mb",             2,      NA,
  "Candidate regions",          "Median length (kb)",                                 "median_region_kb",    1,      NA,
  "Candidate regions",          "Maximum length (kb)",                                "max_region_kb",       1,      NA,
  "Candidate regions",          "Candidate genes",                                    "n_cand_genes",        NA,     NA,
  inv_section,                  "SNPs tested in the inversion",                       "snps_in_inv",         NA,     NA,
  inv_section,                  "SNPs per 10 kb, inversion",                          "snps_per_10kb_inv",   2,      NA,
  inv_section,                  "SNPs per 10 kb, genome-wide",                        "snps_per_10kb_genome", 2,     NA,
  inv_section,                  "SNP density in the inversion, % of genome-wide",     "density_pct",         1,      NA,
  inv_section,                  "Median P in the inversion",                          "median_p_inv",        NA,     2,
  inv_section,                  "Lowest P in the inversion",                          "min_p_inv",           NA,     2,
  inv_section,                  "Expected lowest P under the null",                   "expected_min_p",      NA,     2,
  inv_section,                  "Significant windows in the inversion",               "sig_windows_in_inv",  NA,     NA,
  inv_section,                  "Candidate regions overlapping the inversion",        "cand_regions_in_inv", NA,     NA
)

run_stats <- read_tsv(run_stats_file)
inv_stats <- read_tsv(inversion_file)

summary_values <- run_stats %>%
  left_join(inv_stats, by = "run") %>%
  mutate(density_pct = 100 * snps_per_10kb_inv / snps_per_10kb_genome) %>%
  select(run, all_of(summary_spec$column)) %>%
  mutate(across(-run, as.numeric)) %>%
  pivot_longer(-run, names_to = "column") %>%
  pivot_wider(names_from = run, values_from = value)

supp_summary <- summary_spec %>%
  left_join(summary_values, by = "column") %>%
  mutate(across(all_of(names(ecotype_label)),
                ~ case_when(!is.na(digits)     ~ round(.x, digits),
                            !is.na(sig_digits) ~ signif(.x, sig_digits),
                            TRUE               ~ .x))) %>%
  select(Section, Statistic, all_of(names(ecotype_label))) %>%
  rename(!!!setNames(names(ecotype_label), ecotype_label))

# ============================================================
# Supplementary Table 18: overlap tests
# ============================================================

overlap <- read_tsv(overlap_file)

supp_overlap <- overlap %>%
  transmute(
    Comparison                           = coalesce(comparison_label[comparison], comparison),
    `Set A (Mb)`                         = round(a_Mb, 2),
    `Set B, shifted (Mb)`                = round(b_Mb, 2),
    `Shared (kb)`                        = round(shared_kb, 1),
    `Expected shared (kb)`               = round(expected_kb, 1),
    `Fold enrichment`                    = round(fold, 2),
    `Empirical P`                        = signif(emp_p, 2),
    `Windows in set A`                   = n_windows_a,
    `Windows in set A overlapping set B` = n_windows_a_overlapping_b
  )

# ============================================================
# Supplementary Table 19: candidate regions and genes
# ============================================================

regions <- read_tsv(regions_file) %>% mutate(Chromosome = as.character(Chromosome))
genes   <- read_tsv(genes_file)

# Chromosome number (as on the figures) for each accession, taken from
# the region table written by the overlap script, which has both
chr_number <- read_tsv(chr_number_file) %>%
  transmute(run = if_else(Ecotype == "Summer run", "summer", "fall"),
            start = Start, end = End, chr_num = Chromosome) %>%
  inner_join(regions %>% select(run, start, end, Chromosome), by = c("run", "start", "end")) %>%
  distinct(Chromosome, chr_num)

if (anyDuplicated(chr_number$Chromosome) || anyDuplicated(chr_number$chr_num))
  stop("Accessions and chromosome numbers do not map one to one")
if (!all(regions$Chromosome %in% chr_number$Chromosome))
  stop(chr_number_file, " does not match ", regions_file, ": rerun the overlap script")

region_genes <- genes %>%
  arrange(run, region_id, gene_start) %>%
  group_by(run, region_id) %>%
  summarize(gene_names     = paste(drop_evidence_tags(gene_name), collapse = "; "),
            gene_distances = paste(distance_bp, collapse = "; "),
            .groups = "drop")

# The region file's own `genes` column holds the plot labels (zebrafish
# symbols, truncated names); the table lists the full annotated names
supp_regions <- regions %>%
  select(-any_of(c("genes", "plot_label"))) %>%
  left_join(chr_number, by = "Chromosome") %>%
  left_join(region_genes, by = c("run", "region_id")) %>%
  rename(accession = Chromosome) %>%
  mutate(run = factor(run, levels = names(ecotype_label))) %>%
  arrange(run, chr_num, start) %>%
  transmute(
    Ecotype                                  = unname(ecotype_label[as.character(run)]),
    Region                                   = region_id,
    Chromosome                               = chr_num,
    Accession                                = accession,
    `Start (bp)`                             = start,
    `End (bp)`                               = end,
    `Width (kb)`                             = round(width_kb, 1),
    `Significant windows`                    = n_windows,
    SNPs                                     = n_snps,
    `Lowest window median P`                 = signif(min_window_median_p, 2),
    `Strongest SNP (bp)`                     = top_pos,
    `Strongest SNP -log10 P`                 = round(max_neglog10p, 2),
    `Second strongest SNP (bp)`              = second_pos,
    `Core interval start (bp)`               = core_start,
    `Core interval end (bp)`                 = core_end,
    `Core interval relative to genes`        = unname(case_label[case]),
    `Overlaps significant windows of the other ecotype` = yes_no(shared),
    `Overlaps FST outlier region`            = yes_no(fst_outlier),
    `Candidate genes`                        = gene_names,
    `Distance to core interval (bp)`         = gene_distances
  )

# ============================================================
# Write
# ============================================================

fwrite(supp_summary, out_summary, sep = "\t")
fwrite(supp_overlap, out_overlap, sep = "\t")
fwrite(supp_regions, out_regions, sep = "\t")

wb <- loadWorkbook(supp_file)
replace_sheet(wb, sheet_summary, supp_summary)
replace_sheet(wb, sheet_overlap, supp_overlap)
replace_sheet(wb, sheet_regions, supp_regions)
saveWorkbook(wb, supp_file, overwrite = TRUE)

cat(sprintf("%s: %d statistics | %s: %d tests | %s: %d candidate regions (%d summer, %d fall)\n",
            sheet_summary, nrow(supp_summary), sheet_overlap, nrow(supp_overlap),
            sheet_regions, nrow(supp_regions),
            sum(regions$run == "summer"), sum(regions$run == "fall")))
cat("Written to", supp_file, "and", out_summary, out_overlap, out_regions, "\n")
