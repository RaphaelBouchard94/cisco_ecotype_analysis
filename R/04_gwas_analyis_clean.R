# ============================================================
# Libraries
# ============================================================

library(tidyverse)
library(data.table)
library(GenomicRanges)
library(biomaRt)
library(stringdist)
library(ggrepel)

# ============================================================
# Load & prepare GWAS data
# ============================================================

gwas <- fread("data/gwas_length/gwas_length.lrt0.gz") %>%
  mutate(
    P_recalc   = pchisq(LRT, df = 1, lower.tail = FALSE),
    neglog10p = -log10(P_recalc),
    Chromosome = factor(Chromosome, levels = unique(Chromosome))
  )

# ============================================================
# Equal-width chromosome coordinates (DEFINED ONCE)
# ============================================================

chrom_spacing <- 1e7
chrom_gap     <- 2e6 

chr_map <- gwas %>%
  distinct(Chromosome) %>%
  dplyr::arrange(Chromosome) %>%
  dplyr::mutate(chr_index = row_number() - 1)

gwas <- gwas %>%
  group_by(Chromosome) %>%
  mutate(pos_scaled = Position / max(Position)) %>%
  ungroup() %>%
  left_join(chr_map, by = "Chromosome") %>%
  mutate(pos_equal = chr_index * (chrom_spacing + chrom_gap) +
           pos_scaled * chrom_spacing)

axis_df <- gwas %>%
  group_by(Chromosome) %>%
  summarize(center = mean(pos_equal), .groups = "drop")

# ============================================================
# Significance threshold
# ============================================================

thr <- quantile(gwas$neglog10p, 0.9999, na.rm = TRUE)

# ============================================================
# Gene annotation & GRanges
# ============================================================

gene_annotation <- fread("data/genome_annotation_table.tsv") %>%
  dplyr::select(ScaffoldName, FromPosition, ToPosition, GeneName)

genes_gr <- GRanges(
  seqnames = gene_annotation$ScaffoldName,
  ranges   = IRanges(
    start = gene_annotation$FromPosition,
    end   = gene_annotation$ToPosition
  ),
  gene_name = gene_annotation$GeneName
)

# ============================================================
# Outlier windows
# ============================================================

outlier_gwas <- gwas %>% filter(neglog10p > thr)

outlier_snps_gr <- GRanges(
  seqnames = outlier_gwas$Chromosome,
  ranges   = IRanges(start = outlier_gwas$Position, width = 1),
  neglog10p = outlier_gwas$neglog10p,
  pval      = outlier_gwas$neglog10p
)

windows_reduced_gr <- reduce(
  outlier_snps_gr,
  min.gapwidth = 100000 + 1,
  with.revmap  = TRUE
)

# ============================================================
# Window statistics
# ============================================================

window_stats <- mcols(windows_reduced_gr)$revmap %>%
  lapply(function(idx) {
    snps <- outlier_snps_gr[idx]
    tibble(
      n_snps          = length(snps),
      max_neglog10p   = max(mcols(snps)$neglog10p),
      min_pval        = min(mcols(snps)$pval),
      mean_neglog10p  = mean(mcols(snps)$neglog10p)
    )
  }) %>%
  bind_rows()

mcols(windows_reduced_gr) <- cbind(
  mcols(windows_reduced_gr),
  window_stats
)

# ============================================================
# Overlap windows with genes
# ============================================================

overlaps <- findOverlaps(windows_reduced_gr, genes_gr)

overlap_df <- tibble(
  window_idx = queryHits(overlaps),
  gene_idx   = subjectHits(overlaps)
) %>%
  dplyr::mutate(
    chromosome   = as.character(seqnames(windows_reduced_gr)[window_idx]),
    window_start = start(windows_reduced_gr)[window_idx],
    window_end   = end(windows_reduced_gr)[window_idx],
    n_snps       = windows_reduced_gr$n_snps[window_idx],
    max_neglog10p = windows_reduced_gr$max_neglog10p[window_idx],
    gene_name    = genes_gr$gene_name[gene_idx],
    gene_start   = start(genes_gr)[gene_idx],
    gene_end     = end(genes_gr)[gene_idx],
    overlap_start = pmax(window_start, gene_start),
    overlap_end   = pmin(window_end, gene_end),
    overlap_width = overlap_end - overlap_start + 1
  ) %>%
  filter(overlap_width > 100) %>%
  distinct(window_start, window_end, gene_name, .keep_all = TRUE) %>%
  dplyr::arrange(desc(max_neglog10p))


#####
#Put results in supp table
# 
# wb <- loadWorkbook("res/supplementary_table.xlsx")
# 
# addWorksheet(wb, "GWAS_outlier_genes")
# 
# writeData(
#   wb,
#   sheet = "GWAS_outlier_genes",
#   overlap_df
# )
# 
# saveWorkbook(wb, "res/supplementary_table.xlsx", overwrite = TRUE)
# 



# ============================================================
# Fuzzy gene name mapping (Ensembl)
# ============================================================

clean_string <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all("\\(.*?\\)", "") %>%
    str_replace_all("[^a-z0-9 ]", " ") %>%
    str_squish()
}

mart <- useMart("ensembl", dataset = "drerio_gene_ensembl")

mart_genes <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "description"
  ),
  mart = mart
) %>%
  mutate(target_desc = clean_string(description)) %>%
  filter(!str_detect(target_desc, "prothymosin|karyopherin"))

fuzzy_map <- overlap_df %>%
  dplyr::distinct(gene_name) %>%
  dplyr::mutate(query_desc = clean_string(gene_name)) %>%
  filter(query_desc != "") %>% 
  dplyr::mutate(
    match = purrr::map(
      query_desc,
      ~ mart_genes %>%
        dplyr::mutate(dist = stringdist(.x, target_desc, method = "jw"))  %>%
        dplyr::arrange(dist) %>%
        slice_head(n = 1)
    )
  ) %>%
  unnest(match)


# Gene labels (use GWAS coordinates — no recomputation!)
# ============================================================

window_peaks <- mcols(windows_reduced_gr)$revmap %>%
  lapply(function(idx) {
    outlier_gwas[idx, ] %>% slice_max(neglog10p, n = 1)
  }) %>%
  bind_rows() %>%
  dplyr::mutate(window_idx = row_number())

gene_labels_plot <- overlap_df %>%
  left_join(window_peaks, by = "window_idx") %>%
  left_join(fuzzy_map, by = "gene_name") %>%
  group_by(window_idx) %>%
  slice_max(max_neglog10p, n = 1) %>%
  ungroup() %>%
  left_join(
    gwas %>% dplyr::select(Chromosome, Position, pos_equal),
    by = c("Chromosome", "Position")
  )

# ============================================================
# Plot
# ============================================================

gwas_plot <- ggplot(gwas, aes(pos_equal, neglog10p, color = Chromosome)) +
  geom_point(size = 0.4, alpha = 0.9) +
  geom_hline(yintercept = thr, linetype = "dashed", color = "red4") +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = seq(1,38),
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, max(gwas$neglog10p, na.rm = TRUE) * 1.05),) +
  scale_color_manual(values = rep(c("grey50", "grey70"),
                                  length(unique(gwas$Chromosome)))) +
  labs(x = "Chromosome", y = expression(-log[10](P))) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )

gwas_plot_label <- gwas_plot +
  geom_text_repel(
    data = gene_labels_plot,
    aes(pos_equal.x, neglog10p, label = external_gene_name),
    inherit.aes = FALSE, 
    color = "#0B1F3A",
    size = 4,
    fontface = "italic",
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "#5B6B7A",
    segment.size = 0.4,
    min.segment.length = 0,
    max.overlaps = Inf
  )


ggsave(
  "figures/04_gwas_length_gene_names.png", gwas_plot_label,
  width = 30, height = 7, units = "cm", dpi = 900
)

####################################3
#Compare results with Fst outliers

gwas_fst <- overlap_df %>% filter(gene_name != "-") %>% 
  left_join(results, by = "gene_name")











