#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(magrittr)
library(data.table)
library(patchwork)
library(GenomicRanges)
library(biomaRt)
library(stringdist)
library(ggrepel)

#------------#
#--- DATA ---#
#------------#

#Sliding window

gwas <- fread("data/gwas_length/gwas_length.lrt0.gz")

# Calculate continuous positions across chromosomes
gwas <- gwas %>%
  mutate(
    P_recalc = pchisq(LRT, df = 1, lower.tail = FALSE),
    neglog10p = -log10(P_recalc),
    Chromosome = factor(Chromosome, levels = unique(Chromosome))
  )


#------------------------------------------#
#--- Equal-spaced chromosome coordinates --#
#------------------------------------------#

chrom_spacing <- 1e7

gwas_equal <- gwas %>%
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

chr_map <- gwas %>%
  distinct(Chromosome) %>%
  arrange(Chromosome) %>%
  mutate(chr_index = row_number() - 1)

gwas <- gwas %>%
  group_by(Chromosome) %>%
  mutate(pos_scaled = Position / max(Position)) %>%
  ungroup() %>%
  left_join(chr_map, by = "Chromosome") %>%
  mutate(pos_equal = chr_index + pos_scaled)

axis_equal <- gwas_equal %>%
  group_by(Chromosome) %>%
  summarize(center = mean(pos_equal))


# Calculate threshold
thr <- quantile(gwas$neglog10p, 0.9999)


##############################################################3
#FInd gene in outlier windows


gene_annotation <- fread("data/genome_annotation_table.tsv")

gene_annotation %<>% dplyr::select(ScaffoldName,FromPosition,ToPosition,GeneName)

genes_gr <- GRanges(
  seqnames = gene_annotation$ScaffoldName,
  ranges = IRanges(start = gene_annotation$FromPosition, end = gene_annotation$ToPosition),
  gene_name = gene_annotation$GeneName
)

#Merge outliers in windows
outlier_gwas <- gwas %>% filter(neglog10p > thr)

outlier_snps_gr <- GRanges(
  seqnames = outlier_gwas$Chromosome,
  ranges = IRanges(start = outlier_gwas$Position, width = 1),
  neglog10p = outlier_gwas$neglog10p,
  pval = outlier_gwas$neglog10p)

merge_distance = 100000

windows_reduced_gr <- reduce(outlier_snps_gr, 
                             min.gapwidth = merge_distance + 1,
                             with.revmap = TRUE)


#Find overlap with genes

window_stats <- mcols(windows_reduced_gr)$revmap %>%
  lapply(function(idx) {
    snps_in_window <- outlier_snps_gr[idx]
    tibble(
      n_snps = length(snps_in_window),
      max_neglog10p = max(mcols(snps_in_window)$neglog10p),
      min_pval = min(mcols(snps_in_window)$pval),
      mean_neglog10p = mean(mcols(snps_in_window)$neglog10p)
    )
  }) %>%
  bind_rows()

mcols(windows_reduced_gr) <- cbind(
  mcols(windows_reduced_gr),
  window_stats
)

#Find overlap with genes
overlaps <- findOverlaps(windows_reduced_gr, genes_gr)

# Extract overlapping pairs
overlap_df <- data.frame(
  window_idx = queryHits(overlaps),
  gene_idx = subjectHits(overlaps)
) %>%
  mutate(
    # Window information
    chromosome = as.character(seqnames(windows_reduced_gr)[window_idx]),
    window_start = start(windows_reduced_gr)[window_idx],
    window_end = end(windows_reduced_gr)[window_idx],
    window_name = windows_reduced_gr$window_name[window_idx],
    n_snps = windows_reduced_gr$n_snps[window_idx],
    max_neglog10p = windows_reduced_gr$max_neglog10p[window_idx],
    min_pval = windows_reduced_gr$min_pval[window_idx],
    
    # Gene information
    gene_name = genes_gr$gene_name[gene_idx],
    gene_start = start(genes_gr)[gene_idx],
    gene_end = end(genes_gr)[gene_idx],
    
    # Overlap information
    overlap_start = pmax(window_start, gene_start),
    overlap_end = pmin(window_end, gene_end),
    overlap_width = overlap_end - overlap_start + 1
  ) %>%
  arrange(desc(max_neglog10p), chromosome, window_start)


overlap_df %<>% filter(overlap_width > 100) %>% 
  distinct(window_start,window_end,gene_name,.keep_all = T) 


######################

mart <- useMart(
  biomart = "ensembl",
  dataset = "drerio_gene_ensembl"  # change if needed
)

#get ensembl gene names
mart_genes <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "description",
    "hsapiens_homolog_associated_gene_name"
  ),
  mart = mart
)


clean_string <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all("\\(.*?\\)", "") %>%   # remove parenthesis
    str_replace_all("[^a-z0-9 ]", " ") %>% # remove punctuation
    str_squish()
}

outlier_annotation <- overlap_df %>% 
  dplyr::select(chromosome, gene_name, gene_start, gene_end) %>% 
  mutate(query_desc = clean_string(gene_name))

mart_genes <- mart_genes %>%
  mutate(target_desc = clean_string(description))

mart_genes %<>% filter(!str_detect(target_desc,"prothymosin|karyopherin"))


fuzzy_map <- outlier_annotation %>%
  distinct(gene_name, query_desc) %>%
  mutate(
    match = map(
      query_desc,
      ~ mart_genes %>%
        mutate(dist = stringdist(.x, target_desc, method = "jw")) %>%
        arrange(dist) %>%
        slice_head(n = 1)
    )
  ) %>%
  unnest(match)

fuzzy_map %<>% filter(!description == "")

####################################################
# Add gene names to gwas

# For each reduced window, get the SNP with max signal
window_peaks <- mcols(windows_reduced_gr)$revmap %>%
  lapply(function(idx) {
    snps <- outlier_gwas[idx, ]
    snps %>%
      slice_max(neglog10p, n = 1)
  }) %>%
  bind_rows() %>%
  mutate(window_id = row_number())

gene_labels <- overlap_df %>%
  mutate(window_idx = window_idx) %>%
  left_join(
    window_peaks %>% mutate(window_idx = row_number()),
    by = "window_idx"
  ) %>% left_join(fuzzy_map, by = "gene_name")



gene_labels_plot <- gene_labels %>%
  group_by(window_idx) %>%
  slice_max(max_neglog10p, n = 1) %>%
  ungroup()


gene_labels_plot$external_gene_name

chr_map <- gwas %>%
  distinct(Chromosome) %>%
  arrange(Chromosome) %>%
  mutate(chr_index = row_number() - 1)

gene_labels_plot <- gene_labels_plot %>%
  left_join(
    gwas %>% dplyr::select(Chromosome, Position, pos_equal),
    by = c("Chromosome", "Position")
  )

# gene_labels_plot <- gene_labels_plot %>%
#   left_join(chr_map, by = "Chromosome") %>%
#   group_by(Chromosome) %>%
#   mutate(pos_scaled = Position / max(Position)) %>%
#   ungroup() %>%
#   mutate(pos_equal = chr_index + pos_scaled)

gwas <- gwas %>%
  group_by(Chromosome) %>%
  mutate(pos_scaled = Position / max(Position)) %>%
  ungroup()

gwas <- gwas %>%
  left_join(chr_map, by = "Chromosome") %>%
  mutate(pos_equal = chr_index + pos_scaled)

axis_df <- gwas %>%
  group_by(Chromosome) %>%
  dplyr::summarize(center = mean(pos_equal), .groups = "drop")

gwas_plot <- ggplot(
  gwas,
  aes(x = pos_equal, y = neglog10p, color = Chromosome)
) +
  geom_point(size = 0.5, alpha = 0.9) +
  geom_hline(
    yintercept = thr,
    linetype = "dashed",
    color = "red4",
    linewidth = 0.7
  ) +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = seq(1,38),
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, max(gwas$neglog10p, na.rm = TRUE) * 1.05),
    expand = c(0, 0)
  ) +
  scale_color_manual(
    values = rep(c("grey50", "grey70"),
                 length(unique(gwas$Chromosome)))
  ) +
  labs(
    x = "Chromosome",
    y = expression(-log10
    )
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )


gwas_plot_genes <- gwas_plot +
  geom_text_repel(
    data = gene_labels_plot,
    aes(
      x = pos_equal,
      y = neglog10p,
      label = external_gene_name
    ),
    inherit.aes = FALSE,
    size = 4,
    fontface = "italic",
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey40",
    segment.size = 0.4,
    min.segment.length = 0,
    max.overlaps = Inf
  )


ggsave("figures/04_gwas_length_gene_names.png",gwas_plot,width = 30, height = 10, unit = "cm", dpi = 900)
