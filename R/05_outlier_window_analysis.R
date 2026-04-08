#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(magrittr)
library(data.table)
library(patchwork)
library(GenomicRanges)
library(windowscanr)
library(openxlsx)
library(biomaRt)
library(stringdist)

#------------#
#--- DATA ---#
#------------#


notsum_rup <- fread("data/genome_scan/notsummer_rup_20kb.outlier.fst")

# pi_window_nots <- fread("data/genome_scan/notsummer.thetaswindow20kb.gz.pestPG")
# 
# pi_window_rup <- fread("data/genome_scan/rup.thetaswindow20kb.gz.pestPG")

pi_window_nots <- fread("data/genome_scan/chr36_outlier_notsummer_merged.pestPG")
pi_window_rup <- fread("data/genome_scan/chr36_outlier_rup_merged.pestPG")

dxy <- fread("data/genome_scan/chr36_outlier_dxy_results.txt")

################################################################################
#Find genes in outlier windows

outlier_windows <- fread("data/genome_scan/notsum_rup_merged_top_candidate.tsv")

gene_annotation <- fread("data/genome_annotation_table.tsv")

gene_annotation %<>% dplyr::select(ScaffoldName,FromPosition,ToPosition,GeneName)


# Convert to GRangeso bjects
outlier_gr <- GRanges(
  seqnames = outlier_windows$chrom,
  ranges = IRanges(start = outlier_windows$window_start, end = outlier_windows$window_end))

genes_gr <- GRanges(
  seqnames = gene_annotation$ScaffoldName,
  ranges = IRanges(start = gene_annotation$FromPosition, end = gene_annotation$ToPosition),
  gene_name = gene_annotation$GeneName
)

# Find overlaps
overlaps <- findOverlaps(outlier_gr, genes_gr)

# Results

results <- data.frame(
  outlier_chr = seqnames(outlier_gr[queryHits(overlaps)]),
  outlier_start = start(outlier_gr[queryHits(overlaps)]),
  outlier_end = end(outlier_gr[queryHits(overlaps)]),
  gene_name = mcols(genes_gr)$gene_name[subjectHits(overlaps)],
  gene_start = start(genes_gr[subjectHits(overlaps)]),
  gene_end = end(genes_gr[subjectHits(overlaps)])
)


results %<>% distinct(outlier_chr,outlier_end,gene_name,.keep_all = T) 

head(results)

unique_genes <- unique(results$gene_name)

results %>% filter(outlier_chr == "CM078251.1")

outlier_windows %>% filter(chrom == "CM078251.1")

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

listDatasets(mart) %>%
  dplyr::filter(grepl("salmo|oncorhynchus|salvelinus|coregon|gasterosteus|drerio", dataset))


mart <- useMart("ensembl")

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

fuzzy_map <- results %>%
  distinct(gene_name) %>%
  mutate(query_desc = clean_string(gene_name)) %>%
  filter(query_desc != "") %>% 
  mutate(
    match = map(
      query_desc,
      ~ mart_genes %>%
        mutate(dist = stringdist(.x, target_desc, method = "jw"))  %>%
        arrange(dist) %>%
        slice_head(n = 1)
    )
  ) %>%
  unnest(match)



#####
#Put results in supp table

wb <- loadWorkbook("res/supplementary_table.xlsx")

addWorksheet(wb, "Fst_outwin_merged_genes")

writeData(
  wb,
  sheet = "Fst_outwin_merged_genes",
  results
)

saveWorkbook(wb, "res/supplementary_table.xlsx", overwrite = TRUE)


#############################################################
#Outlier FST chromosome 34 = CM078251.1

targ_win <- c(34800000,35619999)
targ_win_plus <- c(34800000-1000000,35619999+1000000)


chrom_34_fst <- notsum_rup %>% filter(chrom == "CM078251.1")

chrom_34_fst %>% filter(pos > targ_win[1] & pos < targ_win[2]) %>% count(fst > 0.9)


chrom_34_fst %>% mutate(mid = (window_end + window_start)/2) %>% 
  group_by(mid) %>% summarize(fst_avg = mean(fst))

fst_plot <- ggplot(chrom_34_fst %>% filter(pos > targ_win_plus[1] & pos < targ_win_plus[2]), 
                   aes(x = pos , y = fst)) +
  geom_point(alpha = 0.5, color = "black",size = 1) +
  geom_smooth(color = "grey50", span = 0.01)+
  geom_rect(data=outlier_windows %>% filter(chrom == "CM078251.1",window_start==targ_win[1]), 
            inherit.aes=FALSE, 
            aes(xmin=window_start, xmax=window_end,
                ymin=0,ymax=1, 
                group=chrom), 
            color="transparent", 
            fill="indianred", 
            alpha=0.2)+
  geom_rect(xmin=35059641, xmax=35060636+2000, 
            ymin=0,ymax=1, 
            color="transparent",
            fill="skyblue4",
            alpha=0.2)+
  geom_rect(xmin=35589506, xmax=35618151,
            ymin=0,ymax=1, 
            color="transparent",
            fill="skyblue4",
            alpha=0.2)+
  labs(x = "",
       y = "Fst") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        panel.spacing = unit(0,'lines'),
        axis.text.x = element_blank(),
        panel.grid = element_blank())

##########################################33
#Tajima
taj_nots <- pi_window_nots %>% filter(Chr == "CM078251.1")

taj_rup <- pi_window_rup %>% filter(Chr == "CM078251.1")


taj <- cbind(taj_nots %>% dplyr::select(WinCenter,Tajima,tP),taj_rup %>% dplyr::select(Tajima,tP))

colnames(taj) <- c("WinCenter","taj_kuu","pi_kuu","taj_rup","pi_rup")


taj_plot <- ggplot(taj) +
  geom_smooth(aes(x = WinCenter*2 , y = taj_kuu),
              color = "#E69F00", span = 0.1) +
  geom_point(aes(x = WinCenter*2 , y = taj_kuu),
             color = "#E69F00", alpha = 0.5, size = 1) +
  geom_smooth(aes(x = WinCenter*2 , y = taj_rup),
              color = "skyblue1", span = 0.1) +
  geom_point(aes(x = WinCenter*2 , y = taj_rup),
             color = "skyblue1", alpha = 0.5, size = 1) +
  geom_rect(data=outlier_windows %>% filter(chrom == "CM078251.1",window_start==targ_win[1]), 
            inherit.aes=FALSE, 
            aes(xmin=window_start, xmax=window_end,
                ymin=-2.5,ymax=max(taj$taj_rup), 
                group=chrom), 
            color="transparent", 
            fill="indianred", 
            alpha=0.2)+
  geom_rect(xmin=35059641, xmax=35060636+2000, 
            ymin=-2.5,ymax=max(taj$taj_rup),
            color="transparent",
            fill="skyblue4",
            alpha=0.2)+
  geom_rect(xmin=35589506, xmax=35618151,
            ymin=-2.5,ymax=max(taj$taj_rup),
            color="transparent",
            fill="skyblue4",
            alpha=0.2)+
  labs(x = "Chromosome 36 (Mbp)",
       y = "Tajima's D") +
  scale_x_continuous(labels = function(x) x / 1000000) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        panel.spacing = unit(0,'lines'),
        axis.title.y =  element_text(angle = 90, vjust = 0.5),
        panel.grid = element_blank())


# taj_plot <- ggplot(taj %>% 
#                      filter(WinCenter > targ_win_plus[1] & WinCenter < targ_win_plus[2])) +
#   geom_smooth(aes(x = WinCenter , y = taj_kuu),
#              color = "#E69F00", span = 0.1) +
#   geom_point(aes(x = WinCenter , y = taj_kuu),
#               color = "#E69F00", alpha = 0.5, size = 1) +
#   geom_smooth(aes(x = WinCenter , y = taj_rup),
#              color = "skyblue1", span = 0.1) +
#   geom_point(aes(x = WinCenter , y = taj_rup),
#              color = "skyblue1", alpha = 0.5, size = 1) +
#   geom_rect(data=outlier_windows %>% filter(chrom == "CM078251.1",window_start==targ_win[1]), 
#             inherit.aes=FALSE, 
#             aes(xmin=window_start, xmax=window_end,
#                 ymin=-2.5,ymax=max(taj$taj_rup), 
#                 group=chrom), 
#             color="transparent", 
#             fill="indianred", 
#             alpha=0.2)+
#   labs(x = "Chromosome 36 (Mbp)",
#        y = "Tajima's D") +
#   scale_x_continuous(labels = function(x) x / 1000000) +
#   theme_bw(base_size = 16) +
#   theme(legend.position = "none",
#         panel.spacing = unit(0,'lines'),
#         axis.title.y =  element_text(angle = 90, vjust = 0.5),
#         panel.grid = element_blank())


# ggplot(taj %>% filter(WinCenter > targ_win_plus[1] & WinCenter < targ_win_plus[2])) +
#   geom_smooth(aes(x = WinCenter , y = taj_kuu),
#               color = "#E69F00", span = 0.1) +
#   geom_smooth(aes(x = WinCenter , y = taj_rup),
#               color = "skyblue1", span = 0.1) +
#   geom_rect(data=outlier_windows %>% filter(chrom == "CM078251.1",window_start==targ_win[1]), 
#             inherit.aes=FALSE, 
#             aes(xmin=window_start, xmax=window_end,
#                 ymin=-2.5,ymax=max(taj$taj_rup), 
#                 group=chrom), 
#             color="transparent", 
#             fill="indianred", 
#             alpha=0.2)+
#   geom_rect(xmin=35059641, xmax=35060636, ymin=-2.5,ymax=5, 
#             color="transparent", 
#             fill="blue", 
#             alpha=0.2)+
#   geom_rect(xmin=35589506, xmax=35618151, ymin=-2.5,ymax=5, 
#             color="transparent", 
#             fill="blue", 
#             alpha=0.2)+
#   labs(x = "",
#        y = "Tajima's") +
#   theme_bw(base_size = 16) +
#   theme(legend.position = "none",
#         panel.spacing = unit(0,'lines'),
#         axis.title.y =  element_text(angle = 90, vjust = 0.5),
#         strip.text.x = element_text(angle = 90))

#######################
#Analyse Tajima

mean(pi_window_nots$Tajima, na.rm = T)

median(pi_window_nots %>% filter(Chr == "CM078251.1",WinCenter > targ_win[1] & WinCenter < targ_win[2]) %>% pull(Tajima))


####################
#Analyse dxy

dxy %<>% mutate(mid = (start+end)/2)

dxy_plot <- ggplot(dxy) +
  geom_point(aes(x = mid , y = dxy),
             color = "black", alpha = 0.6, size = 1) +
  geom_smooth(aes(x = mid , y = dxy),
             color = "grey50",span = 0.5) +
  geom_rect(data=outlier_windows %>% filter(chrom == "CM078251.1",window_start==targ_win[1]), 
            inherit.aes=FALSE, 
            aes(xmin=window_start, xmax=window_end,
                ymin=min(dxy$dxy),ymax=max(dxy$dxy), 
                group=chrom), 
            color="transparent", 
            fill="indianred", 
            alpha=0.2)+
    geom_rect(xmin=35059641, xmax=35060636+2000, ymin=min(dxy$dxy),ymax=max(dxy$dxy),
              color="transparent",
              fill="skyblue4",
              alpha=0.2)+
    geom_rect(xmin=35589506, xmax=35618151, ymin=min(dxy$dxy), ymax=max(dxy$dxy),
              color="transparent",
              fill="skyblue4",
              alpha=0.2)+
  labs(x = "",
       y = "dxy") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        panel.spacing = unit(0,'lines'),
        axis.title.y =  element_text(angle = 90, vjust = 0.5),
        strip.text.x = element_text(angle = 90),
        axis.text.x = element_blank(),
        panel.grid = element_blank())

################################################3
#Depth

depth <- fread("data/chr34/id_merged_coverage.txt.gz")

colnames(depth) <- c("id","junk","chrom","start","end","coverage")

metadata <- read.table("data/chr34/karyo_info.tsv",header=T)

depth %<>% left_join(metadata, by = "id")

head(depth)

# Define window size
window_size <- 100


depth_avg <- depth %>%
  filter(!is.na(geno_clus)) %>%
  mutate(window_start = floor(start / window_size) * window_size,
         window_end = window_start + window_size) %>%
  group_by(geno_clus, chrom, window_start, window_end) %>%
  dplyr::summarize(mean_cov = mean(coverage, na.rm = T),
            .groups = "drop")

# Rename columns for clarity
colnames(depth_avg) <- c("geno_clus", "id", "window_start", "window_end", "mean_coverage")


targ_win_plus2 <- c(34800000-100000,35619999+100000)

# Plot
depth_plot <- ggplot(depth_avg %>% filter(!is.na(geno_clus), mean_coverage < 50,
                                          window_start > targ_win_plus2[1] & window_start < targ_win_plus2[2])) +
  geom_line(aes(x = window_start, y = mean_coverage, color = geno_clus)) +
  labs(x = "",
       y = "Coverage") +
  theme_bw(base_size = 16) +
  geom_rect(data = outlier_windows %>% filter(chrom == "CM078251.1", 
                                              window_start == targ_win[1]), 
            inherit.aes = FALSE, 
            aes(xmin = window_start, xmax = window_end,
                ymin = 0, ymax = 10), 
            color = "transparent", 
            fill = "indianred", 
            alpha = 0.2) +
  geom_rect(xmin=35059641, xmax=35060636, ymin=0,ymax=10,
            color="transparent",
            fill="blue")+
  scale_color_manual(values = c("#E69F00", "slateblue4", "steelblue"))+
  facet_wrap(~ geno_clus, ncol = 1) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.spacing = unit(0, 'lines'),
        axis.title.y = element_text(angle = 90, vjust = 0.5),
        axis.text.x = element_blank())

ggsave("figures/05_depth_outlier window.png", depth_plot,
       width = 20, height = 10, unit = "cm", dpi = 900)


depth_plot <- ggplot(depth_avg %>% filter(!is.na(geno_clus), mean_coverage < 50,
                                          window_start > 34750000 & window_start < 34850000)) +
  geom_line(aes(x = window_start, y = mean_coverage, color = geno_clus)) +
  labs(x = "",
       y = "Coverage") +
  theme_bw(base_size = 16) +
  geom_rect(data = outlier_windows %>% filter(chrom == "CM078251.1", 
                                              window_start == targ_win[1]), 
            inherit.aes = FALSE, 
            aes(xmin = window_start, xmax = 34850000,
                ymin = 0, ymax = 10), 
            color = "transparent", 
            fill = "indianred", 
            alpha = 0.2) +
  scale_color_manual(values = c("#E69F00", "slateblue4", "steelblue"))+
  facet_wrap(~ geno_clus, ncol = 1) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.spacing = unit(0, 'lines'),
        axis.title.y = element_text(angle = 90, vjust = 0.5),
        axis.text.x = element_blank())


ggsave("figures/05_depth_outlier window_zoom.png", depth_plot,
       width = 20, height = 10, unit = "cm", dpi = 900)


#################################################
#final figure



outlier_chr34 <- fst_plot/ dxy_plot/ taj_plot & 
  theme(plot.margin = margin(t = 1, r = 5, b = 0, l = 5, unit = "pt"))

outlier_chr34

ggsave("figures/05_chr34_outlier.png", outlier_chr34,height =16,width = 12, unit = "cm", dpi = 900)

















