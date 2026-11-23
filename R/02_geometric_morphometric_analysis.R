rm(list = ls())

###############################################################################
# 02c - Geometric morphometrics: shape PCA (Fig. 1e) and shape-based LDA
#
# Needs: tmp/final_data_gpa_projected.Rds (from the arching script)
#
# CHANGES vs the previous version
#  1. Metadata are matched to shapes by fish ID, not by row order.
#  2. The PC1-PC4 wireframes now come from the SAME PCA as the scatter plots.
#     (They used to come from a different PCA, run before the re-alignment,
#     so a PC's min/max could be flipped relative to the plotted axis.)
#  3. Axis labels (% variance) are computed, not typed in by hand.
#  4. The number of PCs for the LDA comes from the broken-stick rule, as the
#     Methods say (the old script only drew a screeplot and used 4).
#  5. LDA: equal priors (the default prior favours the bigger group), balanced
#     accuracy, and a permutation test against chance.
#  6. Dropped unused packages. plyr in particular masks dplyr functions.
###############################################################################

#################
#### Library ####
#################

library(MASS)        # load before tidyverse so dplyr::select() isn't masked
library(geomorph)
library(tidyverse)
library(magrittr)
library(patchwork)
library(openxlsx)

source("99_morpho_helpers.R", local = TRUE)
dir.create("figures", showWarnings = FALSE)

N_PERM <- 1000       # permutations for the LDA test
N_PC   <- NA         # PCs used in the LDA: NA = broken-stick rule; or a number, e.g. 4

#################
### Load data ###
#################

final_data_gpa_projected <- readRDS("tmp/final_data_gpa_projected.Rds")
spec_ids <- dimnames(final_data_gpa_projected$coords)[[3]]

data <- read.csv("data/data_with_gillrakers.csv", header = TRUE, sep = ",")
data %<>% mutate(river = if_else(str_detect(sample_id, 'NOT'), "NOT", "RUP"))
data <- align_to_shapes(data, spec_ids)       ### CHANGED

##############################################
###### PCA of the arching-corrected shapes ###
##############################################

pca_gm <- gm.prcomp(final_data_gpa_projected$coords)
summary(pca_gm)

### CHANGED: % variance computed from the PCA itself
pct <- round(100 * pca_gm$sdev^2 / sum(pca_gm$sdev^2), 2)
lab <- function(i) paste0("PC", i, " (", pct[i], "%)")

pca_ggplot <- data.frame(sex = data$sex,
                         eco = data$true_eco,
                         PC1 = pca_gm$x[, 1],
                         PC2 = pca_gm$x[, 2],
                         PC3 = pca_gm$x[, 3],
                         PC4 = pca_gm$x[, 4])

pc_plot <- function(xpc, ypc) {
  pca_ggplot %>%
    filter(sex %in% c("F", "M")) %>%
    ggplot(aes(x = .data[[xpc]], y = .data[[ypc]], color = eco, shape = sex)) +
    geom_point(size = 2) +
    stat_ellipse(aes(linetype = sex)) +
    xlab(lab(as.integer(str_remove(xpc, "PC")))) +
    ylab(lab(as.integer(str_remove(ypc, "PC")))) +
    scale_color_manual(name = "Ecotype",
                       breaks = c("summer_run", "fall_run"),
                       labels = c("Kâcikâsikumekw", "Nûtamesânîw"),
                       values = c("summer_run" = "#E69F00", "fall_run" = "skyblue1")) +
    scale_shape_manual(name = "Sex", labels = c("Female", "Male"), values = c(1, 16)) +
    scale_linetype_manual(name = "Sex", labels = c("Female", "Male"), values = c("solid", "dotted")) +
    theme_bw(base_size = 14) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          legend.position = "none")
}

geo_pc12 <- pc_plot("PC1", "PC2")
geo_pc34 <- pc_plot("PC3", "PC4")

ggsave("figures/02_geomorph_pc12.png", geo_pc12, unit = "cm", width = 10, height = 10, dpi = 600)
ggsave("figures/02_geomorph_pc34.png", geo_pc34, unit = "cm", width = 10, height = 10, dpi = 600)
ggsave("figures/02_geomorph_pc1234.png", geo_pc12 + geo_pc34, unit = "cm", width = 20, height = 10, dpi = 600)

#############################################################
### Shape change along each PC (wireframes in Fig. 1e)    ###
#############################################################

# Links between landmarks used to draw the outline (same as before).
links <- cbind(c(2, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16),
               c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 2))

### CHANGED: shapes at the PC extremes come from pca_gm (the plotted PCA).
# mag = 2 exaggerates the shape change two-fold: say so in the caption.
for (i in 1:4) {
  png(paste0("figures/02_PC", i, "_morph_change.png"),
      width = 3000, height = 3000, res = 300, bg = "transparent")
  plotRefToTarget(pca_gm$shapes[[paste0("shapes.comp", i)]]$min,
                  pca_gm$shapes[[paste0("shapes.comp", i)]]$max,
                  method = "points", mag = 2, links = links,
                  gridPars = gridPar(pt.bg = "grey40", tar.pt.bg = "black",
                                     pt.size = 3, link.lwd = 3, link.col = "grey40",
                                     tar.link.col = "black", tar.link.lwd = 3,
                                     tar.pt.size = 3))
  dev.off()
}

#############################################################
### Cross-validated LDA on PC scores                      ###
#############################################################

### CHANGED: number of PCs from the broken-stick rule
bs <- broken_stick(pca_gm$sdev^2)
print(head(bs$table, 8))
cat("Broken-stick keeps", bs$n, "PCs\n")
n_pc <- if (is.na(N_PC)) bs$n else N_PC
if (n_pc < 1) stop("Broken-stick kept no PCs; set N_PC at the top of the script.")
if (!is.na(N_PC) && N_PC != bs$n)
  message("Using ", N_PC, " PCs as set by N_PC (broken-stick would keep ", bs$n,
          "). Don't cite the broken-stick rule in the Methods in that case.")

pcs <- pca_gm$x[, 1:n_pc, drop = FALSE]

lda_shape <- lda_loo(pcs, data$true_eco)            # leave-one-out, equal priors
lda_shape$confusion
lda_shape$sensitivity                               # success rate per ecotype
cat("Balanced accuracy:", round(lda_shape$balanced_accuracy, 3), "\n")

p_shape <- lda_perm_test(pcs, data$true_eco, lda_shape$balanced_accuracy, n_perm = N_PERM)
cat("Permutation p-value (vs chance = 0.5):", signif(p_shape, 3), "\n")

write_sheet(data.frame(n_PCs = n_pc,
                       n_fish = nrow(pcs),
                       balanced_accuracy = lda_shape$balanced_accuracy,
                       success_fall_run   = lda_shape$sensitivity[["fall_run"]],
                       success_summer_run = lda_shape$sensitivity[["summer_run"]],
                       permutation_p = p_shape,
                       n_permutations = N_PERM),
            "LDA_shape")
