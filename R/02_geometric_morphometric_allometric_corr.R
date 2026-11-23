rm(list = ls())

###############################################################################
# 02b - Geometric morphometrics: allometry (Fig. 1d, Supplementary Tables)
#
# Needs: tmp/final_data_gpa_projected.Rds and tmp/real_csize.Rds
#        (made by 02_geometric_morphometric_arching_biais_corr.R)
#
# CHANGES vs the previous version
#  1. Size = the fish's real centroid size (tmp/real_csize.Rds), not the Csize
#     of the re-aligned shapes, which was ~1 for every fish. ALL size results
#     change: the size effect, the size x ecotype interaction, Fig. 1d, and the
#     ecotype/sex R-squared values (they are adjusted for size).
#  2. Metadata are matched to shapes by fish ID, not by row order.
#  3. Every test uses 9,999 iterations (= 10,000 permutations including the
#     observed data, as in the Methods). The group-effects table used 1,000.
#  4. Removed `logsz = T` from procD.lm(): it isn't a procD.lm() argument
#     and was silently ignored.
#  5. Removed the duplicated code at the end of the old script.
###############################################################################

#################
#### Library ####
#################

library(geomorph)
library(RRPP)
library(tidyverse)
library(magrittr)
library(openxlsx)

source("99_morpho_helpers.R", local = TRUE)
dir.create("figures", showWarnings = FALSE)

##################################
### Load data                  ###
##################################

final_data_gpa_projected <- readRDS("tmp/final_data_gpa_projected.Rds")
real_csize <- readRDS("tmp/real_csize.Rds")

spec_ids <- dimnames(final_data_gpa_projected$coords)[[3]]
stopifnot(identical(names(real_csize), spec_ids))   # sizes are in the same fish order as the shapes

data <- read.csv("data/data_with_gillrakers.csv", header = TRUE, sep = ",")
data %<>% mutate(river = if_else(str_detect(sample_id, 'NOT'), "NOT", "RUP"))

### CHANGED: match metadata to shapes by ID
data <- align_to_shapes(data, spec_ids)

gdf <- geomorph.data.frame(shape = final_data_gpa_projected$coords,
                           pop   = factor(data$river),
                           sex   = factor(data$sex),
                           eco   = factor(data$true_eco),
                           cs    = as.numeric(real_csize))     ### CHANGED: real size

##################################
### Is allometry significant?  ###
##################################

ciscoAllometry <- procD.lm(shape ~ log(cs), data = gdf, iter = 9999,
                           SS.type = "III", print.progress = FALSE)
summary(ciscoAllometry)

##################################################
### Do allometric slopes differ between ecotypes?
##################################################

common_allo <- procD.lm(shape ~ log(cs) + eco + sex + pop,
                        data = gdf, iter = 9999, SS.type = "II", print.progress = FALSE)
diff_allo   <- procD.lm(shape ~ log(cs) + eco + sex + pop + log(cs):eco,
                        data = gdf, iter = 9999, SS.type = "II", print.progress = FALSE)

onto <- anova(common_allo, diff_allo, print.progress = FALSE)
onto

anova_onto_df <- as.data.frame(onto$table) %>%
  rownames_to_column("Term")

write_sheet(anova_onto_df, "Ontogenic_allometric_slope_diff")

# Pairwise comparison of slopes (length and angle). Not reported in the
# manuscript, but useful if a reviewer asks how the slopes differ.
PW.u <- pairwise(diff_allo, groups = gdf$eco,
                 covariate = log(gdf$cs), print.progress = FALSE)
summary(PW.u, test = "dist")
summary(PW.u, test = "VC", angle.type = "deg")

#####################################################################
### Effects of size, ecotype, sex and river on shape (full model) ###
#####################################################################

# Same model as diff_allo. The R-squared column gives the share of shape
# variance per term (the old 1.05% ecotype / 3.49% size / 4.38% sex values).
anova_tab <- anova(diff_allo)
anova_tab

anova_df <- as.data.frame(anova_tab$table) %>%
  rownames_to_column("Effect")

write_sheet(anova_df, "Allometry_group_effects")

##############################
### Fig. 1d               ###
##############################

data %<>% mutate(color_allo = case_when(true_eco == "summer_run" ~ "#E69F00",
                                        true_eco == "fall_run"   ~ "skyblue1"),
                 shape_allo = case_when(sex == "F" ~ 1,
                                        sex == "M" ~ 19,
                                        TRUE ~ 4))    # any other sex code (e.g. R) = x

png("figures/02_allometry_fit.png",
    width = 1200, height = 1200, res = 300, bg = "transparent")

par(cex.lab = 1.2, cex.axis = 1, mar = c(4.5, 5.5, 2, 2),
    col.axis = "grey25", col.lab = "black", col = "black",
    bty = "n", mgp = c(3.3, 0.3, 0))   # axis titles 3.3 lines out, clear of the tick labels

# plotAllometry() logs the size you give it (logsz = TRUE by default), so the
# x-axis is log(centroid size).
plotAllometry(fit    = diff_allo,
              size   = gdf$cs,
              method = "PredLine",
              pch    = data$shape_allo,
              col    = data$color_allo,
              xlab   = "log(centroid size)",
              ylab   = "",
              xaxt   = "n",
              yaxt   = "n")

# signif() keeps 2-3 significant digits, so small values like -0.005 don't
# all round to "0.00" (and the minus sign is kept).
yticks <- axTicks(2)
axis(2, at = yticks, labels = signif(yticks, 2),
     las = 1, col.axis = "grey25", lwd = 0.5, tck = -0.015)

xticks <- axTicks(1)
axis(1, at = xticks, labels = signif(xticks, 3),
     las = 1, col.axis = "grey25", lwd = 0.5, tck = -0.015)

usr <- par("usr")
rect(usr[1], usr[3], usr[2], usr[4], border = "grey25", lwd = 1.5, xpd = FALSE)

dev.off()
