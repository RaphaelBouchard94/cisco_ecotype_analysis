rm(list = ls())

###############################################################################
# 02a - Geometric morphometrics: remove body arching (Burnaby projection)
#
# Run this first: it creates the files in tmp/ used by the other 02 scripts.
#
# CHANGES vs the previous version
#  1. Saves each fish's REAL centroid size (from the first Procrustes fit).
#     The old scripts took centroid size from gpagen() run on coordinates that
#     were already scaled to size 1, so every fish had a "size" of ~1.0001-1.0007.
#     That is why the Fig. 1d x-axis ran from 0.0001 to 0.0005 (it is log size).
#  2. Uses geomorph's two.d.array() / arrayspecs() for the reshaping. They do
#     what the old loops did, but keep the specimen IDs attached to the shapes.
#  3. Tests the angle between uncorrected PC1 and the first FOUR corrected PCs
#     (as the Methods say), not just PC1 and PC2, and saves the table.
###############################################################################

#################
#### Library ####
#################

library(geomorph)
library(tidyverse)
library(magrittr)

source("99_ProjectOrthogonal.R", local = TRUE)
source("99_TestOfAngle.R", local = TRUE)
source("99_critical_angle.R", local = TRUE)
source("99_morpho_helpers.R", local = TRUE)

rad2deg <- function(rad) {(rad * 180) / (pi)}   # used inside TestOfAngle()

dir.create("tmp", showWarnings = FALSE)
dir.create("res", showWarnings = FALSE)

#################
### Load Data ###
#################

data_tps <- readland.tps("data/cisco_rup_bay.TPS", specID = "ID")

n_landmarks <- dim(data_tps)[1]   # 16 in your file
n_fish      <- dim(data_tps)[3]   # 190 in your file
cat("Landmarks:", n_landmarks, " Fish:", n_fish, "\n")

data_gpa <- gpagen(data_tps, print.progress = FALSE)

spec_ids <- dimnames(data_gpa$coords)[[3]]

### CHANGED: keep the real centroid size of every fish ########################
# gpagen() reports the centroid size of the coordinates you give it, BEFORE it
# scales them to unit size. Here those are the raw digitised coordinates, so
# this is the fish's actual size (in mm if your TPS file has SCALE= lines,
# otherwise in pixels; either way it is fine on a log scale).
real_csize <- setNames(data_gpa$Csize, spec_ids)
write_rds(real_csize, "tmp/real_csize.Rds")
summary(real_csize)

PCA <- gm.prcomp(data_gpa$coords)

#############################
## Remove body arching     ##
#############################

### CHANGED: simpler reshaping that keeps the IDs ##############################
# two.d.array() turns the 16 x 2 x n array into an n x 32 matrix with columns
# x1, y1, x2, y2, ... which is the same order as the rows of PCA$rotation.
X <- two.d.array(data_gpa$coords)

# Burnaby's orthogonal projection: remove all variation along uncorrected PC1,
# which we treat as the body-arching axis (Valentin et al. 2008).
X_proj <- ProjectOrthogonal(X, PCA$rotation[, 1])

# arrayspecs() turns the n x 32 matrix back into a 16 x 2 x n array.
data_gpa_proj <- arrayspecs(X_proj, p = n_landmarks, k = 2)
dimnames(data_gpa_proj)[[3]] <- spec_ids

write_rds(data_gpa_proj, "tmp/data_gpa_proj.Rds")

# Re-align the corrected configurations (Methods: "re-aligned using GPA").
# NB: the Csize stored in this object is NOT body size (the input is already
# unit size). Always use tmp/real_csize.Rds for size.
final_data_gpa_projected <- gpagen(data_gpa_proj, print.progress = FALSE)
write_rds(final_data_gpa_projected, "tmp/final_data_gpa_projected.Rds")

######################################################
## Check that arching was removed (angle tests)     ##
######################################################

# If the projection worked, every corrected PC should be ~90 degrees from the
# uncorrected PC1 (the arching axis).
pca_final <- gm.prcomp(final_data_gpa_projected$coords)

angle_tab <- map_dfr(1:4, function(i) {
  a <- TestOfAngle(PCA$rotation[, 1], pca_final$rotation[, i], flip = FALSE)
  tibble(comparison = paste0("uncorrected PC1 vs corrected PC", i),
         angle_deg  = a$angle_deg,
         p_value    = a$p.value)
})

print(angle_tab)

write_sheet(angle_tab, "Arching_correction_angles")
