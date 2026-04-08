rm(list=ls())

#################
#### Library ####
#################

library(tidyverse)
library(magrittr)
library(geomorph)
library(abind)
library(Morpho)
library(MASS)
library(ggpubr)
library(plyr)
library(cowplot)
library(ggExtra)
library(png)
library(patchwork)
library(car)
library(multcomp)
library(RRPP)
library(caret)
library(klaR)
library(openxlsx)

##################################

#Allometric relationship analysis

final_data_gpa_projected <- readRDS("tmp/final_data_gpa_projected.Rds")

data <- read.csv("data/data_with_gillrakers.csv", header = T,sep = ",")

data %<>% mutate(river = if_else(str_detect(sample_id,'NOT'),"NOT","RUP"))

gdf <- geomorph.data.frame(shape = final_data_gpa_projected$coords,
                           pop = data$river,
                           sex = data$sex,
                           eco = data$true_eco,
                           cs = final_data_gpa_projected$Csize)

#Is allometry significant?

ciscoAllometry <- procD.lm(shape ~ cs, data = gdf, iter = 10000,SS.type = "III")

summary(ciscoAllometry)

plot(ciscoAllometry)

#Yes! Effect of size on shape

#Test allometric slope heterogeneity

common_allo <- procD.lm(shape ~ log(cs) + eco + sex + pop, logsz = T, data = gdf, iter = 9999,SS.type = "II")
diff_allo <- procD.lm(shape ~ log(cs) + eco + sex + log(cs):eco + pop, logsz = T, data = gdf, iter = 9999,SS.type = "II")

anova(common_allo, diff_allo)

#############################################################
#Save results of ontogenic difference in supplementary table
onto <- anova(common_allo, diff_allo)

anova_onto_df <- as.data.frame(onto$table)

# Add the model/term names as a column
anova_onto_df$Term <- rownames(anova_onto_df)

# Reorder columns for readability
anova_onto_df <- anova_onto_df %>%
  dplyr::select(Term, ResDf, Df, RSS, SS, MS, Rsq, "F", Z, P, `Pr(>F)`)

# Create workbook and add a single sheet
wb <- loadWorkbook("res/supplementary_table.xlsx")
addWorksheet(wb, "Ontogenic_allometric_slope_diff")
writeData(wb, sheet = "Ontogenic_allometric_slope_diff", anova_onto_df)

# Save Excel file
saveWorkbook(wb, "res/supplementary_table.xlsx", overwrite = T)

#Hypothesis tests; vector lengths
PW.u <- pairwise(diff_allo, groups = gdf$eco,
                 covariate = log(gdf$cs), print.progress = FALSE)

summary(PW.u, test = "dist")

#Hypothesis tests: angles between vectors
summary(PW.u, test = "VC", angle.type = "deg")

#Suggests allometric slopes are not parallel and significantly diverge

####################################################################
#Is there difference between group on this allometric relationship?

ciscoAllometry_group_effect <- procD.lm(shape ~ log(cs) + eco + sex + log(cs):eco + pop, logsz = T, data = gdf, iter = 1000,SS.type = "II")

summary(ciscoAllometry_group_effect)

anova_tab <- anova(ciscoAllometry_group_effect)

anova_df <- as.data.frame(anova_tab$table)

# Keep term names explicit
anova_df$Effect <- rownames(anova_df)
rownames(anova_df) <- NULL

# Optional: reorder columns nicely
anova_df <- anova_df %>%
  dplyr::select(
    Effect, Df, SS, MS, Rsq, "F", Z, "Pr(>F)"
  )

###TO ADJUST

data %<>% mutate(color_allo = case_when(true_eco == "summer_run" ~ "#E69F00",
                                        true_eco == "fall_run" ~ "skyblue1"),
                 shape_allo = case_when(sex == "F" ~ 1,
                                        sex == "M" ~ 19))

##############################
plotAllometry(f = ciscoAllometry_group_effect, 
              size = gdf$cs, 
              method = "RegScore",
              pch = data$shape_allo,
              col = data$color_allo,
              xlab = "Centroid Size",
              ylab = "Shape (Predicted)")

##########
png("figures/02_allometry_fit.png",
    width = 1200,
    height = 1200,
    res = 300,
    bg = "transparent")

# ---- Graphics parameters ----
par(
  cex.lab  = 1.2,
  cex.axis = 1,
  mar      = c(4.5, 5.5, 2, 2),
  col.axis = "grey25",
  col.lab  = "black",
  col      = "black",
  bty      = "n",
  mgp      = c(2, 0.3, 0)   # prevent automatic box
)

# ---- Main plot (suppress default axes) ----
plotAllometry(
  fit    = ciscoAllometry_group_effect,
  size   = gdf$cs,
  method = "PredLine",
  pch    = data$shape_allo,
  col    = data$color_allo,
  xlab   = "Centroid Size",
  ylab   = "",
  xaxt   = "n",
  yaxt   = "n"
)

# ---- Custom Y axis ----
yticks <- axTicks(2)
yticks_sub <- yticks[seq(1, length(yticks), by = 2)]

axis(2,
     at     = yticks_sub,
     labels = sprintf("%.2f", yticks_sub),
     las    = 1,
     col.axis = "grey25",
     lwd = 0.5,
     tck = -0.015)

# ---- Custom X axis ----
xticks <- axTicks(1)
xticks_sub <- xticks[seq(1, length(xticks), by = 2)]

axis(1,
     at     = xticks_sub,
     labels = sprintf("%.4f", xticks_sub),
     las    = 1,
     col.axis = "grey25",
     lwd = 0.5,
     tck = -0.015)

# ---- Thin full panel border ----
usr <- par("usr")
rect(usr[1], usr[3], usr[2], usr[4],
     border = "grey25",
     lwd = 1.5,
     xpd = FALSE)

dev.off()







# Add the model/term names as a column
anova_onto_df$Term <- rownames(anova_onto_df)

# Reorder columns for readability
anova_onto_df <- anova_onto_df %>%
  dplyr::select(Term, ResDf, Df, RSS, SS, MS, Rsq, "F", Z, P, `Pr(>F)`)

# Create workbook and add a single sheet

wb <- loadWorkbook("res/supplementary_table.xlsx")

addWorksheet(wb, "Allometry_group_effects")

writeData(
  wb,
  sheet = "Allometry_group_effects",
  anova_df
)

saveWorkbook(wb, "res/supplementary_table.xlsx", overwrite = TRUE)
