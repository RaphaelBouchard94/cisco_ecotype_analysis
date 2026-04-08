rm(list = ls())

#################
#### Library ####
#################

library(tidyverse)
library(magrittr)
library(MASS)
library(vegan)
library(effectsize)
library(reshape2)
library(ggridges)
library(patchwork)
library(forcats)
library(ggrepel)
library(openxlsx)

#################
### Load Data ###
#################

data <- read.csv("data/data_gillrak_geomorph_age.csv", header = TRUE)

data %<>%
  mutate(gill_raker_count = upper_gill_raker + lower_gill_raker) %>%
  filter(sex %in% c("F", "M"))

data$true_eco <- as.factor(data$true_eco)
data$sex <- as.factor(data$sex)

data %<>% mutate(river = if_else(str_detect(sample_id,'NOT'),"NOT","RUP"))

################################################
# Correlation with fork length (diagnostic only)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "complete.obs"))
  txt <- formatC(r, format = "f", digits = digits)
  if (missing(cex.cor)) cex.cor <- 0.8
  text(0.5, 0.5, txt, cex = cex.cor)
}

pairs(
  data[, c(4, 12:30)],
  upper.panel = panel.cor,
  lower.panel = panel.smooth
)

###################################################
# Length standardisation (residuals)

traits_to_correct <- 13:29

for (i in traits_to_correct) {
  m <- lm(log(data[[i]]) ~ log(data$long_fourche))
  data[[paste0(names(data)[i], "_res")]] <- rstandard(m)
}

#######################################
# MANOVA

data_manova <- data %>%
  dplyr::select(
    river, sex, true_eco,
    AMO,
    MXL_res:HDD_res,
    -AMO_res,
    middle_gill_raker_length,
    gill_raker_count
  ) %>%
  drop_na()

manova_fit <- manova(
  cbind(
    AMO, MXL_res, MDL_res, POL_res, OOL_res, PSL_res,
    TTL_res, DOL_res, LUL_res, ANL_res, CPL_res, CPD_res,
    PCL_res, PVL_res, BDD_res, HDD_res,
    gill_raker_count, middle_gill_raker_length
  ) ~ true_eco * sex + river,
  data = data_manova
)

summary(manova_fit, test = "Wilks")


aov_list <- summary.aov(manova_fit)

# Combine into a single data frame
manova_combined <- bind_rows(
  lapply(names(aov_list), function(trait) {
    df <- as.data.frame(aov_list[[trait]])
    df$Trait <- trait             # add column for trait name
    df$Term <- rownames(df)       # rownames = model terms
    df <- df[, c("Trait", "Term", "Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")]
    return(df)
  })
)

# Create workbook and add a single sheet
# wb <- createWorkbook()
# addWorksheet(wb, "MANOVA_linear_traits")
# writeData(wb, sheet = "MANOVA_linear_traits", manova_combined)
# 
# # Save Excel file
# saveWorkbook(wb, "res/supplementary_table.xlsx", overwrite = TRUE)
# 

sign_var_eco <- c("DOL_res","ANL_res","CPL_res","gill_raker_count","middle_gill_raker_length")

#Plot individual trait

prop_var_grk  <- aov_list$` Response gill_raker_count`$`Sum Sq`[1] / (aov_list$` Response gill_raker_count`$`Sum Sq`[1] + aov_list$` Response gill_raker_count`$`Sum Sq`[5]) * 100

f_grk <- aov_list$` Response gill_raker_count`$`F value`[1]

label_text_grk <- paste0("F = ", round(f_grk,3) , "\nPVE = ", round(prop_var_grk,3))

grk <- ggplot(data %>% mutate(true_eco = fct_relevel(true_eco, "summer_run","fall_run")), aes(x = true_eco, y = gill_raker_count)) +
  geom_boxplot(fill = "white", width = 0.3, notch = T,outlier.alpha = 0) +
  geom_jitter(aes(color = true_eco), width = 0.1, alpha = 0.6)+
  scale_color_manual(values = c("#E69F00","skyblue1"))+
  scale_x_discrete(labels = c("KUU","NUT"))+
  xlab("")+
  ylab("Total rakers counts")+
  theme_bw(base_size = 14)+
  theme(panel.border = element_rect(colour = "black",fill=NA),
        panel.grid = element_blank(),
        legend.position = "none")+
  annotate("text", 
           x = 0.5, # position on x-axis (adjust if needed)
           y = 32, 
           label = label_text_grk, 
           hjust = 0, vjust = -0.2, size = 4)

grk

ggsave("figures/01_totgrc.png",width = 10, height = 10, unit = "cm", dpi = 600)

prop_var_mgr  <- aov_list$` Response middle_gill_raker_length`$`Sum Sq`[1] / (aov_list$` Response middle_gill_raker_length`$`Sum Sq`[1] + aov_list$` Response middle_gill_raker_length`$`Sum Sq`[5]) * 100

f_mgr <- aov_list$` Response middle_gill_raker_length`$`F value`[1]

label_text_mgr <- paste0("F = ", round(f_mgr,3) , "\nPVE = ", round(prop_var_mgr,3))

mgrk <- ggplot(data %>% mutate(true_eco = fct_relevel(true_eco, "summer_run","fall_run")), aes(x = true_eco, y = middle_gill_raker_length)) +
  geom_boxplot(fill = "white", width = 0.3, notch = T,outlier.alpha = 0) +
  geom_jitter(aes(color = true_eco), width = 0.1, alpha = 0.6)+
  scale_color_manual(values = c("#E69F00","skyblue1"))+
  scale_x_discrete(labels = c("KUU","NUT"))+
  xlab("")+
  ylab("Middle gill raker length")+
  theme_bw(base_size = 14)+
  theme(panel.border = element_rect(colour = "black",fill=NA),
        panel.grid = element_blank(),
        legend.position = "none")+
  annotate("text", 
           x = 0.5, 
           y = 33, 
           label = label_text_mgr, 
           hjust = 0, vjust = -0.2, size = 4)

mgrk

ggsave("figures/01_mgrc.png",width = 10, height = 10, unit = "cm", dpi = 600)

prop_var_anl  <- aov_list$` Response ANL_res`$`Sum Sq`[1] / (aov_list$` Response ANL_res`$`Sum Sq`[1] + aov_list$` Response ANL_res`$`Sum Sq`[5]) * 100

f_anl <- aov_list$` Response ANL_res`$`F value`[1]

label_text_anl <- paste0("F = ", round(f_anl,3) , "\nPVE = ", round(prop_var_anl,3))

anl <- ggplot(data %>% mutate(true_eco = fct_relevel(true_eco, "summer_run","fall_run")), aes(x = true_eco, y = ANL/long_fourche)) +
  geom_boxplot(fill = "white", width = 0.3, notch = T, outlier.alpha = 0) +
  geom_jitter(aes(color = true_eco), width = 0.1, alpha = 0.6)+
  scale_color_manual(values = c("#E69F00","skyblue1"))+
  scale_x_discrete(labels = c("KUU","NUT"))+
  xlab("")+
  ylab("Anal")+
  theme_bw(base_size = 14)+
  theme(panel.border = element_rect(colour = "black",fill=NA),
        panel.grid = element_blank(),
        legend.position = "none")+
  annotate("text", 
           x = 0.5, 
           y = 0.007, 
           label = label_text_anl, 
           hjust = 0, vjust = -0.2, size = 4)

anl

ggsave("figures/01_anl.png",width = 10, height = 10, unit = "cm", dpi = 600)

prop_var_dol  <- aov_list$` Response DOL_res`$`Sum Sq`[1] / (aov_list$` Response DOL_res`$`Sum Sq`[1] + aov_list$` Response DOL_res`$`Sum Sq`[5]) * 100

f_dol <- aov_list$` Response DOL_res`$`F value`[1]

label_text_dol <- paste0("F = ", round(f_dol,3) , "\nPVE = ", round(prop_var_dol,3))

dol <- ggplot(data %>% mutate(true_eco = fct_relevel(true_eco, "summer_run","fall_run")), aes(x = true_eco, y = DOL/long_fourche)) +
  geom_boxplot(fill = "white", width = 0.3, notch = T, outlier.alpha = 0) +
  geom_jitter(aes(color = true_eco), width = 0.1, alpha = 0.6)+
  scale_color_manual(values = c("#E69F00","skyblue1"))+
  scale_x_discrete(labels = c("KUU","NUT"))+
  xlab("")+
  ylab("Dorsal")+
  theme_bw(base_size = 14)+
  theme(panel.border = element_rect(colour = "black",fill=NA),
        panel.grid = element_blank(),
        legend.position = "none")+
  annotate("text", 
           x = 0.5, 
           y = 0.008, 
           label = label_text_dol, 
           hjust = 0, vjust = -0.2, size = 4)

dol

ggsave("figures/01_dol.png",width = 10, height = 10, unit = "cm", dpi = 600)

prop_var_cpl  <- aov_list$` Response CPL_res`$`Sum Sq`[1] / (aov_list$` Response CPL_res`$`Sum Sq`[1] + aov_list$` Response CPL_res`$`Sum Sq`[5]) * 100

f_cpl <- aov_list$` Response CPL_res`$`F value`[1]

label_text_cpl <- paste0("F = ", round(f_cpl,3) , "\nPVE = ", round(prop_var_cpl,3))

cpl <- ggplot(data %>% mutate(true_eco = fct_relevel(true_eco, "summer_run","fall_run")), aes(x = true_eco, y = CPL/long_fourche)) +
  geom_boxplot(fill = "white", width = 0.3, notch = T, outlier.alpha = 0) +
  geom_jitter(aes(color = true_eco), width = 0.1, alpha = 0.6)+
  scale_color_manual(values = c("#E69F00","skyblue1"))+
  scale_x_discrete(labels = c("KUU","NUT"))+
  xlab("")+
  ylab("Caudal Peduncle")+
  theme_bw(base_size = 14)+
  theme(panel.border = element_rect(colour = "black",fill=NA),
        panel.grid = element_blank(),
        legend.position = "none")+
  annotate("text", 
           x = 0.5, 
           y = 0.0055, 
           label = label_text_cpl, 
           hjust = 0, vjust = -0.2, size = 4)

cpl

ggsave("figures/01_cpl.png",width = 10, height = 10, unit = "cm", dpi = 600)


#######################################
# LDA 

lda_data <- data_manova %>%
  dplyr::select(
    true_eco,
    ANL_res,
    DOL_res,
    CPL_res,
    gill_raker_count,
    middle_gill_raker_length
  ) %>%
  drop_na()

lda_data$true_eco <- droplevels(lda_data$true_eco)

set.seed(123)  # for reproducibility

train_prop <- 0.7  # 70% training, 30% testing

train_index <- lda_data %>%
  group_by(true_eco) %>%
  group_map(~ sample(seq_len(nrow(.)), size = floor(train_prop * nrow(.)))) %>%
  unlist()

train_data <- lda_data[train_index, ]
test_data  <- lda_data[-train_index, ]

lda_fit <- lda(true_eco ~ ., data = train_data)

lda_pred <- predict(lda_fit, newdata = test_data)

conf_mat <- table(
  Observed  = test_data$true_eco,
  Predicted = lda_pred$class
)

conf_mat

accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
accuracy

prop.table(table(test_data$true_eco))
#My classifier performs worst that a naive majority-class classifier

table(test_data$true_eco) / nrow(test_data)


n_perm <- 1000
perm_acc <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  perm_labels <- sample(train_data$true_eco)
  
  lda_perm <- lda(perm_labels ~ ., data = train_data)
  pred_perm <- predict(lda_perm, newdata = test_data)$class
  
  perm_acc[i] <- mean(pred_perm == test_data$true_eco)
}

p_value <- mean(perm_acc >= accuracy)

p_value

# -----------------------------
# Settings
# -----------------------------
terms_keep <- c("true_eco", "sex", "river", "true_eco:sex")

term_to_pred <- c(
  "true_eco"     = "Ecotype",
  "sex"          = "Sex",
  "river"        = "River",
  "true_eco:sex" = "Ecotype×Sex"
)

pred_levels <- c("Ecotype", "Sex", "River", "Ecotype×Sex")

trait_label <- function(x) {
  x %>%
    str_remove("_res$") %>%
    recode(
      "middle_gill_raker_length" = "MGRL",
      "gill_raker_count" = "GRC"
    )
}

dodge_width <- 0.7

# -----------------------------
# Build eta_df from MANOVA univariate tables
# -----------------------------
aov_list <- summary.aov(manova_fit)

eta_df <- purrr::imap_dfr(aov_list, function(tab, resp_name) {
  df <- as.data.frame(tab) %>%
    mutate(
      Term  = str_squish(rownames(tab)),                       # trim trailing spaces
      Trait = str_replace(resp_name, "^ Response ", "")
    ) %>%
    filter(Term %in% c(terms_keep, "Residuals")) %>%
    dplyr::select(Trait, Term, `Sum Sq`, `Df`, `F value`, `Pr(>F)`)
  
  ss_resid <- df$`Sum Sq`[df$Term == "Residuals"][1]
  
  df %>%
    filter(Term != "Residuals") %>%
    mutate(
      eta2_partial = `Sum Sq` / (`Sum Sq` + ss_resid),
      Predictor    = unname(term_to_pred[Term])
    )
})

alpha_sig <- 0.05

eta_df <- eta_df %>%
  mutate(
    Significance = factor(
      if_else(`Pr(>F)` < alpha_sig,
              "p < 0.05",
              "p ≥ 0.05"),
      levels = c("p < 0.05", "p ≥ 0.05")
    )
  )



# -----------------------------
# Order traits by decreasing Ecotype eta^2
# -----------------------------
eco_order <- eta_df %>%
  filter(Predictor == "Ecotype") %>%
  arrange(desc(eta2_partial)) %>%
  pull(Trait)

# -----------------------------
# Prep plotting data (clean labels + manual x positions)
# -----------------------------

pred_levels <- c("Ecotype", "Ecotype×Sex", "Sex", "River")

plot_df <- eta_df %>%
  mutate(
    Trait     = factor(Trait, levels = eco_order),
    Trait_lab = factor(trait_label(as.character(Trait)),
                       levels = trait_label(eco_order)),
    Predictor = factor(Predictor, levels = pred_levels),
    x_base    = as.numeric(Trait_lab),
    n_pred    = nlevels(Predictor),
    x_pos     = x_base + ((as.numeric(Predictor) - (n_pred + 1)/2) *
                            (dodge_width / n_pred))
  )

# -----------------------------
# Plot: lollipop (vertical sticks)
# -----------------------------
# ggplot(plot_df, aes(y = eta2_partial, color = Predictor)) +
#   
#   geom_segment(aes(x = x_pos,
#                    xend = x_pos,
#                    y = 0,
#                    yend = eta2_partial,
#                    linetype = Significance),
#                linewidth = 1) +
#   
#   geom_point(aes(x = x_pos), size = 3) +
#   
#   scale_x_continuous(
#     breaks = seq_along(levels(plot_df$Trait_lab)),
#     labels = levels(plot_df$Trait_lab),
#     expand = expansion(mult = c(0.02, 0.02))
#   ) +
#   
#   scale_color_manual(values = c(
#     "Ecotype"     = "#D18F00",   # deeper amber
#     "Ecotype×Sex" = "#1B6F6A",   # dark teal
#     "Sex"         = "grey40",
#     "River"       = "grey65"
#   )) +
#   
#   scale_linetype_manual(
#     name = "Significance",
#     values = c(
#       "p < 0.05" = "solid",
#       "p ≥ 0.05" = "22"
#     )
#   ) +
#   
#   theme_bw(base_size = 14) +
#   theme(
#     panel.grid = element_blank(),
#     axis.text.x = element_text(angle = 60, hjust = 1),
#     legend.spacing.y = unit(2, "pt"),
#     legend.key.height = unit(10, "pt"),
#     legend.box.spacing = unit(3, "pt"),
#     legend.background = element_blank(),
#     legend.position = "none"
#   ) +
#   ylab(expression("Partial " * eta^2)) +
#   xlab("")

# ggsave("figures/01_linear_trait_effect.png", unit = "cm", width = 22, height = 8, dpi = 600, )


ggplot(plot_df, aes(y = eta2_partial, color = Predictor)) +
  
  geom_segment(aes(x = x_pos,
                   xend = x_pos,
                   y = 0,
                   yend = eta2_partial),
               linewidth = 1) +
  
  geom_point(aes(x = x_pos, shape = Significance), size = 3) +
  
  scale_x_continuous(
    breaks = seq_along(levels(plot_df$Trait_lab)),
    labels = levels(plot_df$Trait_lab),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  scale_color_manual(values = c(
    "Ecotype"     = "#D18F00",   # deeper amber
    "Ecotype×Sex" = "#1B6F6A",   # dark teal
    "Sex"         = "grey40",
    "River"       = "grey65"
  )) +
  
  scale_shape_manual(
    name = "Significance",
    values = c(
      "p < 0.05" = 8,
      "p ≥ 0.05" = 16
    )
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.spacing.y = unit(2, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.box.spacing = unit(3, "pt"),
    legend.background = element_blank(),
    legend.position = "none"
  ) +
  ylab(expression("Partial " * eta^2)) +
  xlab("")

ggsave("figures/01_linear_trait_effect.png", unit = "cm", width = 22, height = 8, dpi = 600, )
