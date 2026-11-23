rm(list = ls())


#################
#### Library ####
#################

library(MASS)        # load before tidyverse so dplyr::select() isn't masked
library(car)
library(tidyverse)
library(magrittr)
library(openxlsx)

source("99_morpho_helpers.R", local = TRUE)
dir.create("figures", showWarnings = FALSE)

N_PERM <- 1000

#################
### Load Data ###
#################

data <- read.csv("data/data_gillrak_geomorph_age.csv", header = TRUE)

data %<>%
  mutate(gill_raker_count = upper_gill_raker + lower_gill_raker,
         river = if_else(str_detect(sample_id, 'NOT'), "NOT", "RUP")) %>%
  filter(sex %in% c("F", "M")) %>%
  mutate(true_eco = factor(true_eco),
         sex      = factor(sex),
         river    = factor(river))

### CHANGED: traits by name ###################################################
# Lengths that scale with body size -> size-corrected
length_traits <- c("MXL", "MDL", "POL", "OOL", "PSL", "TTL", "DOL", "LUL",
                   "ANL", "CPL", "CPD", "PCL", "PVL", "BDD", "HDD",
                   "middle_gill_raker_length")

missing_cols <- setdiff(c(length_traits, "AMO", "long_fourche", "gill_raker_count"), names(data))
if (length(missing_cols) > 0) stop("Columns not found in the CSV: ", paste(missing_cols, collapse = ", "))

###################################################
# Size correction (residuals on log fork length)
###################################################

for (tr in length_traits) {
  m <- lm(log(data[[tr]]) ~ log(data$long_fourche), na.action = na.exclude)
  data[[paste0(tr, "_res")]] <- rstandard(m)
}

#######################################
# MANOVA (Type II)
#######################################

responses <- c("AMO",
               paste0(length_traits, "_res"),
               "gill_raker_count")

data_manova <- data %>%
  dplyr::select(sample_id, river, sex, true_eco, long_fourche,
                all_of(length_traits), all_of(responses)) %>%
  drop_na(all_of(responses))

cat("Fish in the MANOVA:", nrow(data_manova), "\n")
print(table(data_manova$true_eco, data_manova$sex))

Y <- as.matrix(data_manova[, responses])
mlm_fit <- lm(Y ~ true_eco * sex + river, data = data_manova)


manova_II <- car::Manova(mlm_fit, type = "II")
manova_multi <- wilks_table(manova_II)
print(manova_multi)
write_sheet(manova_multi, "MANOVA_multivariate")


uni <- map_dfr(responses, function(tr) {
  fit <- lm(reformulate(c("true_eco * sex", "river"), response = tr), data = data_manova)
  tab <- car::Anova(fit, type = "II")
  ss_res <- tab["Residuals", "Sum Sq"]
  tibble(Trait = tr,
         Term  = rownames(tab),
         SS    = tab[["Sum Sq"]],
         Df    = tab[["Df"]],
         F_value = tab[["F value"]],
         p_value = tab[["Pr(>F)"]]) %>%
    filter(Term != "Residuals") %>%
    mutate(eta2_partial = SS / (SS + ss_res))     # partial eta-squared
}) %>%
  group_by(Term) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%   # across the 18 traits
  ungroup()

print(uni %>% filter(Term == "true_eco") %>% arrange(desc(eta2_partial)), n = Inf)
write_sheet(uni, "MANOVA_linear_traits")

#####################################################
# Fig. 1c: partial eta-squared per trait and effect
#####################################################

term_to_pred <- c("true_eco" = "Ecotype", "sex" = "Sex",
                  "river" = "River", "true_eco:sex" = "Ecotype×Sex")
pred_levels <- c("Ecotype", "Ecotype×Sex", "Sex", "River")

trait_label <- function(x) {
  x %>%
    str_remove("_res$") %>%
    dplyr::recode("middle_gill_raker_length" = "MGRL",
                  "gill_raker_count" = "GRC")
}

alpha_sig   <- 0.05
dodge_width <- 0.7


eta_df <- uni %>%
  mutate(Predictor = unname(term_to_pred[Term]),
         Significance = factor(if_else(p_value < alpha_sig, "p < 0.05", "p ≥ 0.05"),
                               levels = c("p < 0.05", "p ≥ 0.05")))

eco_order <- eta_df %>%
  filter(Predictor == "Ecotype") %>%
  arrange(desc(eta2_partial)) %>%
  pull(Trait)

plot_df <- eta_df %>%
  mutate(Trait_lab = factor(trait_label(Trait), levels = trait_label(eco_order)),
         Predictor = factor(Predictor, levels = pred_levels),
         x_base    = as.numeric(Trait_lab),
         n_pred    = nlevels(Predictor),
         # nudge the 4 predictors side by side around each trait position
         x_pos     = x_base + ((as.numeric(Predictor) - (n_pred + 1) / 2) * (dodge_width / n_pred)))

fig1c <- ggplot(plot_df, aes(y = eta2_partial, color = Predictor)) +
  geom_segment(aes(x = x_pos, xend = x_pos, y = 0, yend = eta2_partial), linewidth = 1) +
  geom_point(aes(x = x_pos, shape = Significance), size = 3) +
  scale_x_continuous(breaks = seq_along(levels(plot_df$Trait_lab)),
                     labels = levels(plot_df$Trait_lab),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = c("Ecotype" = "#D18F00", "Ecotype×Sex" = "#1B6F6A",
                                "Sex" = "grey40", "River" = "grey65")) +
  scale_shape_manual(name = "Significance", values = c("p < 0.05" = 8, "p ≥ 0.05" = 16), drop = FALSE) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "none") +
  ylab(expression("Partial " * eta^2)) +
  xlab("")

fig1c
ggsave("figures/01_linear_trait_effect.png", fig1c, unit = "cm", width = 22, height = 8, dpi = 600)

#####################################################
# Single-trait plots for the ecotype-significant traits
#####################################################

sig_eco <- uni %>% filter(Term == "true_eco", p_value < alpha_sig) %>% pull(Trait)
cat("Traits with a significant ecotype effect:", paste(trait_label(sig_eco), collapse = ", "), "\n")

y_for_trait <- function(tr) {
  raw <- str_remove(tr, "_res$")
  if (raw %in% c("gill_raker_count", "AMO")) data_manova[[raw]] else data_manova[[raw]] / data_manova$long_fourche
}

for (tr in sig_eco) {
  st <- uni %>% filter(Trait == tr, Term == "true_eco")
  df <- data_manova %>%
    mutate(y = y_for_trait(tr),
           true_eco = fct_relevel(true_eco, "summer_run", "fall_run"))
  p <- ggplot(df, aes(x = true_eco, y = y)) +
    geom_boxplot(fill = "white", width = 0.3, notch = TRUE, outlier.alpha = 0) +
    geom_jitter(aes(color = true_eco), width = 0.1, alpha = 0.6) +
    scale_color_manual(values = c("#E69F00", "skyblue1")) +
    scale_x_discrete(labels = c("KUU", "NUT")) +
    xlab("") +
    ylab(if (tr %in% c("gill_raker_count", "AMO")) trait_label(tr) else paste(trait_label(tr), "/ fork length")) +
    annotate("text", x = 0.5, y = Inf, hjust = 0, vjust = 1.3, size = 4,
             label = paste0("F = ", round(st$F_value, 2),
                            "\npartial eta² = ", round(st$eta2_partial, 3))) +
    theme_bw(base_size = 14) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(),
          legend.position = "none")
  ggsave(paste0("figures/01_", trait_label(tr), ".png"), p,
         width = 10, height = 10, unit = "cm", dpi = 600)
}

#######################################
# LDA on the ecotype-significant traits
#######################################

if (length(sig_eco) == 0) stop("No trait has a significant ecotype effect, so there is nothing to put in the LDA.")
lda_X <- data_manova %>% dplyr::select(all_of(sig_eco))

lda_lin <- lda_loo(lda_X, data_manova$true_eco)
lda_lin$confusion
lda_lin$sensitivity
cat("Balanced accuracy:", round(lda_lin$balanced_accuracy, 3), "\n")

p_lin <- lda_perm_test(lda_X, data_manova$true_eco, lda_lin$balanced_accuracy, n_perm = N_PERM)
cat("Permutation p-value (vs chance = 0.5):", signif(p_lin, 3), "\n")

write_sheet(data.frame(traits = paste(trait_label(sig_eco), collapse = ", "),
                       n_fish = nrow(lda_X),
                       balanced_accuracy  = lda_lin$balanced_accuracy,
                       success_fall_run   = lda_lin$sensitivity[["fall_run"]],
                       success_summer_run = lda_lin$sensitivity[["summer_run"]],
                       permutation_p = p_lin,
                       n_permutations = N_PERM),
            "LDA_linear_traits")
