#------------------------------------------------------------------------------#
# 06_binomial_model_prob_kuukamek.R
# Seasonal variation in the probability of harvesting kâcikâsikumekw in Cree
# subsistence fisheries (binomial GLMM with community as a random effect).
#------------------------------------------------------------------------------#

rm(list = ls())

library(tidyverse)
library(magrittr)
library(lubridate)
library(lme4)
library(splines)
library(DHARMa)
library(performance)
library(ggeffects)
library(lattice)
library(openxlsx)

#--------------------#
#---- PARAMETERS ----#
#--------------------#

pofz_min <- 0.99
min_week <- 26      # earlier weeks contain too few individuals
ref_year <- 2021    # non-leap year used to display weeks as dates

#-------------------#
#---- FUNCTIONS ----#
#-------------------#

clean_id <- function(x) {
  x <- str_replace_all(x, "-", "_")
  str_replace(x, "([A-Za-z]*s)_(\\d+)_", function(m) {
    g <- str_match(m, "([A-Za-z]*s)_(\\d+)_")
    sprintf("%s_%03d_", g[, 2], as.integer(g[, 3]))
  })
}

# Mid-point date of a week as defined by lubridate::week()
# (week 1 = 1-7 January, week 2 = 8-14 January, ...)
week_to_date <- function(week, year = ref_year) {
  as.Date(paste0(year, "-01-01")) + 7 * (week - 1) + 3
}

week_range_label <- function(week, year = ref_year) {
  start <- as.Date(paste0(year, "-01-01")) + 7 * (week - 1)
  paste(format(start, "%d %B"), "-", format(start + 6, "%d %B"))
}

write_supp_sheet <- function(data, sheet, file = "res/supplementary_table.xlsx") {
  wb <- if (file.exists(file)) loadWorkbook(file) else createWorkbook()
  if (sheet %in% names(wb)) removeWorksheet(wb, sheet)
  addWorksheet(wb, sheet)
  writeData(wb, sheet, data)
  saveWorkbook(wb, file, overwrite = TRUE)
}

# Marginal and conditional R2 (Nakagawa & Schielzeth 2013) including the
# random-slope variance (Johnson 2014). The random-effect variance is the mean of
# z' Sigma z over observations, for each grouping term.
r2_random_slopes <- function(model, data) {
  z_all <- cbind("(Intercept)" = 1, week_c = data$week_c)
  var_fixed <- var(predict(model, re.form = NA))
  var_random <- sum(map_dbl(VarCorr(model), function(S) {
    z <- z_all[, colnames(S), drop = FALSE]
    mean(rowSums((z %*% S) * z))
  }))
  var_resid <- pi^2 / 3
  c(R2_marginal = var_fixed / (var_fixed + var_random + var_resid),
    R2_conditional = (var_fixed + var_random) / (var_fixed + var_random + var_resid))
}

# Most parsimonious model among those within 2 AIC units of the lowest AIC
select_parsimonious <- function(tab) {
  tab %>%
    filter(delta_AIC < 2) %>%
    arrange(df, AIC) %>%
    slice(1) %>%
    pull(model)
}

#--------------#
#---- DATA ----#
#--------------#

assignment_res <- read.table("data/migrating_individuals_assignment_rubias.txt") %>%
  filter(PofZ > pofz_min) %>%
  mutate(community = factor(community, levels = c("Whapmagoostui", "Chisasibi", "Wemindji",
                                                  "Eastmain", "Waskaganish")),
         collection = factor(collection, levels = c("V1", "V2", "V3"),
                             labels = c("Fall-run north", "Fall-run south", "Kâcikâsikumekw")))

metadata <- read.csv("data/migrating_samples_metadata.csv") %>%
  transmute(indiv = clean_id(fishes_id), date = as.Date(date, format = "%d-%m-%Y")) %>%
  distinct(indiv, .keep_all = TRUE)

assignment_res_date <- assignment_res %>%
  left_join(metadata, by = "indiv") %>%
  mutate(week_of_year = week(date))

message("Assigned individuals without sampling date: ", sum(is.na(assignment_res_date$date)))
message(sprintf("Individuals caught before week %d (excluded): %d", min_week,
                sum(assignment_res_date$week_of_year < min_week, na.rm = TRUE)))

assignment_res_date_stats <- assignment_res_date %>%
  filter(!is.na(week_of_year), week_of_year >= min_week) %>%
  mutate(is_kacik = as.integer(collection == "Kâcikâsikumekw"),
         week_c = week_of_year - mean(week_of_year))

mean_week <- mean(assignment_res_date_stats$week_of_year)

message("Individuals in the model: ", nrow(assignment_res_date_stats))
assignment_res_date_stats %>% count(collection)
assignment_res_date_stats %>% count(community)

fit_glmm <- function(formula) {
  glmer(formula, data = assignment_res_date_stats, family = binomial,
        control = glmerControl(optimizer = "bobyqa"))
}

#----------------------------------------#
#---- FIXED EFFECTS: EFFECT OF WEEK ----#
#----------------------------------------#

models <- list(
  null      = fit_glmm(is_kacik ~ 1 + (1 | community)),
  linear    = fit_glmm(is_kacik ~ week_c + (1 | community)),
  quadratic = fit_glmm(is_kacik ~ week_c + I(week_c^2) + (1 | community))
)

for (k in 2:5) {
  models[[paste0("spline_df", k)]] <-
    fit_glmm(as.formula(sprintf("is_kacik ~ ns(week_c, df = %d) + (1 | community)", k)))
}

# Residual diagnostics for each candidate model
diagnostics <- map_dfr(names(models), function(m) {
  sim <- simulateResiduals(models[[m]], n = 1000, plot = FALSE)
  tibble(model = m,
         ks_p = testUniformity(sim, plot = FALSE)$p.value,
         dispersion_p = testDispersion(sim, plot = FALSE)$p.value,
         outlier_p = testOutliers(sim, plot = FALSE)$p.value)
})

aic_comparison <- tibble(
  model = names(models),
  df = map_int(models, ~ as.integer(attr(logLik(.x), "df"))),
  AIC = map_dbl(models, AIC),
  BIC = map_dbl(models, BIC)) %>%
  mutate(delta_AIC = AIC - min(AIC),
         delta_BIC = BIC - min(BIC)) %>%
  left_join(diagnostics, by = "model") %>%
  arrange(AIC)

best_fixed <- select_parsimonious(aic_comparison)
aic_comparison %<>% mutate(selected = model == best_fixed)

print(aic_comparison)
message("Selected fixed-effect structure: ", best_fixed)

write_supp_sheet(aic_comparison, "binomial_model_mixed_stock")

mod_fixed <- models[[best_fixed]]

#--------------------------------------------------------#
#---- RANDOM EFFECTS: COMMUNITY-SPECIFIC TIMING ----#
#--------------------------------------------------------#

# Random intercept only, plus a random slope for week by community with and
# without the intercept-slope correlation.
fixed_formula <- update(formula(mod_fixed), . ~ . - (1 | community))

random_models <- list(
  intercept          = mod_fixed,
  slope_uncorrelated = fit_glmm(update(fixed_formula, . ~ . + (1 | community) + (0 + week_c | community))),
  slope_correlated   = fit_glmm(update(fixed_formula, . ~ . + (1 + week_c | community)))
)

walk2(names(random_models), random_models, function(name, m) {
  cat("\n", name, "\n", sep = "")
  print(VarCorr(m))
})

lrt_vs_intercept <- function(m) {
  if (identical(m, random_models$intercept)) return(tibble(LRT_chisq = NA_real_, LRT_df = NA_real_, LRT_p = NA_real_))
  a <- anova(random_models$intercept, m)
  tibble(LRT_chisq = a$Chisq[2], LRT_df = a$Df[2], LRT_p = a$`Pr(>Chisq)`[2])
}

random_comparison <- tibble(
  model = names(random_models),
  random_effects = c("(1 | community)",
                     "(1 | community) + (0 + week_c | community)",
                     "(1 + week_c | community)"),
  df = map_int(random_models, ~ as.integer(attr(logLik(.x), "df"))),
  AIC = map_dbl(random_models, AIC),
  BIC = map_dbl(random_models, BIC),
  singular = map_lgl(random_models, isSingular)) %>%
  bind_cols(map_dfr(random_models, lrt_vs_intercept)) %>%
  mutate(delta_AIC = AIC - min(AIC))

best_random <- select_parsimonious(random_comparison)
random_comparison %<>% mutate(selected = model == best_random)

print(random_comparison)
message("Selected random-effect structure: ", best_random)

write_supp_sheet(random_comparison, "binomial_random_effects")

mod_final <- random_models[[best_random]]

#-----------------------#
#---- FINAL MODEL ----#
#-----------------------#

summary(mod_final)

sim_resid_final <- simulateResiduals(mod_final, n = 1000, plot = TRUE)
testUniformity(sim_resid_final)
testDispersion(sim_resid_final)
testOutliers(sim_resid_final)
testZeroInflation(sim_resid_final)
plotResiduals(sim_resid_final, form = assignment_res_date_stats$week_c)
plotResiduals(sim_resid_final, form = assignment_res_date_stats$community)

check_model(mod_final)
model_performance(mod_final)

r2_final <- r2_random_slopes(mod_final, assignment_res_date_stats)
cat("\nMarginal R2 (fixed effects):", round(r2_final[["R2_marginal"]], 3))
cat("\nConditional R2 (fixed + random):", round(r2_final[["R2_conditional"]], 3), "\n")

dotplot(ranef(mod_final, condVar = TRUE))

#-----------------------------------------#
#---- POPULATION-LEVEL PREDICTIONS ----#
#-----------------------------------------#

# Population-level predictions with 95% confidence intervals (fixed effects)
pred_df <- ggpredict(mod_final, terms = "week_c [all]") %>%
  as.data.frame() %>%
  mutate(week = x + mean_week,
         date = week_to_date(week))

week_summary <- assignment_res_date_stats %>%
  group_by(week_of_year) %>%
  summarise(prop_kacik = mean(is_kacik),
            n = n(),
            .groups = "drop") %>%
  mutate(date = week_to_date(week_of_year))

peak <- pred_df %>% slice_max(predicted, n = 1, with_ties = FALSE)
message(sprintf("Peak: week %d (%s), predicted probability = %.2f (95%% CI %.2f-%.2f)",
                round(peak$week), week_range_label(round(peak$week)),
                peak$predicted, peak$conf.low, peak$conf.high))

p1 <- ggplot() +
  geom_ribbon(data = pred_df, aes(x = date, ymin = conf.low, ymax = conf.high),
              fill = "#E69F00", alpha = 0.2) +
  geom_line(data = pred_df, aes(x = date, y = predicted),
            color = "#E69F00", linewidth = 1.3) +
  geom_point(data = week_summary, aes(x = date, y = prop_kacik, size = n),
             alpha = 0.5, color = "gray30") +
  scale_size_continuous(range = c(2, 6), name = "Sample size") +
  scale_x_date(date_labels = "%B", date_breaks = "1 month", expand = c(0.02, 0.02)) +
  labs(x = "Month", y = "Probability of harvesting kâcikâsikumekw") +
  theme_bw(base_size = 16) +
  theme(legend.position = c(0.8, 0.8),
        legend.background = element_rect(fill = NA),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black"),
        axis.text.x = element_text(hjust = 0.5))

print(p1)

ggsave("figures/06_kuukamek_probability_by_week.png", p1, width = 5, height = 5, dpi = 600)

#--------------------------------------------#
#---- COMMUNITY-LEVEL PEAK OF HARVEST ----#
#--------------------------------------------#

# Community-specific predictions, restricted to the weeks sampled in each community
community_weeks <- assignment_res_date_stats %>%
  group_by(community) %>%
  summarise(first_week = min(week_of_year),
            last_week = max(week_of_year),
            n = n(),
            .groups = "drop")

week_grid <- community_weeks %>%
  reframe(week_of_year = seq(first_week, last_week, by = 0.1), .by = community) %>%
  mutate(week_c = week_of_year - mean_week)

week_grid$predicted <- predict(mod_final, newdata = week_grid, type = "response")

community_peaks <- week_grid %>%
  group_by(community) %>%
  slice_max(predicted, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(peak_week = round(week_of_year),
         peak_dates = week_range_label(peak_week)) %>%
  left_join(community_weeks, by = "community") %>%
  dplyr::select(community, n, first_week, last_week, peak_week, peak_dates,
                peak_probability = predicted)

print(community_peaks)

write_supp_sheet(community_peaks, "binomial_community_peaks")
