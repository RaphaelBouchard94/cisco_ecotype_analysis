rm(list = ls())

###############################################################################
# 03 - Life history: age comparison, back-calculated growth, von Bertalanffy
#      mixed models (Fig. 2c, Supplementary Table 4)
#
# CHANGES vs your original script
#  1. Back-calculation now uses the BIOLOGICAL INTERCEPT method (Campana 1990)
#     with literature values for cisco (see "Biological intercept" below).
#     Your previous body-proportional method (v = 0.23, estimated from adults
#     only) gave age-1 lengths of ~250-300 mm. Set BACKCALC <- "BPH" to get it back.
#  2. The last otolith measurement of each fish is checked: on your data it is
#     the otolith EDGE (n measurements = age + 1). It is now used as the radius
#     at capture and no longer plotted/modelled as an extra year of growth.
#  3. Sensitivity analysis: refits the growth model across the published range
#     of length at hatch and a range of otolith radius at hatch, so you can show
#     the results don't depend on the exact values.
#  4. Best model: chosen by a stated rule (BEST_MODEL = "auto": among models
#     within 2 AIC units of the lowest, the one with fewest parameters), or
#     set by hand. Curves, CIs and phi' work for any of the 8 models.
#  5. Growth performance phi' = log10(K) + 2*log10(Linf). The old code used
#     log10(K + K) for the summer run (0.30 too high). Fixed.
#  6. Likelihood-ratio tests (as the Methods say) + Wald tests from the full model.
#  7. Age and length-at-age models on the real scale (years, mm).
#  8. Starting values are checked and fall back to safe values if FSA's are bad.
###############################################################################

########################
# Library
########################

library(MASS)     # mvrnorm(); load before tidyverse so dplyr::select() isn't masked
library(nlme)
library(FSA)
library(tidyverse)
library(magrittr)
library(openxlsx)

source("99_morpho_helpers.R", local = TRUE)
dir.create("figures", showWarnings = FALSE)

########################
# Settings
########################

# Is the last otolith measurement of each fish the otolith EDGE (radius at
# capture)? NA = decide from the data (compares number of measurements with
# otolith age). On your data this comes out TRUE.
LAST_RING_IS_EDGE <- NA

# ---- Back-calculation method ----
#   "BI"  = biological intercept (Campana 1990)   <- default
#   "BPH" = body proportional hypothesis, Li = (Si/S)^v * L (previous method)
BACKCALC <- "BI"

# ---- Biological intercept: length and otolith radius at hatch ----
# L0 = total length at hatch (mm). At hatch the caudal fin isn't forked, so
#      total length = fork length.
#   * 10.69 +/- 0.42 mm (mean +/- SD): Lake Superior cisco reared in the lab
#     (Oyadomari & Auer 2007, Trans. Am. Fish. Soc. 136: 766-777).
#   * 8.68-10.21 mm (Lake Superior) and 9.75-11.32 mm (Lake Ontario) cisco
#     incubated at 9 to 2 degC; larvae from colder incubation hatch larger
#     (Stewart et al. 2021, "Influence of warming temperatures on coregonine
#     embryogenesis within and among species").
#   James Bay eggs incubate under ice near 0-2 degC, i.e. at the cold end, so
#   we use 10.69 mm and test 8.68-11.32 mm in the sensitivity analysis.
BI_L0_MM <- 10.69

# R0 = otolith radius at hatch (same units as the X..._um columns).
#   * Otoliths at hatch were 65.6 +/- 4.6 um long in Oyadomari & Auer (2007)
#     (same Lake Superior cisco), so radius ~ 65.6 / 2 = 32.8 um. No hatch check
#     forms in cisco (the first check appears ~28 d after hatch, enclosing a
#     region 63.1 +/- 5.6 um long), so you can't measure this on your sections.
#   This is half the otolith length, not necessarily the radius along your
#   reading axis, so the sensitivity analysis tests 20-50 um.
BI_R0_UM <- 32.8

# Ranges for the sensitivity analysis (only run when BACKCALC == "BI")
SENS_L0 <- c(8.68, 10.69, 11.32)
SENS_R0 <- c(20, 32.8, 50)

# ---- Growth model choice ----
# "auto" = among models within 2 AIC units of the lowest AIC, take the one with
# the fewest parameters (the usual parsimony rule). Or name one: "mL", "mLK", ...
BEST_MODEL <- "auto"

N_BOOT <- 1000    # parametric bootstrap draws for the confidence bands
set.seed(123)

########################
# Load Data
########################

data <- read.table("data/data_gillrak_geomorph_age_morphres.csv", header = TRUE)

data$ecotype <- factor(data$true_eco)
data$river   <- if_else(str_detect(data$sample_id, "NOT"), "NOT", "RUP")

#############################################
# 1. Age at reproduction: Welch's t-test
#############################################

age_summary <- data %>%
  filter(!is.na(true_age)) %>%
  group_by(ecotype) %>%
  summarise(n = n(), mean_age = mean(true_age), sd_age = sd(true_age))
print(age_summary)

age_test <- t.test(true_age ~ ecotype, data = data)   # Welch's test by default
age_test

age_out <- tibble(mean_fall   = age_test$estimate[[1]],
                  mean_summer = age_test$estimate[[2]],
                  difference  = age_test$estimate[[2]] - age_test$estimate[[1]],
                  t = age_test$statistic[[1]], df = age_test$parameter[[1]],
                  p_value = age_test$p.value)
write_sheet(age_out, "Age_comparison")

#############################################
# 2. Are summer-run fish larger at a given age?
#############################################

# This result needs no back-calculation: otolith age + length at capture.
len_age     <- lm(long_fourche ~ true_age + ecotype, data = data)
len_age_int <- lm(long_fourche ~ true_age * ecotype, data = data)
summary(len_age)              # ecotypesummer_run = extra mm at the same age
anova(len_age, len_age_int)   # does the age slope differ between ecotypes?

write_sheet(as.data.frame(summary(len_age)$coefficients) %>% rownames_to_column("Term"),
            "Length_given_age")

##########################
# 3. Otolith measurements
##########################

# Otolith radius columns X1_um, X2_um, ... sorted by their number
ring_cols <- grep("^X[0-9]+_um$", names(data), value = TRUE)
ring_cols <- ring_cols[order(as.integer(str_extract(ring_cols, "[0-9]+")))]

bc <- data %>%
  filter(!is.na(long_fourche)) %>%
  mutate(across(all_of(ring_cols), ~ replace_na(.x, 0)))

R <- as.matrix(bc[, ring_cols])
bc$n_rings   <- rowSums(R > 0)                                    # number of measurements
bc$last_ring <- apply(R > 0, 1, function(r) if (any(r)) max(which(r)) else 0)
bc$total_rad <- do.call(pmax, as.data.frame(R))                   # largest = outermost
bc <- bc %>% filter(total_rad > 0)

# Edge or last annulus?
ring_check <- bc %>% filter(!is.na(true_age)) %>% count(measurements_minus_age = n_rings - true_age)
print(ring_check)

if (is.na(LAST_RING_IS_EDGE)) {
  d <- bc$n_rings - bc$true_age
  LAST_RING_IS_EDGE <- if (all(is.na(d))) FALSE else
    mean(d == 1, na.rm = TRUE) > mean(d == 0, na.rm = TRUE)
}
if (LAST_RING_IS_EDGE) {
  message("Last measurement treated as the otolith EDGE (radius at capture).")
} else {
  warning("Last measurement treated as the last ANNULUS: the edge radius isn't in the data, ",
          "so growth after the last annulus is ignored.")
}

# Fish length vs otolith radius. v is only used by the BPH method, but the
# no-interaction test is worth reporting either way.
lm_eco <- lm(log(long_fourche) ~ log(total_rad) + ecotype, data = bc)
lm_int <- lm(log(long_fourche) ~ log(total_rad) * ecotype, data = bc)
anova(lm_eco, lm_int)
v <- unname(coef(lm_eco)["log(total_rad)"])
cat("Body-proportional exponent v =", round(v, 3), "(used only if BACKCALC = \"BPH\")\n")

##########################
# 4. Back-calculation
##########################

# backcalc() returns one row per fish per annulus with the back-calculated length.
#   Si = radius at annulus i (radius), S = radius at capture (total_rad),
#   L  = fork length at capture (long_fourche)
#   BI : Li = L + (Si - S) * (L - L0) / (S - R0)
#        = a straight line, for each fish, through its capture point (S, L)
#          and the hatch point (R0, L0) shared by all fish.
#   BPH: Li = (Si / S)^v * L
backcalc <- function(bc, method = BACKCALC, L0 = BI_L0_MM, R0 = BI_R0_UM) {
  bc %>%
    dplyr::select(sample_id, true_eco, river, long_fourche, total_rad, last_ring, all_of(ring_cols)) %>%
    pivot_longer(all_of(ring_cols), names_to = "ring", values_to = "radius") %>%
    mutate(age_at_length = as.integer(str_extract(ring, "[0-9]+"))) %>%
    filter(radius > 0) %>%
    filter(!(LAST_RING_IS_EDGE & age_at_length == last_ring)) %>%    # drop the edge point
    mutate(length_at_age = if (method == "BI") {
                             long_fourche + (radius - total_rad) * (long_fourche - L0) / (total_rad - R0)
                           } else {
                             (radius / total_rad)^v * long_fourche
                           },
           true_eco = factor(true_eco))
}

cat("Back-calculation method:", BACKCALC,
    if (BACKCALC == "BI") paste0("(L0 = ", BI_L0_MM, " mm, R0 = ", BI_R0_UM, " um)") else "", "\n")
length_at_age <- backcalc(bc)

# Sanity check: mean back-calculated length at ages 1-3
laa_check <- length_at_age %>%
  filter(age_at_length <= 3) %>%
  group_by(true_eco, age_at_length) %>%
  summarise(n_fish = n(), mean_length = mean(length_at_age),
            min = min(length_at_age), max = max(length_at_age), .groups = "drop")
print(laa_check)

write_sheet(data.frame(method = BACKCALC, L0_mm = BI_L0_MM, R0_um = BI_R0_UM, v = v,
                       last_measurement_is_edge = LAST_RING_IS_EDGE,
                       p_radius_x_ecotype = anova(lm_eco, lm_int)$`Pr(>F)`[2]),
            "Backcalc_settings")
write_sheet(laa_check, "Backcalc_length_age1to3")

cat("Fish in the growth models (for the Fig. 2c caption):\n")
print(length_at_age %>% distinct(sample_id, true_eco) %>% count(true_eco))

###########################
# 5. von Bertalanffy mixed models
###########################

vb <- function(age, Linf, K, t0) Linf * (1 - exp(-K * (age - t0)))

# ---- Starting values (only used to start the optimiser) ----
# Each candidate is checked (Linf, K > 0 and all finite); on your data
# findGrowthStarts() returned K < 0 and t0 = NaN, so we try several sources.
starts_ok <- function(s) {
  !is.null(s) && all(c("Linf", "K", "t0") %in% names(s)) &&
    all(is.finite(s[c("Linf", "K", "t0")])) && s[["Linf"]] > 0 && s[["K"]] > 0
}
get_starts <- function(df) {
  f <- length_at_age ~ age_at_length
  cand <- list(
    vbStarts = function() unlist(suppressWarnings(FSA::vbStarts(f, data = df))),
    findGrowthStarts = function() FSA::findGrowthStarts(f, data = df, type = "von Bertalanffy"),
    nls_means = function() {
      m <- df %>% group_by(age_at_length) %>% summarise(L = mean(length_at_age), .groups = "drop")
      coef(nls(L ~ Linf * (1 - exp(-K * (age_at_length - t0))), data = m,
               start = list(Linf = max(m$L) * 1.1, K = 0.3, t0 = -1)))
    },
    generic = function() c(Linf = max(df$length_at_age) * 1.1, K = 0.3, t0 = -1)
  )
  for (nm in names(cand)) {
    s <- tryCatch(cand[[nm]](), error = function(e) NULL)
    if (starts_ok(s)) {
      message("Starting values from ", nm, ".")
      return(s[c("Linf", "K", "t0")])
    }
  }
}

# ---- One model = ecotype affects the parameters listed in eco_on ----
#   fixed: list of 3 formulas, e.g. Linf ~ true_eco, K ~ 1, t0 ~ 1
#   start: in the same order as the coefficients:
#          Linf, [Linf ecotype effect], K, [K effect], t0, [t0 effect]
#          (an ecotype effect starts at 0 = "no difference")
#   random: Linf varies among rivers and among fish within rivers (Methods).
# nlme() fits by maximum likelihood by default, so AIC and LRTs are valid.
# The VB equation is written out (not vb()) because nlme() can't find
# functions defined inside a script run with 00_run_all.R.
vb_model <- function(eco_on, df, sv) {
  pars  <- c("Linf", "K", "t0")
  fixed <- lapply(pars, function(p) as.formula(paste(p, if (p %in% eco_on) "~ true_eco" else "~ 1")))
  start <- unlist(lapply(pars, function(p) if (p %in% eco_on) c(sv[[p]], 0) else sv[[p]]))
  tryCatch(
    nlme(length_at_age ~ Linf * (1 - exp(-K * (age_at_length - t0))),
         data = df, fixed = fixed, random = Linf ~ 1 | river/sample_id,
         start = start, na.action = na.omit),
    error = function(e) { warning("Model with ecotype on {", paste(eco_on, collapse = ", "),
                                  "} did not converge: ", conditionMessage(e)); NULL })
}

model_specs <- list(m0 = character(0), mL = "Linf", mK = "K", mt0 = "t0",
                    mLK = c("Linf", "K"), mLt0 = c("Linf", "t0"),
                    mKt0 = c("K", "t0"), mLKt0 = c("Linf", "K", "t0"))

sv0 <- get_starts(length_at_age)
sv0
models <- compact(map(model_specs, vb_model, df = length_at_age, sv = sv0))

if (length(models) == 0) stop("None of the growth models converged; see warnings().")
if (!"mLKt0" %in% names(models)) warning("The full model (mLKt0) didn't converge.")

# ---- AIC model selection (Supplementary Table 4) ----
model_summary <- map_dfr(names(models), function(m) {
  ll <- logLik(models[[m]])
  tibble(Model = m, n_par = attr(ll, "df"), AIC = AIC(models[[m]]), BIC = BIC(models[[m]]),
         logLik = as.numeric(ll))
}) %>%
  mutate(delta_AIC = AIC - min(AIC)) %>%
  arrange(AIC)
print(model_summary)
write_sheet(model_summary, "Model_selection_AIC")

# Parsimony rule: simplest model among those within 2 AIC units of the best
best <- if (BEST_MODEL == "auto") {
  model_summary %>% filter(delta_AIC <= 2) %>% arrange(n_par, AIC) %>% slice(1) %>% pull(Model)
} else BEST_MODEL
cat("Model used for Fig. 2c and parameters:", best,
    "(ecotype affects:", paste(model_specs[[best]], collapse = " + "), ")\n")

# ---- Wald tests in the full model ----
if ("mLKt0" %in% names(models)) {
  full_tab <- as.data.frame(summary(models$mLKt0)$tTable) %>% rownames_to_column("Parameter")
  print(full_tab)
  write_sheet(full_tab, "Full_model_Wald")
}

# ---- Likelihood-ratio tests (Methods: "compared ... using likelihood-ratio tests") ----
lrt <- function(simple, complex, label) {
  if (!all(c(simple, complex) %in% names(models))) return(NULL)
  a <- anova(models[[simple]], models[[complex]])
  tibble(Comparison = label, L_ratio = a$L.Ratio[2], df = diff(a$df), p_value = a$`p-value`[2])
}
lrt_tab <- bind_rows(
  lrt("m0", "mL",   "Ecotype effect on Linf (m0 vs mL)"),
  lrt("mL", "mLK",  "Adding ecotype effect on K (mL vs mLK)"),
  lrt("mL", "mLt0", "Adding ecotype effect on t0 (mL vs mLt0)")
)
print(lrt_tab)
write_sheet(lrt_tab, "LRT_tests")

# ---- Best model parameters ----
m_best <- models[[best]]
fixef_tab <- as.data.frame(summary(m_best)$tTable) %>%
  rownames_to_column("Parameter") %>%
  dplyr::rename(Estimate = Value, SE = Std.Error, t_value = `t-value`, p_value = `p-value`) %>%
  mutate(across(Estimate:p_value, ~ round(.x, 4)))
print(fixef_tab)
write_sheet(fixef_tab, "Best_model_parameters")

#########################################
# 6. Population growth curves (Fig. 2c)
#########################################

# par_value(): the value of parameter p (Linf, K or t0) for one ecotype, from a
# vector of fixed effects (one fit) or a matrix (one row per bootstrap draw).
# nlme names a parameter "K" when it is shared, and "K.(Intercept)" +
# "K.true_ecosummer_run" when it differs between ecotypes.
par_value <- function(beta, p, eco) {
  nms <- if (is.matrix(beta)) colnames(beta) else names(beta)
  col <- function(nm) if (is.matrix(beta)) beta[, nm] else beta[[nm]]
  base <- if (p %in% nms) col(p) else col(paste0(p, ".(Intercept)"))
  eff  <- paste0(p, ".true_ecosummer_run")
  if (eco == "summer_run" && eff %in% nms) base + col(eff) else base
}

ecos <- c("fall_run", "summer_run")
fe <- fixef(m_best)
pars <- map_dfr(ecos, function(e) tibble(ecotype = e,
                                         Linf = par_value(fe, "Linf", e),
                                         K    = par_value(fe, "K", e),
                                         t0   = par_value(fe, "t0", e)))
print(pars)
cat("Delta Linf (summer - fall):", round(diff(pars$Linf), 1), "mm\n")
cat("Predicted length at age 1:", paste(ecos, round(vb(1, pars$Linf, pars$K, pars$t0)), "mm", collapse = "; "), "\n")

max_age  <- max(length_at_age$age_at_length)
age_seq  <- seq(0, max_age, length.out = 200)
pred_pop <- expand.grid(age = age_seq, ecotype = ecos, stringsAsFactors = FALSE) %>%
  left_join(pars, by = "ecotype") %>%
  mutate(length = vb(age, Linf, K, t0))

# Parametric bootstrap: draw fixed effects from their estimated sampling
# distribution (multivariate normal), recompute the curve for each draw and
# take the 2.5% / 97.5% quantiles. Population curve uncertainty only.
beta_sim <- mvrnorm(N_BOOT, mu = fe, Sigma = vcov(m_best))
sim_pars <- map(purrr::set_names(ecos), function(e)     # purrr:: because magrittr masks set_names()
  list(Linf = par_value(beta_sim, "Linf", e), K = par_value(beta_sim, "K", e), t0 = par_value(beta_sim, "t0", e)))

ci_pop <- expand.grid(age = age_seq, ecotype = ecos, stringsAsFactors = FALSE) %>%
  mutate(LCI = NA_real_, UCI = NA_real_)
for (i in seq_len(nrow(ci_pop))) {
  sp <- sim_pars[[ci_pop$ecotype[i]]]
  sims <- vb(ci_pop$age[i], sp$Linf, sp$K, sp$t0)
  ci_pop$LCI[i] <- quantile(sims, 0.025, na.rm = TRUE)
  ci_pop$UCI[i] <- quantile(sims, 0.975, na.rm = TRUE)
}

fig2c <- ggplot() +
  geom_point(data = length_at_age, aes(age_at_length, length_at_age, color = true_eco), alpha = 0.3) +
  geom_ribbon(data = ci_pop, aes(age, ymin = LCI, ymax = UCI, fill = ecotype), alpha = 0.25, color = NA) +
  geom_line(data = pred_pop, aes(age, length, color = ecotype), linewidth = 1.2) +
  scale_color_manual(values = c("fall_run" = "skyblue", "summer_run" = "#E69F00")) +
  scale_fill_manual(values = c("fall_run" = "skyblue", "summer_run" = "#E69F00")) +
  scale_x_continuous(breaks = seq(0, max_age, by = 1)) +
  theme_bw(base_size = 16) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(),
        legend.position = "none") +
  labs(x = "Age", y = "Fork length (mm)")
fig2c
ggsave("figures/03_vb_growth_rates.png", fig2c, width = 10, height = 10, unit = "cm", dpi = 900)

##############################################
# 7. Growth performance index phi' (FIXED)
##############################################

# phi' = log10(K) + 2 * log10(Linf)   (Pauly & Munro 1984; unitless).
# If K is shared in the best model, delta-phi' = 2*log10(Linf ratio): it just
# restates the Linf difference.
phi_sim <- tibble(
  fall   = log10(sim_pars$fall_run$K)   + 2 * log10(sim_pars$fall_run$Linf),
  summer = log10(sim_pars$summer_run$K) + 2 * log10(sim_pars$summer_run$Linf)
) %>%
  mutate(delta = summer - fall)

phi_out <- map_dfr(names(phi_sim), function(nm) {
  q <- quantile(phi_sim[[nm]], c(0.025, 0.5, 0.975), na.rm = TRUE)
  tibble(quantity = nm, median = q[[2]], lower_95 = q[[1]], upper_95 = q[[3]])
})
print(phi_out)
write_sheet(phi_out, "Growth_performance_phi")

##############################################
# 8. Sensitivity to length and radius at hatch (BI only)
##############################################

# Refit the selected model for every combination of L0 and R0 in SENS_L0 x
# SENS_R0 and compare Linf, K, t0 and length at age 1. If these barely move,
# the results don't hinge on the exact literature values.
if (BACKCALC == "BI") {
  sens <- expand.grid(L0 = SENS_L0, R0 = SENS_R0) %>%
    pmap_dfr(function(L0, R0) {
      laa <- backcalc(bc, "BI", L0, R0)
      fit <- vb_model(model_specs[[best]], laa, get_starts(laa))
      if (is.null(fit)) return(tibble(L0_mm = L0, R0_um = R0, converged = FALSE))
      f <- fixef(fit)
      tibble(L0_mm = L0, R0_um = R0, converged = TRUE,
             Linf_fall = par_value(f, "Linf", "fall_run"), Linf_summer = par_value(f, "Linf", "summer_run"),
             K_fall = par_value(f, "K", "fall_run"), K_summer = par_value(f, "K", "summer_run"),
             t0_fall = par_value(f, "t0", "fall_run"), t0_summer = par_value(f, "t0", "summer_run"),
             L1_mean_fall   = mean(laa$length_at_age[laa$age_at_length == 1 & laa$true_eco == "fall_run"]),
             L1_mean_summer = mean(laa$length_at_age[laa$age_at_length == 1 & laa$true_eco == "summer_run"]))
    }) %>%
    mutate(delta_Linf = Linf_summer - Linf_fall)
  print(sens, width = Inf)
  write_sheet(sens, "BI_sensitivity")
}
