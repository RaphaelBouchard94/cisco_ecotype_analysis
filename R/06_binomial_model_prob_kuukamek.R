
#-----------------------#
#------- LIBRARY -------#
#-----------------------#

library(tidyverse)
library(lme4)
library(DHARMa)
library(performance)
library(car)
library(lattice)
library(ggplot2)
library(splines)
library(ggeffects)


#--------------------#
#------- DATA -------#
#--------------------#

#-Format data
#---> Import assignment results from rubias

assignment_res <- read.table("data/migrating_individuals_assignment_rubias.txt")

#Remove individuals with low assignment confidence and from potential unsampled
#populations

assignment_res %<>% filter(PofZ > 0.99)

#---> Mixed-stock proportion per fishing location

assignment_res %>% dplyr::group_by(latitude,longitude) %>% dplyr::count() 

#42 different fishing LOCATIONS

#---> Format data for plotting
assignment_res %<>% mutate(year = case_when(str_detect(indiv,"_20") ~ 2020,
                                            str_detect(indiv,"_21") ~ 2021,
                                            str_detect(indiv,"_22") ~ 2022,
                                            str_detect(indiv,"_23") ~ 2023,
                                            str_detect(indiv,"NA") ~ 2021))

assignment_res$community <- factor(assignment_res$community , 
                                   levels=c("Whapmagoostui", "Chisasibi", "Wemindji", "Eastmain","Waskaganish"))

# assignment_res %<>% filter(community != "Whapmagoostui")

assignment_res$collection <- factor(assignment_res$collection , 
                                    levels=c("V1", "V2","V3"),
                                    labels = c("Fall-run north","Fall-run south","Kuukamek"))


metadata <- read.csv(file = "data/migrating_samples_metadata.csv")

metadata %<>% distinct() %>% 
  dplyr::select(fishes_id,date) %>% 
  mutate(indiv = str_replace_all(fishes_id,"-","_")) %>% 
  dplyr::select(-fishes_id)

assignment_res_date <- left_join(assignment_res,metadata, c("indiv"))

assignment_res_date$date <- as.Date(assignment_res_date$date,format = "%d-%m-%Y")

# Extract the month from the date
assignment_res_date$day_month <- as.Date(format(assignment_res_date$date, "%d-%m"), format = "%d-%m")

#---------------------------------------#
#--- TEST THE EFFECT OF DATE AND LAT ---#
#--- ON PROBABILITY TO HARVEST STOCK ---#
#---------------------------------------#

# Calculate week since the beginning of the year
assignment_res_date <- assignment_res_date %>%
  mutate(week_of_year = week(date))  # Gets the week number as a numeric variable

sort(unique(assignment_res_date$week_of_year))

#Week 19 removed because outlier with few data points

assignment_res_date %>% group_by(collection) %>% dplyr::count()

assignment_res_date %>% group_by(community) %>% dplyr::count()

#Remove month with two few data points
assignment_res_date_stats <- assignment_res_date %>% filter(week_of_year > 25) 

assignment_res_date_stats <- assignment_res_date_stats %>%
  mutate(
    is_kuukamek = ifelse(collection == "Kuukamek", 1, 0),
    week_c = scale(week_of_year, center = TRUE, scale = FALSE)
  )

#===============================================
# COMPREHENSIVE COMPARISON OF MODEL OPTIONS
#===============================================

# List to store models
models <- list()
diagnostics <- list()

# Model 0: Null
models$null <- glmer(
  is_kuukamek ~ 1 + (1 | community),
  data = assignment_res_date_stats,
  family = binomial
)


# Model 1: Linear
models$linear <- glmer(
  is_kuukamek ~ week_c + (1 | community),
  data = assignment_res_date_stats,
  family = binomial
)

# Model 2: Quadratic
models$quadratic <- glmer(
  is_kuukamek ~ week_c + week_c^2 + (1 | community),
  data = assignment_res_date_stats,
  family = binomial
)

# Model 3: Natural spline
models$spline <- glmer(
  is_kuukamek ~ ns(week_c, df = 3) + (1 | community),
  data = assignment_res_date_stats,
  family = binomial
)

# Run diagnostics for each
for(model_name in names(models)) {
  cat("\n=== Testing:", model_name, "===\n")
  
  sim <- simulateResiduals(models[[model_name]], n = 1000, plot = FALSE)
  
  # Store results
  diagnostics[[model_name]] <- list(
    ks_test = testUniformity(sim, plot = FALSE),
    disp_test = testDispersion(sim, plot = FALSE),
    outlier_test = testOutliers(sim, plot = FALSE)
  )
  
  # Print results
  cat("KS test p-value:", diagnostics[[model_name]]$ks_test$p.value, "\n")
  cat("Dispersion p-value:", diagnostics[[model_name]]$disp_test$p.value, "\n")
  cat("Outlier p-value:", diagnostics[[model_name]]$outlier_test$p.value, "\n")
}

# Compare AICs
aic_comparison <- AIC(models$null,models$linear, models$quadratic, models$spline)
aic_comparison$model <- rownames(aic_comparison)
aic_comparison <- aic_comparison %>% arrange(AIC)

print(aic_comparison)

# Choose best model
best_model_name <- aic_comparison$model[1]
cat("\nBest model by AIC:", best_model_name, "\n")


#####
#Put results in supp table

wb <- loadWorkbook("res/supplementary_table.xlsx")

addWorksheet(wb, "binomial_model_mixed_stock")

writeData(
  wb,
  sheet = "binomial_model_mixed_stock",
  aic_comparison
)

saveWorkbook(wb, "res/supplementary_table.xlsx", overwrite = TRUE)


##########################################
#Check if other df parameters fit better than df=3

# Compare different df values using AIC
df_comparison <- data.frame(
  df = integer(),
  AIC = numeric(),
  BIC = numeric()
)

for(df_val in 2:5) {
  mod_temp <- glmer(
    is_kuukamek ~ ns(week_c, df = df_val) + (1 | community),
    data = assignment_res_date_stats,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa")
  )
  
  df_comparison <- rbind(df_comparison, 
                         data.frame(
                           df = df_val,
                           AIC = AIC(mod_temp),
                           BIC = BIC(mod_temp)
                         ))
}

df_comparison$delta_AIC <- df_comparison$AIC - min(df_comparison$AIC)
df_comparison$delta_BIC <- df_comparison$BIC - min(df_comparison$BIC)

print(df_comparison)

# Lowest AIC is not df=3 but AIC close to 4 and BIC is lower for df = 3
#Will go with most parcimonious df = 3

#########################################
#Spline is the best model

mod_final <- glmer(
  is_kuukamek ~ ns(week_c, df = 3) + (1 | community),
  data = assignment_res_date_stats,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

summary(mod_final)

r2_final <- r2_nakagawa(mod_final)

# DHARMa diagnostics
sim_resid_final <- simulateResiduals(mod_final, n = 1000, plot = TRUE)

# All tests
testUniformity(sim_resid_final)    # Should be better now!
testDispersion(sim_resid_final)    # Check overdispersion
testOutliers(sim_resid_final)      # Check outliers
testZeroInflation(sim_resid_final) # Check zero-inflation

# Check residuals vs predictors
plotResiduals(sim_resid_final, form = assignment_res_date_stats$week_c)
plotResiduals(sim_resid_final, form = assignment_res_date_stats$community)

# Performance metrics
check_model(mod_final)
model_performance(mod_final)

# Overdispersion
check_overdispersion(mod_final)

# R-squared
r2_final <- r2_nakagawa(mod_final)
cat("\nMarginal R² (fixed effects):", round(r2_final$R2_marginal, 3))
cat("\nConditional R² (fixed + random):", round(r2_final$R2_conditional, 3), "\n")

# Random effects
icc(mod_final)
dotplot(ranef(mod_final, condVar = TRUE))

#-------------------------------------------------------#
# Get predictions across the range of weeks
pred_final <- ggpredict(mod_final, terms = "week_c [all]")

# Back-transform to original week scale
mean_week <- mean(assignment_res_date_stats$week_of_year)
pred_df <- as.data.frame(pred_final)
pred_df$week <- pred_df$x + mean_week


# Calculate observed proportions for comparison
week_summary <- assignment_res_date_stats %>%
  dplyr::group_by(week_of_year) %>%
  dplyr::summarise(
    prop_kuukamek = mean(is_kuukamek),
    n = n(),
    se = sqrt(prop_kuukamek * (1 - prop_kuukamek) / n),
    .groups = "drop"
  )


#################################################################

ref_year <- 2021

pred_df <- as.data.frame(pred_final)

# Back-transform week
pred_df$week <- pred_df$x + mean_week

# Convert week -> date
pred_df$date <- as.Date(
  paste(ref_year, round(pred_df$week), 2, sep = "-"),
  format = "%Y-%U-%u"
)


week_summary <- assignment_res_date_stats %>%
  dplyr::group_by(week_of_year) %>%
  dplyr::summarise(
    prop_kuukamek = mean(is_kuukamek),
    n = n(),
    se = sqrt(prop_kuukamek * (1 - prop_kuukamek) / n),
    .groups = "drop"
  ) %>%
  mutate(
    date = as.Date(
      paste(ref_year, week_of_year, 1, sep = "-"),
      format = "%Y-%U-%u"
    )
  )

p1 <- ggplot() +
  geom_point(
    data = week_summary,
    aes(x = date, y = prop_kuukamek, size = n),
    alpha = 0.5, color = "gray30"
  ) +
  geom_line(
    data = pred_df,
    aes(x = date, y = predicted),
    color = "#E69F00", linewidth = 1.3
  ) +
  geom_ribbon(
    data = pred_df,
    aes(x = date, ymin = conf.low, ymax = conf.high),
    fill = "#E69F00", alpha = 0.2
  ) +
  scale_size_continuous(range = c(2, 6), name = "Sample size") +
  scale_x_date(
    date_labels = "%B",
    date_breaks = "1 month",
    expand = c(0.02, 0.02)
  ) +
  labs(
    x = "Month",
    y = "Probability of catching Kuukamek"
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill = NA),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black"),
    axis.text.x = element_text(hjust = 0.5)
  )

print(p1)


ggsave("figures/06_kuukamek_probability_by_week.png", p1, 
       width = 5, height = 5, dpi = 600)


##########################################################
#Some information based on model

# Get week range and convert to dates
week_range <- range(assignment_res_date_stats$week_of_year)
year <- 2024  # Adjust if needed

date_start <- ymd(paste(year, "01", "01")) + weeks(week_range[1] - 1)
date_end <- ymd(paste(year, "01", "01")) + weeks(week_range[2] - 1)


# Find peak week
pred_final <- ggpredict(mod_final, terms = "week_c [all]")
pred_df <- as.data.frame(pred_final)
pred_df$week <- pred_df$x + mean(assignment_res_date_stats$week_of_year)
pred_df$date <- ymd(paste(year, "01", "01")) + weeks(round(pred_df$week) - 1)

peak_week <- pred_df$week[which.max(pred_df$predicted)]
peak_date <- pred_df$date[which.max(pred_df$predicted)]
peak_prob <- max(pred_df$predicted)

min_prob <- min(pred_df$predicted)
prob_range <- peak_prob - min_prob

low_date <- pred_df$date[which.min(pred_df$predicted)]



