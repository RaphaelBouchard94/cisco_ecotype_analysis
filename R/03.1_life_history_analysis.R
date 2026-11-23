rm(list = ls())

###############################################################################
# 03.1 - Age vs fork length (Fig. 2b) and age-length keys for unaged fish
#
# CHANGES vs the previous version
#  1. Fixed the typo that stopped the script (cc.unage -> cc.unagedS).
#  2. Fig. 2b: coord_cartesian() instead of ylim(). ylim() doesn't zoom: it
#     DELETES points outside 250-425 mm before drawing the boxplots.
#  3. The age-length key now covers the full length range of the fish, so
#     alkIndivAge() can't fail on a fish outside 260-420 mm.
#  4. Removed check_model() (diagnostic plots needing extra packages) and
#     packages that weren't used.
#
# NB: the ages inferred here (infered_age) are not used anywhere in the
# manuscript. The output goes to 00_data/ as before; script 01 reads
# data/data_gillrak_geomorph_age.csv, so it isn't affected by this script.
# If any other analysis (e.g. the GWAS age covariate) uses infered_age, the
# Methods need a sentence on age-length keys.
###############################################################################

#################
#### Library ####
#################

library(FSA)       # lencat(), alkIndivAge(), alkPlot()
library(nnet)      # multinom()
library(tidyverse)
library(magrittr)

dir.create("figures", showWarnings = FALSE)
dir.create("00_data", showWarnings = FALSE)

#################
### Load Data ###
#################

data <- read.csv("data/data_with_gillrakers.csv", header = TRUE, sep = ",")
age  <- read.csv("data/data_agegrowth_cisco.csv")

########################
# Fig. 2b: fork length vs otolith age
########################

age %<>% left_join(data, by = c("Id" = "sample_id"))

fig2b <- age %>%
  filter(!is.na(age)) %>%
  ggplot(aes(x = as.factor(age), y = fork_length, color = true_eco)) +
  geom_boxplot() +
  geom_jitter(position = position_jitterdodge(0.3), alpha = 0.5) +
  ylab("Fork Length (mm)") +
  xlab("Age") +
  coord_cartesian(ylim = c(250, 425)) +          ### CHANGED (was ylim())
  scale_color_manual(labels = c("Fall run", "Summer run"), values = c("skyblue", "#E69F00")) +
  theme_bw(base_size = 16) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(),
        legend.position = "none")
fig2b
ggsave("figures/03_age_fl.png", fig2b, width = 10, height = 10, unit = "cm", dpi = 900)

# How many otolith-aged fish per ecotype and age (useful for the caption)
age %>% filter(!is.na(age)) %>% count(true_eco, age) %>% as_tibble() %>% print(n = Inf)

########################################
# Age-length keys (one per ecotype)
########################################

# NB: columns are picked by position as in your original script (column 2 = Id,
# columns 11-21 = age and otolith columns). Check this if the CSV changes.
data %<>% left_join(age[, c(2, 11:21)], by = c("sample_id" = "Id"))

### CHANGED: length categories covering every fish (10-mm bins)
lens <- seq(floor(min(data$long_fourche, na.rm = TRUE) / 10) * 10,
            floor(max(data$long_fourche, na.rm = TRUE) / 10) * 10, by = 10)

make_alk <- function(aged) {
  aged %<>% mutate(lcat10 = lencat(long_fourche, w = 10))
  print(xtabs(~ lcat10 + age, data = aged))
  mlr <- multinom(age ~ lcat10, data = aged, maxit = 500, trace = FALSE)
  alk <- predict(mlr, data.frame(lcat10 = lens), type = "probs")
  row.names(alk) <- lens
  alk
}

### SUMMER-RUN
cc.agedS <- data %>% filter(!is.na(age), true_eco == "summer_run")
alk.sm   <- make_alk(cc.agedS)

### FALL-RUN
cc.agedF <- data %>% filter(!is.na(age), true_eco == "fall_run")
alk.f    <- make_alk(cc.agedF)

png("figures/ALK_summer.jpeg", units = "cm", width = 10, height = 10, res = 200)
alkPlot(alk.sm, type = "area", pal = "gray", showLegend = TRUE, leg.cex = 0.7, xlab = "Fork Length (mm)")
dev.off()

png("figures/ALK_fall.jpeg", units = "cm", width = 10, height = 10, res = 200)
alkPlot(alk.f, type = "area", pal = "gray", showLegend = TRUE, leg.cex = 0.7, xlab = "Fork Length (mm)")
dev.off()

#############################
## Infer age of unaged fish
#############################

cc.unagedS <- data %>% filter(is.na(age), !is.na(long_fourche), true_eco == "summer_run")   ### CHANGED (typo)
cc.unagedS.mod <- alkIndivAge(alk.sm, age ~ long_fourche, data = cc.unagedS)

cc.unagedF <- data %>% filter(is.na(age), !is.na(long_fourche), true_eco == "fall_run")
cc.unagedF.mod <- alkIndivAge(alk.f, age ~ long_fourche, data = cc.unagedF)

cc.fnl <- bind_rows(cc.unagedS.mod, cc.agedS,
                    cc.unagedF.mod, cc.agedF)

ggplot(cc.fnl, aes(x = as.factor(age), y = long_fourche, color = true_eco)) +
  geom_boxplot() +
  geom_jitter(position = position_jitterdodge(0.3)) +
  ylab("Fork length (mm)") + xlab("Age") +
  scale_color_manual(labels = c("Fall run", "Summer run"), values = c("#56B4E9", "#D55E00")) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        legend.position = c(0.2, 0.9))
ggsave("figures/Age-length_relationship.jpeg", width = 6, height = 6)

ggplot(cc.fnl, aes(x = age, fill = true_eco)) +
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single")) +
  ylab("Count") + xlab("Age") +
  scale_fill_manual(labels = c("Fall run", "Summer run"), values = c("#56B4E9", "#D55E00")) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        legend.position = c(0.8, 0.8))
ggsave("figures/Age_distribution.jpeg", width = 6, height = 6)

#################################################
# Effect of age and ecotype on length (all fish, incl. inferred ages)
#################################################

summary(lm(long_fourche ~ age, data = cc.fnl))
summary(lm(long_fourche ~ age + true_eco, data = cc.fnl))

data %<>% left_join(cc.fnl[, c("sample_id", "age")], by = "sample_id")
data %<>% dplyr::rename(infered_age = age.y, true_age = age.x)

###############################################
# Export whole dataset with age data
###############################################

write.csv(data, "00_data/data_gillrak_geomorph_age.csv", row.names = FALSE, quote = FALSE)
