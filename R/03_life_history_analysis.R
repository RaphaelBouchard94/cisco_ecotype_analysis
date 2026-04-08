rm(list = ls())

########################
# 1 - Library  #
########################

library(tidyverse)
library(magrittr)
library(performance)
library(FSA)
library(FSAsim)
library(FSAmisc)
library(nlstools)
library(nlme)
library(lme4)
library(MASS)
library(ggpubr)
library(grid)
library(patchwork)
library(scales)
library(openxlsx)


########################
# 2 - Load Data  #
########################

data <- read.table("data/data_gillrak_geomorph_age_morphres.csv", header = T)

data$ecotype <- as.factor(data$true_eco)

str(data)

data %>% filter(!is.na(true_age)) %>% ggplot(aes(x=as.factor(true_age),fill = ecotype))+
  geom_bar(position=position_dodge2(width = 0.9, preserve = "single"))+
  xlab("Age")+
  ylab("Count")+
  scale_fill_manual(labels=c("Fall run","Summer run"),values = c("skyblue1","#E69F00"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        legend.position = c(0.8,0.8))

# ggsave("figures/age_dist.jpeg",width = 10,height=10)  




data %>% dplyr::group_by(ecotype) %>% dplyr::summarise(mean_eco = mean(true_age,na.rm=T),
                                                       sd_eco = sd(true_age,na.rm=T))

data$rescale_age <- rescale(data$true_age)

t.test(rescale_age~ecotype,data=data)

#

data$rescale_fl <- rescale(data$long_fourche)

size_diff <- lm(rescale_fl ~ rescale_age,data=data)

summary(size_diff)

size_diff_eco <- lm(rescale_fl ~ rescale_age*ecotype,data=data)

summary(size_diff_eco)

#
data %>% group_by(is.na(true_age),ecotype) %>% count()


##########################
##Back-calculated length##
##########################

## 1- Calculate total radius size

back_calc_tmp <- mutate_all(data %>% dplyr::select(X1_um:X9_um), ~replace_na(.,0))

back_calc <- bind_cols(data %>% dplyr::select(-(X1_um:X9_um)),back_calc_tmp)


back_calc %<>% group_by(sample_id) %>% dplyr::mutate(total_rad = max(across(X1_um:X9_um)))

#Remove individuals without data, total_rad = 0

back_calc %<>% filter(total_rad > 0)

back_calc %>%
  ggplot(aes(x=long_fourche,y=total_rad,color=ecotype))+
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("Fork Length (mm)")+
  ylab("Otholith radius")+
  scale_color_manual(values = c("skyblue","#E69F00"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18))


##2 - Regress length against radius

##without ecotype

lm1 <- lm(log(long_fourche)~log(total_rad), data=back_calc)

summary(lm1)

##with ecotype

lm2 <- lm(log(long_fourche)~log(total_rad)+ecotype, data=back_calc)

summary(lm2)

#Model assumptions are respected, normality of residuals not perfect but still ok

##with ecotype + interaction

lm3 <- lm(log(long_fourche)~log(total_rad)*ecotype, data=back_calc)

summary(lm3)

#Pas d'interaction!

##3 - Back calculate fish length using regression results

#Using model lm2 slope 

lm2 <- lm(log(long_fourche)~log(total_rad)+ecotype, data=back_calc)

summary(lm2)

v = lm2$coefficients[2]

#v = 

#Calculate length at age using equation
##Li=[(Si/S)]^v*L

back_mutate <- back_calc %>% mutate(across(X1_um:X9_um,
                             ~ (.x / total_rad)^v * long_fourche,
                             .names = "age{col}"))

##Change 0 to NA
length_at_age <- back_mutate %>% dplyr::select(ageX1_um:ageX9_um) %>% dplyr::mutate_all(~na_if(., 0))

##Join length at age data to the dataset
length_at_age %<>% left_join(back_calc %>% dplyr::select(sample_id:total_rad),by = "sample_id")

##Tidy dataframe 
length_at_age %<>% pivot_longer(ageX1_um:ageX9_um,names_to = "back_age",values_to = "length_at_age")

##Create column "back age" with age as numeric

length_at_age %<>% mutate(age_at_length = case_when(back_age == "ageX1_um" ~ 1,
                                                    back_age == "ageX2_um" ~ 2,
                                                    back_age == "ageX3_um" ~ 3,
                                                    back_age == "ageX4_um" ~ 4,
                                                    back_age == "ageX5_um" ~ 5,
                                                    back_age == "ageX6_um" ~ 6,
                                                    back_age == "ageX7_um" ~ 7,
                                                    back_age == "ageX8_um" ~ 8,
                                                    back_age == "ageX9_um" ~ 9))

str(length_at_age)

##Plot length at age relationship

length_at_age %>% 
  ggplot(aes(x=as.factor(age_at_length),y=length_at_age,color=ecotype,group=sample_id))+
  geom_line()+
  geom_point(alpha=0.4)+
  xlab("Age")+
  ylab("Fork Length (mm)")

###########################
#Between-group comparaison#
###########################

# Implementing model with random effect

#create river of origin column

length_at_age %<>% mutate(river = if_else(str_detect(sample_id, "NOT"),"NOT","RUP"))

length_at_age$true_eco <- as.factor(length_at_age$true_eco)

#list random effects

random = list(
  river = pdDiag(Linf ~ 1),
  sample_id = pdDiag(Linf ~ 1)
)

#Provide starting value for VB model

sv0 <- vbStarts(length_at_age ~ age_at_length, data = length_at_age)
sv0

#

vb <- function(age, Linf, K, t0) {
  Linf * (1 - exp(-K * (age - t0)))
}


#Full model
fit_mix_Linf_river <- nlme(
  length_at_age ~ vb(age_at_length, Linf, K, t0),
  data = length_at_age,
  
  fixed = Linf + K + t0 ~ true_eco,
  random = Linf ~ 1 | river/sample_id,
  
  start = list(fixed = c(
    Linf = c(sv0["Linf"], 0),  # intercept and slope for true_eco
    K    = c(sv0["K"], 0),
    t0   = c(sv0["t0"], 0)
  )),
  
  na.action = na.omit
)

(fit_mix_Linf_river)

#Rivers do not meaningfully differ in asymptotic length once ecotype is accounted for.
#Most unexplained variation is among individuals

#Linf fall-run : 311.64
#Linf summer-run: 360.34

#No evidence for ecotype differences in growth rate.
#No ecotype effect on t0.

####
#Model selection

#model 1 - no ecotype effect
start0 <- c(
  Linf = as.numeric(sv0["Linf"]),
  K    = as.numeric(sv0["K"]),
  t0   = as.numeric(sv0["t0"])
)

m0 <- nlme(
  length_at_age ~ vb(age_at_length, Linf, K, t0),
  data = length_at_age,
  fixed = Linf + K + t0 ~ 1,
  random = Linf ~ 1 | river/sample_id,
  start = start0,
  na.action = na.omit
)

#model2: Ecoype affects Linf only
start_m2 <- c(
  "Linf.(Intercept)" = start0["Linf"],
  "Linf.true_ecosummer_run" = 0,   # difference from baseline
  "K.(Intercept)" = start0["K"],
  "t0.(Intercept)" = start0["t0"]
)

m2 <- nlme(
  length_at_age ~ vb(age_at_length, Linf, K, t0),
  data = length_at_age,
  fixed = list(
    Linf ~ true_eco,  # Linf varies by ecotype
    K    ~ 1,
    t0   ~ 1
  ),
  random = Linf ~ 1 | river/sample_id,
  start = start_m2,
  na.action = na.omit
)

#model 3 : Ecoype affects K only
start_m3 <- c(
  "Linf.(Intercept)" = start0["Linf"],
  "K.(Intercept)" = start0["K"],
  "K.true_ecosummer_run" = 0,   # difference from baseline
  "t0.(Intercept)" = start0["t0"]
)

m3 <- nlme(
  length_at_age ~ vb(age_at_length, Linf, K, t0),
  data = length_at_age,
  fixed = list(
    Linf ~ 1,  # Linf varies by ecotype
    K    ~ true_eco,
    t0   ~ 1
  ),
  random = Linf ~ 1 | river/sample_id,
  start = start_m3,
  na.action = na.omit
)

#model 4 : Ecoype affects t0 only
start_m4 <- c(
  "Linf.(Intercept)" = start0["Linf"],
  "K.(Intercept)" = start0["K"],
  "t0.(Intercept)" = start0["t0"],
  "t0.true_ecosummer_run" = 0  # difference from baseline
)

m4 <- nlme(
  length_at_age ~ vb(age_at_length, Linf, K, t0),
  data = length_at_age,
  fixed = list(
    Linf ~ 1,  # Linf varies by ecotype
    K    ~ 1,
    t0   ~ true_eco
  ),
  random = Linf ~ 1 | river/sample_id,
  start = start_m4,
  na.action = na.omit
)

#model 5: ecotype affects Linf + K

start_m5 <- c(
  "Linf.(Intercept)" = start0["Linf"],
  "Linf.true_ecosummer_run" = 0,
  "K.(Intercept)" = start0["K"],
  "K.true_ecosummer_run" = 0,
  "t0.(Intercept)" = start0["t0"]
)

m5 <- nlme(
  length_at_age ~ vb(age_at_length, Linf, K, t0),
  data = length_at_age,
  fixed = list(
    Linf ~ true_eco,
    K    ~ true_eco,  # K now varies by ecotype
    t0   ~ 1
  ),
  random = Linf ~ 1 | river/sample_id,
  start = start_m5,
  na.action = na.omit
)

#model 6: ecotype affects Linf + t0

start_m6 <- c(
  "Linf.(Intercept)" = start0["Linf"],
  "Linf.true_ecosummer_run" = 0,
  "K.(Intercept)" = start0["K"],
  "t0.(Intercept)" = start0["t0"],
  "t0.true_ecosummer_run" = 0
)

m6 <- nlme(
  length_at_age ~ vb(age_at_length, Linf, K, t0),
  data = length_at_age,
  fixed = list(
    Linf ~ true_eco,
    K    ~ 1,  # K now varies by ecotype
    t0   ~ true_eco
    ),
  random = Linf ~ 1 | river/sample_id,
  start = start_m6,
  na.action = na.omit
)

#model 7: ecotype affects K + t0

start_m7 <- c(
  "Linf.(Intercept)" = start0["Linf"],
  "K.(Intercept)" = start0["K"],
  "K.true_ecosummer_run" = 0,
  "t0.(Intercept)" = start0["t0"],
  "t0.true_ecosummer_run" = 0
)

m7 <- nlme(
  length_at_age ~ vb(age_at_length, Linf, K, t0),
  data = length_at_age,
  fixed = list(
    Linf ~ 1,
    K    ~ true_eco,  # K now varies by ecotype
    t0   ~ true_eco
  ),
  random = Linf ~ 1 | river/sample_id,
  start = start_m7,
  na.action = na.omit
)

#model 8 : Ecotype affects Linf+K+t0

start_m8 <- c(
  "Linf.(Intercept)" = start0["Linf"],
  "Linf.true_ecosummer_run" = 0,
  "K.(Intercept)" = start0["K"],
  "K.true_ecosummer_run" = 0,
  "t0.(Intercept)" = start0["t0"],
  "t0.true_ecosummer_run" = 0
)

m8 <- nlme(
  length_at_age ~ vb(age_at_length, Linf, K, t0),
  data = length_at_age,
  fixed = list(
    Linf ~ true_eco,
    K    ~ true_eco,
    t0   ~ true_eco   #
  ),
  random = Linf ~ 1 | river/sample_id,
  start = start_m8,
  na.action = na.omit
)

#AIC comparison

models <- list(
  m0 = m0,
  mL = m2,
  mK = m3,
  mt0 = m4,
  mLK = m5,
  mLt0 = m6,
  mKt0 = m7,
  mLKt0 = m8
)

model_summary <- do.call(rbind, lapply(names(models), function(mname) {
  mod <- models[[mname]]
  data.frame(
    Model = mname,
    AIC   = AIC(mod),
    BIC   = BIC(mod),
    logLik = logLik(mod)
  )
}))

model_summary

#Save in Supplementary table:

wb <- loadWorkbook("res/supplementary_table.xlsx")

addWorksheet(wb, "Model_selection_AIC")
writeData(wb, sheet = "Model_selection_AIC", model_summary)

saveWorkbook(wb, "res/supplementary_table.xlsx", overwrite = TRUE)


#Save best model

# Extract fixed effects
fixef_tab <- summary(models$mL)$tTable %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Parameter") %>%
  dplyr::rename(
    Estimate = Value,
    SE = Std.Error,
    DF = DF,
    t_value = `t-value`,
    p_value = `p-value`
  )

fixef_tab %<>% mutate(across(.cols = Estimate:p_value, ~ round(.x , 4)))

addWorksheet(wb, "Best_model_parameters")

writeData(wb, sheet = "Best_model_parameters", fixef_tab )

# Save workbook
saveWorkbook(wb, "res/supplementary_table.xlsx", overwrite = TRUE)

#########################################
#Model-based predictions based on mL

fe <- fixed.effects(m2)
fe

# Population-level parameter table
pars <- data.frame(
  ecotype = c("fall_run", "summer_run"),
  Linf = c(
    fe["Linf.(Intercept)"],
    fe["Linf.(Intercept)"] + fe["Linf.true_ecosummer_run"]
  ),
  K = rep(fe["K"], 2),   # K constant
  t0 = rep(fe["t0"], 2)  # t0 constant
)

pars

# Predict population-level growth curves
age_seq <- seq(0, 9, length.out = 200)

pred_pop <- expand.grid(
  age = age_seq,
  ecotype = pars$ecotype
) %>%
  merge(pars, by = "ecotype") %>%
  mutate(length = Linf * (1 - exp(-K * (age - t0))))

# Confidence intervals via parametric bootstrap
vc <- vcov(m2)
nsim <- 1000
beta_sim <- mvrnorm(nsim, mu = fe, Sigma = vc)

get_curve <- function(beta_row, eco, age) {
  # Linf differs by ecotype; K and t0 are constant
  if (eco == "fall_run") {
    Linf <- beta_row["Linf.(Intercept)"]
  } else {
    Linf <- beta_row["Linf.(Intercept)"] + beta_row["Linf.true_ecosummer_run"]
  }
  K  <- beta_row["K"]
  t0 <- beta_row["t0"]
  
  Linf * (1 - exp(-K * (age - t0)))
}

ci_pop <- expand.grid(
  age = age_seq,
  ecotype = pars$ecotype
) %>%
  rowwise() %>%
  mutate(
    sim = list(apply(beta_sim, 1, get_curve, eco = ecotype, age = age)),
    LCI = quantile(unlist(sim), 0.025, na.rm = T),
    UCI = quantile(unlist(sim), 0.975, na.rm = T)
  ) %>%
  ungroup()

# Plot per population growth rates
ggplot() +
  geom_point(data = length_at_age,
             aes(age_at_length, length_at_age, color = true_eco),
             alpha = 0.3) +
  geom_line(data = pred_pop,
            aes(age, length, color = ecotype),
            linewidth = 1.2) +
  geom_ribbon(data = ci_pop,
              aes(age, ymin = LCI, ymax = UCI, fill = ecotype),
              alpha = 0.25, color = NA) +
  scale_color_manual(values = c("skyblue","#E69F00"))+
  scale_fill_manual(values = c("skyblue","#E69F00"))+
  scale_x_continuous(breaks = seq(0, 9, by = 1)) +
  theme_bw(base_size = 16)+
  theme(panel.border = element_rect(colour = "black",fill=NA),
        panel.grid = element_blank(),
        legend.position = "none")+
  labs(x = "Age", y = "Fork length (mm)")

ggsave("figures/03_vb_growth_rates.png",
       width = 10, height = 10, unit = "cm", dpi = 900)


##############################################
#Test for difference in omega

#Compute for each ecotype

pars$omega <- with(pars, log10(K) + 2 * log10(Linf))
pars

omega_sim <- data.frame(
  fall = log10(
    beta_sim[, "K"]
  ) + 2 * log10(
    beta_sim[, "Linf.(Intercept)"]
  ),
  summer = log10(
    beta_sim[, "K"] + beta_sim[, "K"]
  ) + 2 * log10(
    beta_sim[, "Linf.(Intercept)"] + beta_sim[, "Linf.true_ecosummer_run"]
  )
)

apply(omega_sim, 2, quantile, probs = c(0.025, 0.5, 0.975))

#Test delta mega via parametric bootstrap

delta_omega <- with(
  omega_sim,
  summer - fall
)

quantile(delta_omega, probs = c(0.025, 0.5, 0.975))

#There is a significant difference in growth performance among ecotypes


