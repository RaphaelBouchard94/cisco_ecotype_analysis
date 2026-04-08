
#-----------------------#
#------- LIBRARY -------#
#-----------------------#

library(tidyverse)
library(magrittr)
library(nnet)
library(caret)
library(reshape2)
library(lmtest)
library(car)
library(effects)
library(jtools)
library(performance)
library("rnaturalearth")
library("rnaturalearthdata")
library("maps")
library(qrmtools)
library(sp)
library(sf)
library(paletteer)
library(plyr)
library(grid)
library(mapdata)
library(mapplots)
library(data.table)
library(sfheaders)
library(scatterpie)
library(shapefiles)
library(scales)
library(corrplot)

#--------------------#
#------- DATA -------#
#--------------------#

#---> Import assignment results from rubias

assignment_res <- read.table("pop_assign_rubias.txt")

assignment_res %<>% filter(PofZ > 0.8,z_score > -2.32)


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

assignment_res %<>% filter(community != "Whapmagoostui")

assignment_res$collection <- factor(assignment_res$collection , 
                                    levels=c("LAG_fall", "SAB_fall","South_fall","NOT_summer"),
                                    labels = c("Fall-run La Grande R", "Fall-run Sabascunica R","Fall-run South","Summer-run"))

#-----------------------#
#Let's plot the respective contribution of each stock 


assignment_res %>% 
  ggplot(aes(x=as.factor(year),fill=collection))+
  geom_bar(position = "fill",color = "black")+
  geom_text(
    aes(label = stat(count)),
    stat = "count",
    position = position_fill(vjust = 0.3),
    vjust = -0.1,
    size = 6)+
  scale_fill_manual(values = c("#0072B2", "#56B4E9" ,"#D55E00"))+
  facet_grid(~community, scale = "free_x")+
  xlab("Year")+
  ylab("Proportion")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.text =  element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18,angle = 45,vjust = 0.6),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

ggsave("figures/prop_harvest_per_community_stacked_bar.png", dpi=600,width = 10,height = 8)



%>% 
  group_by(community) %>% 
  count(collection) %>% 
  mutate(freq = n / sum(n) * 100) %>% 
  select(-n)

#------------------------------------#
#--- PROPORTION HARVEST PER MONTH ---#
#------------------------------------#

metadata <- read.csv(file = "../02_genotyping/GTscore_results/cisco_to_sequence.csv")

metadata %<>% distinct() %>% 
  dplyr::select(fishes_id,date) %>% 
  mutate(indiv = str_replace_all(fishes_id,"-","_")) %>% 
  dplyr::select(-fishes_id)

assignment_res_date <- left_join(assignment_res,metadata, c("indiv"))

assignment_res_date$date <- as.Date(assignment_res_date$date,format = "%d-%m-%Y")

# Extract the month from the date
assignment_res_date$day_month <- as.Date(format(assignment_res_date$date, "%d-%m"), format = "%d-%m")

# Plot the proportion of samples by month for each community
ggplot(assignment_res_date, aes(x = day_month, fill = collection)) +
  geom_bar(position = "stack",width = 8) +
  scale_fill_manual(values = c("#0072B2", "#56B4E9" ,"#D55E00"))+
  labs(x = "Month", y = "Proportion of Samples") +
  facet_wrap(~community, scale = "free_y")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18,angle=90),
        axis.text.y = element_text(size = 18),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_blank())

ggsave("figures/prop_catch_per_community.jpeg")

#---------------------------------------#
#--- TEST THE EFFECT OF DATE AND LAT ---#
#--- ON PROBABILITY TO HARVEST STOCK ---#
#---------------------------------------#

library(lubridate)

# Calculate week since the beginning of the year
assignment_res_date <- assignment_res_date %>%
  mutate(week_of_year = week(date))  # Gets the week number as a numeric variable

sort(unique(assignment_res_date$week_of_year))

#Week 19 removed because outlier with few data points

assignment_res_date %>% group_by(collection) %>% dplyr::count()

assignment_res_date %>% group_by(community) %>% dplyr::count()

#Chisasibi and LAG_fall and SAB_fall removed because of to few data

assignment_res_date_stats <- assignment_res_date %>% filter(week_of_year > 20, 
                                                            community != "Chisasibi",
                                                            collection != "LAG_fall") 



assignment_res_date_stats$collection <- factor(assignment_res_date_stats$collection , 
                                               levels=c("South_fall","NOT_summer"))


assignment_trial <- assignment_res_date_stats %>% group_by(week_of_year,community) %>% count()

colnames(assignment_trial) <- c("week_of_year","community","total_n")

assignment_success <- assignment_res_date_stats %>% 
  group_by(week_of_year, community) %>% 
  count(collection) %>% filter(collection == "NOT_summer")

binom_assign <- left_join(assignment_trial, assignment_success, by = c("week_of_year","community")) %>% select(-collection)

binom_assign[is.na(binom_assign)] <- 0

binom_assign %<>% mutate(fail = total_n - n) %>% ungroup()

trials = cbind(binom_assign$n,binom_assign$fail)

model.log = glm(trials ~ week_of_year + community ,
                data = binom_assign,
                family = binomial(link="logit"))

summary(model.log)


# Fit a null model (intercept-only model) for comparison
null.model = glm(trials ~ 1 ,
                 data = binom_assign,
                 family = binomial(link="logit"))


lrtest(null.model,model.log)

# Calculate McFadden's pseudo R2
with(summary(model.log), 1 - deviance/null.deviance)

#McFadden’s R-squared = 32%

plot(residuals(model.log , type = "deviance"))

# Perform Hosmer-Lemeshow test

performance_hosmer(model.log, n_bins = 6)

#---------------#
#Model selection

model.log2 = glm(trials ~ week_of_year ,
                 data = binom_assign,
                 family = binomial(link="logit"))

summary(model.log2)

model.log3 = glm(trials ~ community ,
                 data = binom_assign,
                 family = binomial(link="logit"))

summary(model.log3)

AIC(null_model,model.log,model.log2,model.log3)

#Model with both variable is better

#--------------------------#
#Let's make some prediction

data_sim <- data.frame(week_of_year = rep(seq(25,38, length = 20),3),
                       community = rep(levels(binom_assign$community),each = 20))

# Get logit predictions of whitefish
logit.predictions <- predict(object = model.log, newdata = data_sim , se=T)

# Apply inverse logit to transform to probabilities
prob.predictions <- plogis(logit.predictions$fit)

LL <- plogis(logit.predictions$fit - 1.96*logit.predictions$se.fit)
UL <- plogis(logit.predictions$fit + 1.96*logit.predictions$se.fit)


data_sim$pred <- prob.predictions
data_sim$LL <- LL
data_sim$UL <- UL

data_sim$day_month <- as.Date(format(data_sim$week_of_year, "%w"), format = "%d-%m")

ggplot(data_sim)+
  geom_line(aes(x = week_of_year, y = pred, color = community))+
  geom_ribbon(aes(x = week_of_year, ymin = LL, ymax = UL, fill = community), alpha = 0.2)+ 
  xlab("Week of the year")+
  ylab("Probability of catching summer-run")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = c(0.8,0.9),
        legend.background = element_blank())

ggsave("figures/prob_catch.png", width = 6, height = 6)


#---------------------------#
#---> Plot pie chart on map

counts <- as.data.table(assignment_res)[, .N, by = .(latitude,longitude, collection)]

test <- counts %>% group_by(latitude,longitude) %>% pivot_wider(names_from = collection, values_from=N) %>% ungroup()

test %<>% mutate(region = as.factor(group_indices(., latitude,longitude)))

test <- as.data.frame(test)

test[is.na(test)] <- 0

test %<>% relocate(longitude,.before = latitude) %>% relocate(region,.before = longitude)

colnames(test) <- c("region", "longitude", "latitude","FALLS","FALLN","SUM")

tot<- test %>% dplyr::select(FALLS:SUM) %>% mutate(tot_ind = rowSums(.)) %>% pull(tot_ind)

test$tot_ind <- tot

test_grouped <- test %>% arrange(latitude) %>%
  group_by(diff = cumsum(c(1,diff(latitude)) >= 0.1) ) %>%
  dplyr::summarise(ID = paste0(region, collapse = "/"),
                   lat = mean(latitude),
                   long = mean(longitude),
                   FALLS = sum(FALLS),
                   FALLN = sum(FALLN),
                   SUM = sum(SUM),
                   tot_ind = sum(tot_ind)) %>%
  ungroup() %>% dplyr::select(-diff)

test_grouped$radius <- scales::rescale(test_grouped$tot_ind, to = c(0.05,0.15))

lag <- read.shapefile("~/Documents/00_current_projects/Baie_James_Paper/01_population_genomics/comparative_analysis/data/lagrande/lagrande")
sab <- read.shapefile("~/Documents/00_current_projects/Baie_James_Paper/01_population_genomics/comparative_analysis/data/sabacunica/sabacunica")
sabl <- read.shapefile("~/Documents/00_current_projects/Baie_James_Paper/01_population_genomics/comparative_analysis/data/sabacunica_l/sabacunica_linea")
sabl <- read.shapefile("~/Documents/00_current_projects/Baie_James_Paper/01_population_genomics/comparative_analysis/data/sabacunica_l/sabacunica_linea")
eas_p <- read.shapefile("~/Documents/00_current_projects/Baie_James_Paper/01_population_genomics/comparative_analysis/data/eastmain/eastmain")
rup <- read.shapefile("~/Documents/00_current_projects/Baie_James_Paper/01_population_genomics/comparative_analysis/data/rupert/rupert")
not <- read.shapefile("~/Documents/00_current_projects/Baie_James_Paper/01_population_genomics/comparative_analysis/data/nottaway/nottaway")


pdf("figures/map_mixed_stock.pdf",width = 5,height = 6.5)

basemap(xlim=c(-81, -76), ylim=c(51, 55), xlab = "Longitude",bg=NA)
map('worldHires', xlim = c(-81, -76), ylim = c(51, 54.8), col="gray90",fill = T,add = T)
draw.shape(lag,type = "l",col="#0072B2")
draw.shape(sab,type = "l",col="#0072B2")
draw.shape(sabl,type = "l",col="#0072B2")
draw.shape(eas_p,type = "l",col="#56B4E9")
draw.shape(rup,type = "l",col="#56B4E9")
draw.shape(not,type = "l",col="#D55E00")
for (i in 1:nrow(test_grouped)){
  add.pie(x=test_grouped$long[i],y=test_grouped$lat[i],z=test_grouped %>% dplyr::select(FALLS:SUM) %>% .[i,] %>% unlist(.) %>% unname(.),
          radius=test_grouped$radius[i],
          labels="",
          col = c(
            "#56B4E9", #EAS-C
            "#0072B2", #BIA-L
            "#D55E00" #RUPERT)
          ))
}
mapplots::legend.bubble(x= -80.5, y = 54.7,z=c(100,400),maxradius = 0.2,n=2)

dev.off()

