
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
library(lme4)
library(performance)
library(lubridate)

#--------------------#
#------- DATA -------#
#--------------------#

#---> Import assignment results from rubias

assignment_res <- read.table("data/migrating_individuals_assignment_rubias.txt")

#Remove individuals with low assignment confidence and from potential unsampled
#populations

assignment_res %<>% filter(PofZ > 0.99)

#---> Mixed-stock proportion per fishing location

assignment_res %>% dplyr::group_by(latitude,longitude) %>% dplyr::count() 

#correct weird latitude
assignment_res %<>% mutate(longitude = if_else(longitude < -100, -78.9,longitude))


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
  scale_fill_manual(values = c("#0072B2","skyblue1","#E69F00"))+
  facet_grid(~community, scale = "free_x")+
  xlab("Year")+
  ylab("Proportion")+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 0.6),
        legend.title = element_blank(),
        legend.position = "none")

ggsave("figures/06_prop_harvest_per_community_stacked_bar.png", dpi=600,width = 14,height = 8)


#mean proportion of harvest of each stock per community

assignment_res %>% 
  group_by(community,year) %>% 
  dplyr::count(collection) %>% 
  pivot_wider(values_from = n, names_from = collection) %>% 
  replace(is.na(.),0) %>% 
  mutate(prop_kuuk = Kuukamek/(`Fall-run north` + `Fall-run south`+ Kuukamek)) %>% 
  ungroup() %>% 
  group_by(community) %>% 
  dplyr::summarise(mean_prop_kuuk = mean(prop_kuuk))
  

#------------------------------------#
#--- PROPORTION HARVEST PER MONTH ---#
#------------------------------------#

metadata <- read.csv(file = "data/migrating_samples_metadata.csv")

metadata %<>% distinct() %>% 
  dplyr::select(fishes_id,date) %>% 
  mutate(indiv = str_replace_all(fishes_id,"-","_")) %>% 
  dplyr::select(-fishes_id)

assignment_res_date <- left_join(assignment_res,metadata, c("indiv"))

assignment_res_date$date <- as.Date(assignment_res_date$date,format = "%d-%m-%Y")

# Extract the month from the date
assignment_res_date$day_month <- as.Date(format(assignment_res_date$date, "%d-%m"), format = "%d-%m")

# Plot the proportion of samples by month for each community
ggplot(assignment_res_date %>% filter(week_of_year > 25), aes(x = day_month, fill = collection)) +
  geom_bar(position = "stack",width = 7) +
  scale_fill_manual(values = c("#0072B2","skyblue1","#E69F00"))+
  labs(x = "Month", y = "Subsistence fisheries samples") +
  theme_bw(base_size = 20)+
  theme(legend.position = "none",)

ggsave("figures/06_prop_catch_per_community.png", dpi=600,width = 5,height = 5)

#---------------------------#
#---> Plot pie chart on map

counts <- as.data.table(assignment_res)[, .N, by = .(latitude,longitude, collection)]

test <- counts %>% group_by(latitude,longitude) %>% pivot_wider(names_from = collection, values_from=N) %>% ungroup()

test %<>% mutate(region = as.factor(group_indices(., latitude,longitude)))

test <- as.data.frame(test)

test[is.na(test)] <- 0

test %<>% relocate(longitude,.before = latitude) %>% relocate(region,.before = longitude)

colnames(test) <- c("region", "longitude", "latitude","FALLS","FALLN","KUU")

tot<- test %>% dplyr::select(FALLS:KUU) %>% mutate(tot_ind = rowSums(.)) %>% pull(tot_ind)

test$tot_ind <- tot

test_grouped <- test %>% arrange(latitude) %>%
  group_by(diff = cumsum(c(1,diff(latitude)) >= 0.01) ) %>%
  dplyr::summarise(ID = paste0(region, collapse = "/"),
                   lat = mean(latitude),
                   long = mean(longitude),
                   FALLS = sum(FALLS),
                   FALLN = sum(FALLN),
                   SUM = sum(KUU),
                   tot_ind = sum(tot_ind)) %>%
  ungroup() %>% dplyr::select(-diff)

test_grouped$radius <- scales::rescale(test_grouped$tot_ind, to = c(0.1,0.2))

lag <- read.shapefile("data/map/lagrande/lagrande")
sab <- read.shapefile("data/map/sabacunica/sabacunica")
sabl <- read.shapefile("data/map/sabacunica_l/sabacunica_linea")
sabl <- read.shapefile("data/map/sabacunica_l/sabacunica_linea")
eas_p <- read.shapefile("data/map/eastmain/eastmain")
rup <- read.shapefile("data/map/rupert/rupert")
not <- read.shapefile("data/map/nottaway/nottaway")

#######################

# Improved overlap resolution with better convergence
resolve_overlaps <- function(coords, radii, max_iter = 200, repel_force = 1.5) {
  n <- nrow(coords)
  adjusted <- coords
  original <- coords  # Keep track of originals
  
  for (iter in 1:max_iter) {
    forces_x <- rep(0, n)
    forces_y <- rep(0, n)
    
    # Calculate repulsion forces between all pairs
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        dx <- adjusted$long[j] - adjusted$long[i]
        dy <- adjusted$lat[j] - adjusted$lat[i]
        dist <- sqrt(dx^2 + dy^2)
        
        # Minimum distance (sum of radii + buffer)
        min_dist <- (radii[i] + radii[j]) * repel_force + 0.02
        
        if (dist < min_dist && dist > 0.0001) {
          # Calculate repulsion force (stronger when closer)
          force <- (min_dist - dist) / dist
          
          # Apply forces to both points
          forces_x[i] <- forces_x[i] - force * dx
          forces_y[i] <- forces_y[i] - force * dy
          forces_x[j] <- forces_x[j] + force * dx
          forces_y[j] <- forces_y[j] + force * dy
        }
      }
    }
    
    # Add spring force back to original positions (keeps pies near true location)
    spring_strength <- 0.1
    for (i in 1:n) {
      forces_x[i] <- forces_x[i] + spring_strength * (original$long[i] - adjusted$long[i])
      forces_y[i] <- forces_y[i] + spring_strength * (original$lat[i] - adjusted$lat[i])
    }
    
    # Apply forces with dampening
    dampening <- 0.5
    adjusted$long <- adjusted$long + forces_x * dampening
    adjusted$lat <- adjusted$lat + forces_y * dampening
    
    # Check for convergence
    max_force <- max(sqrt(forces_x^2 + forces_y^2))
    if (max_force < 0.001) {
      message(paste("Converged after", iter, "iterations"))
      break
    }
  }
  
  # Add column showing how far each pie moved
  adjusted$moved <- sqrt((adjusted$long - original$long)^2 + 
                           (adjusted$lat - original$lat)^2)
  
  return(adjusted)
}

# Apply with adjusted parameters
coords_adjusted <- resolve_overlaps(
  data.frame(long = test_grouped$long, lat = test_grouped$lat),
  radii = test_grouped$radius,
  max_iter = 300,
  repel_force = 1.5  # Increase this if still overlapping (try 1.8 or 2.0)
)

# Complete plotting code with leader lines
png("figures/06_map_mixed_stock.png", 
    width = 5, height = 6.5, 
    units = "in",           # Specify units as inches
    res = 600,              # 600 DPI for publication quality
    pointsize = 8,          # Smaller point size for crisper text
    bg = "white")   
basemap(xlim=c(-81, -76), ylim=c(51, 57), xlab = "Longitude", bg=NA)
map('worldHires', xlim = c(-81, -76), ylim = c(51, 57), col="gray90", fill = T, add = T)

# Draw river lines
draw.shape(lag, type = "l", col="#0072B2")
draw.shape(sab, type = "l", col="#0072B2")
draw.shape(sabl, type = "l", col="#0072B2")
draw.shape(eas_p, type = "l", col="skyblue1")
draw.shape(rup, type = "l", col="skyblue1")
draw.shape(not, type = "l", col="#E69F00")

# Draw leader lines first (so they appear behind pies)
for (i in 1:nrow(test_grouped)){
  if (coords_adjusted$moved[i] > 0.02) {  # Only draw if moved significantly
    segments(x0 = test_grouped$long[i], 
             y0 = test_grouped$lat[i],
             x1 = coords_adjusted$long[i], 
             y1 = coords_adjusted$lat[i],
             col = "gray40", lty = 2, lwd = 0.8)
    # Original position marker
    points(test_grouped$long[i], test_grouped$lat[i], 
           pch = 20, cex = 0.5, col = "black")
  }
}

# Draw pies
for (i in 1:nrow(test_grouped)){
  add.pie(x=coords_adjusted$long[i], 
          y=coords_adjusted$lat[i],
          z=test_grouped %>% dplyr::select(FALLS:SUM) %>% .[i,] %>% unlist(.) %>% unname(.),
          radius=test_grouped$radius[i],
          labels="",
          col = c("#56B4E9", "#0072B2", "#E69F00"))
}

mapplots::legend.bubble(x= -81, y = 56.5, z=c(50,130), maxradius = 0.2, n=2)
maps::map.scale(x = -77.8, y = 51.2, ratio = FALSE, relwidth = 0.3, cex = 0.8)

dev.off()



##############

png("figures/06_map_pop_gen.png", 
    width = 2.4, 
    height = 3,  # keep aspect ratio
    units = "in",
    res = 600,
    pointsize = 6.5,  # will print at 6.5 pt
    bg = "white")    
basemap(xlim=c(-81, -76), ylim=c(51, 57), xlab = "Longitude", bg=NA)
map('worldHires', xlim = c(-81, -76), ylim = c(51, 57), col="gray90", fill = T, add = T)

# Draw river lines
draw.shape(lag, type = "l", col="black")
draw.shape(sab, type = "l", col="black")
draw.shape(sabl, type = "l", col="black")
draw.shape(eas_p, type = "l", col="black")
draw.shape(rup, type = "l", col="black")
draw.shape(not, type = "l", col="black")

maps::map.scale(x = -77.8, y = 51.2, ratio = FALSE, relwidth = 0.3, cex = 0.8)

dev.off()
 

