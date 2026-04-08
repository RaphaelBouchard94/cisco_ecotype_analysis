rm(list = ls())

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

source(file = "ProjectOrthogonal.R")
source(file = "TestOfAngle.R")
source(file = "critical_angle.R")

rad2deg <- function(rad) {(rad * 180) / (pi)}

#################
### Load Data ###
#################
data_tps <- readland.tps("data/cisco_rup_bay.TPS", specID = "ID")

data_gpa <- gpagen(data_tps)

plot(data_gpa)

PCA <- gm.prcomp(data_gpa$coords)

##Remove body arching from fish 

data_tps_matrix_tmp <- apply(data_gpa$coords,2,c)

N_Landmark=16
N_ind=190

result_matrix <-list()
test_vite_fait_per_ind =split(data_tps_matrix_tmp,sort(rep(1:N_ind,N_Landmark)))

for(i in 1:N_ind){
  
  result_matrix[[i]] <- test_vite_fait_per_ind[[i]] %>% 
    t() %>% 
    as.data.frame() %>%
    unlist
  
}

result_matrix <- do.call(rbind,result_matrix)

colnames(result_matrix) 

v <- c(1:N_Landmark)
colnames(result_matrix) <- c(paste("X",v,sep = ""),paste("Y",v,sep = ""))

order_col <- paste(rep(c("X","Y"),N_Landmark),sort(rep(1:N_Landmark,2)), sep = "")

result_matrix <- result_matrix[,order_col]

#Make data corrected projection
data_tps_proj <- ProjectOrthogonal(result_matrix,PCA$rotation[,1])

data_tps_proj

result_matrix_reverse <- list()

for(i in 1:nrow(data_tps_proj)){
  result_matrix_reverse [[i]] <- data_tps_proj[i,]%>% matrix(data = .,ncol = 2, byrow = T)
}

result_matrix_reverse <- as.data.frame(do.call(rbind,result_matrix_reverse))

lvl <- (nrow(result_matrix_reverse)/N_Landmark)

data_gpa_proj <- array(unlist(by(result_matrix_reverse, 
                                 INDICES = rep(1:lvl, each = N_Landmark), 
                                 FUN = function(x) x)),
                       dim = c(N_Landmark, ncol(result_matrix_reverse), lvl))

write_rds(data_gpa_proj,"tmp/data_gpa_proj.Rds")

#Compute PCA with Barnaby corrected data
PCA_proj <- gm.prcomp(data_gpa_proj)

write_rds(PCA_proj,"tmp/PCA_proj.Rds")

#Test if corrected PC1 is significantly similar to PC1 uncorrected
TestOfAngle(PCA$rotation[,1],PCA_proj$rotation[,1], flip=F)
TestOfAngle(PCA$rotation[,1],PCA_proj$rotation[,2], flip=F)

final_data_gpa_projected <- gpagen(data_gpa_proj)

write_rds(final_data_gpa_projected,"tmp/final_data_gpa_projected.Rds")






