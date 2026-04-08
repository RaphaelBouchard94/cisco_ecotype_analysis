rm(list = ls())

#################
#### Library ####
#################

library(tidyverse)
library(magrittr)
library(performance)
library(FSA)
library(FSAsim)
library(FSAmisc)
library(nlstools)
library(lme4)
library(rstatix)
library(ggpubr)
library(arm)
library(grid)
library(patchwork)
library(nnet)


#################
### Load Data ###
#################

data <- read.csv("data/data_with_gillrakers.csv", header = T,sep = ",")

age <- read.csv("data/data_agegrowth_cisco.csv")

########################
#Get length at age data for regression

age %<>% left_join(data,by=c("Id"="sample_id"))


#Visualize length vs age

age %>% filter(!is.na(age)) %>% ggplot(aes(x=as.factor(age),y = fork_length,color=true_eco))+
  geom_boxplot()+
  geom_jitter(position = position_jitterdodge(0.3), alpha = 0.5)+
  ylab("Fork Length (mm)")+
  xlab("Age")+
  ylim(c(250,425))+
  scale_color_manual(labels=c("Fall run","Summer run"),values = c("skyblue","#E69F00"))+
  theme_bw(base_size = 16)+
  theme(panel.border = element_rect(colour = "black",fill=NA),
        panel.grid = element_blank(),
        legend.position = "none")

ggsave("figures/03_age_fl.png",
       width = 10, height = 10, unit = "cm", dpi = 900)


########################################
#Comparison between Summer and Fall run#
########################################

data %<>% left_join(age[,c(2,11:21)],by=c("sample_id"="Id"))

#Remove individuals without age data
cc.agedS <- data %>% filter(!is.na(age),true_eco=="summer_run")

cc.agedS %<>% mutate(lcat10 = lencat(long_fourche, w=10))

xtabs(~lcat10+age,data=cc.agedS)

cc.mlrS <- multinom(age ~ lcat10, data = cc.agedS, maxit = 500)

lens <- seq(260,420,10)

alk.sm <- predict(cc.mlrS,data.frame(lcat10=lens),type="probs")
row.names(alk.sm) <- lens

###FALL-RUN
cc.agedF <- data %>% filter(!is.na(age),true_eco=="fall_run")

cc.agedF %<>% mutate(lcat10 = lencat(long_fourche, w=10))

xtabs(~lcat10+age,data=cc.agedF)

cc.mlrF <- multinom(age ~ lcat10, data = cc.agedF, maxit = 500)

lens <- seq(260,420,10)

alk.f <- predict(cc.mlrF,data.frame(lcat10=lens),type="probs")
row.names(alk.f) <- lens

png("figures/ALK_summer.jpeg",units="cm",width=10,height=10,res=200)
alkPlot(alk.sm,type="area",pal="gray",showLegend=T,leg.cex = 0.7, xlab = "Fork Length (mm)")
dev.off()

png("figures/ALK_fall.jpeg",units="cm",width=10,height=10,res=200)
alkPlot(alk.f,type="area",pal="gray",showLegend=T,leg.cex = 0.7, xlab = "Fork Length (mm)")
dev.off()

#############################
##Infer age for each ecotype

#Summer-run

cc.unage <- data %>% filter(is.na(age),true_eco == "summer_run") 
cc.unagedS.mod <- alkIndivAge(alk.sm, age ~ long_fourche,data=cc.unagedS)

#Fall-run

cc.unagedF <- data %>% filter(is.na(age),true_eco == "fall_run") 
cc.unagedF.mod <- alkIndivAge(alk.f, age ~ long_fourche,data=cc.unagedF)

cc.fnl <- rbind(cc.unagedS.mod,cc.agedS %>% dplyr::select(-lcat10),
                cc.unagedF.mod,cc.agedF %>% dplyr::select(-lcat10))


#Age-length relationship

ggplot(cc.fnl, aes(x=as.factor(age),y=long_fourche,color=true_eco))+
  geom_boxplot()+
  geom_jitter(position = position_jitterdodge(0.3))+
  ylab("Fork length (mm)")+
  xlab("Age")+
  scale_color_manual(labels=c("Fall run","Summer run"),values = c("#56B4E9","#D55E00"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        legend.position = c(0.2,0.9))

ggsave("figures/Age-length_relationship.jpeg",width=6,height=6)

#Age distribution
ggplot(cc.fnl,aes(x=age,fill=true_eco))+
  geom_bar(position=position_dodge2(width = 0.9, preserve = "single"))+
  ylab("Count")+
  xlab("Age")+
  scale_fill_manual(labels=c("Fall run","Summer run"),values = c("#56B4E9","#D55E00"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        legend.position = c(0.8,0.8))

ggsave("figures/Age_distribution.jpeg",width=6,height=6)
#################################################
#Significant effect of age AND ecotype on length

lm1 <- lm(long_fourche~age,data=cc.fnl)
summary(lm1)
check_model(lm1)

lm2 <- lm(long_fourche ~ age + true_eco,data=cc.fnl)
summary(lm2)

#tvalue of ecotype = 14.246
#tvalue of age = 6.158

str(cc.fnl)

data %<>% left_join(cc.fnl[,c("sample_id","age")],by="sample_id") 

data %<>% dplyr::rename(infered_age = age.y,true_age = age.x)

###############################################
#Export whole dataset with age data

write.csv(data,"00_data/data_gillrak_geomorph_age.csv",row.names = F,quote=F)

