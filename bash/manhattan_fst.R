setwd("~/2024-11-01_cisco_james_bay/01_scripts/")

#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(magrittr)
library(data.table)
library(patchwork)

#------------#
#--- DATA ---#
#------------#

#Sliding window

NOT_fall_sum <- fread("../12_fst/NOT_fall_NOT_summer.slidingwindow")[,-1]

colnames(NOT_fall_sum) <- c("chr","pos","nsites","fst")

RUP_NOT_summer <- fread("../12_fst/RUP_NOT_summer.slidingwindow")[,-1]

colnames(RUP_NOT_summer) <- c("chr","pos","nsites","fst")

EAS_NOT_summer <- fread("../12_fst/EAS_NOT_summer.slidingwindow")[,-1]

colnames(EAS_NOT_summer) <- c("chr","pos","nsites","fst")

RUP_NOT_fall <- fread("../12_fst/RUP_NOT_fall.slidingwindow")[,-1]

colnames(RUP_NOT_fall) <- c("chr","pos","nsites","fst")

NOT_sum_LAG <- fread("../12_fst/LAG_NOT_summer.slidingwindow")[,-1]

colnames(NOT_sum_LAG) <- c("chr","pos","nsites","fst")

NOT_fall_LAG <- fread("../12_fst/LAG_NOT_fall.slidingwindow")[,-1]

colnames(NOT_fall_LAG) <- c("chr","pos","nsites","fst")

#----------------#
#--- Analysis ---#
#----------------#

#For plotting purpose:

axis_man <- NOT_fall_sum %>% 
  group_by(chr) %>%
  summarize(center = mean(pos))

#NOT fall - NOT summer

a <- ggplot(NOT_fall_sum, aes(x = pos , y = fst, color = as_factor(chr))) +
  geom_point(alpha = 0.2, size=0.5) +
  geom_smooth(color = "black",linewidth = .5)+
  scale_x_continuous(
    label = axis_man$chr,
    breaks = axis_man$center,
    expand = c(0.05,0)) +
  scale_color_manual(values = rep(
    c("grey70", "grey90"),
    unique(length(axis_man$chr)))) +
  facet_grid(~chr, scales="free_x") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  labs(
    x = "NOT_summer - NOT_fall",
    y = "Fst") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0,'lines'),
        axis.title.y =  element_text(angle = 90, size = 12, vjust = 0.5),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90))

#NOT summer - RUP

b <- ggplot(RUP_NOT_summer, aes(x = pos , y = fst, color = as_factor(chr))) +
  geom_point(alpha = 0.2, size=0.5) +
  geom_smooth(color = "black",linewidth = .5)+
  scale_x_continuous(
    label = axis_man$chr,
    breaks = axis_man$center,
    expand = c(0.05,0)) +
  scale_color_manual(values = rep(
    c("grey70", "grey90"),
    unique(length(axis_man$chr)))) +
  facet_grid(~chr, scales="free_x") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  labs(
    x = "NOT_summer - RUP",
    y = "Fst") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0,'lines'),
        axis.title.y =  element_text(angle = 90, size = 12, vjust = 0.5),
        axis.text.x = element_blank(),
        strip.text.x = element_blank())

#NOT summer - EAS

c <- ggplot(EAS_NOT_summer, aes(x = pos , y = fst, color = as_factor(chr))) +
  geom_point(alpha = 0.2, size=0.5) +
  geom_smooth(color = "black",linewidth = .5)+
  scale_x_continuous(
    label = axis_man$chr,
    breaks = axis_man$center,
    expand = c(0.05,0)) +
  scale_color_manual(values = rep(
    c("grey70", "grey90"),
    unique(length(axis_man$chr)))) +
  facet_grid(~chr, scales="free_x") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  labs(
    x = "NOT_summer - EAS",
    y = "Fst") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0,'lines'),
        axis.title.y =  element_text(angle = 90, size = 12, vjust = 0.5),
        axis.text.x = element_blank(),
        strip.text.x = element_blank())


#Plot with NOT summer vs South fall populations

NOT_SOUTH <- a/b/c

NOT_SOUTH

ggsave("../12.1_fst_manhattan/NOT_summer_SOUTH.jpeg", plot = NOT_SOUTH)

#NOT fall - NOT summer

d <- ggplot(NOT_fall_sum, aes(x = pos , y = fst, color = as_factor(chr))) +
  geom_point(alpha = 0.2, size=0.5) +
  geom_smooth(color = "black",linewidth = .5)+
  scale_x_continuous(
    label = axis_man$chr,
    breaks = axis_man$center,
    expand = c(0.05,0)) +
  scale_color_manual(values = rep(
    c("grey70", "grey90"),
    unique(length(axis_man$chr)))) +
  facet_grid(~chr, scales="free_x") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  labs(
    x = "NOT_summer - NOT_South",
    y = "Fst") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0,'lines'),
        axis.title.y =  element_text(angle = 90, size = 12, vjust = 0.5),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90))


#LAG - NOT summer

e <- ggplot(NOT_sum_LAG, aes(x = pos , y = fst, color = as_factor(chr))) +
  geom_point(alpha = 0.2, size=0.5) +
  geom_smooth(color = "black",linewidth = .5)+
  scale_x_continuous(
    label = axis_man$chr,
    breaks = axis_man$center,
    expand = c(0.05,0)) +
  scale_color_manual(values = rep(
    c("grey70", "grey90"),
    unique(length(axis_man$chr)))) +
  facet_grid(~chr, scales="free_x") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  labs(
    x = "NOT_summer - LAG",
    y = "Fst") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0,'lines'),
        axis.title.y =  element_text(angle = 90, size = 12, vjust = 0.5),
        axis.text.x = element_blank(),
        strip.text.x = element_blank())

#LAG - NOT fall

f <- ggplot(NOT_fall_LAG, aes(x = pos , y = fst, color = as_factor(chr))) +
  geom_point(alpha = 0.2, size=0.5) +
  geom_smooth(color = "black",linewidth = .5)+
  scale_x_continuous(
    label = axis_man$chr,
    breaks = axis_man$center,
    expand = c(0.05,0)) +
  scale_color_manual(values = rep(
    c("grey70", "grey90"),
    unique(length(axis_man$chr)))) +
  facet_grid(~chr, scales="free_x") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  labs(
    x = "NOT_fall - LAG",
    y = "Fst") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0,'lines'),
        axis.title.y =  element_text(angle = 90, size = 12, vjust = 0.5),
        axis.text.x = element_blank(),
        strip.text.x = element_blank())

#RUP - NOT fall

g <- ggplot(RUP_NOT_fall, aes(x = pos , y = fst, color = as_factor(chr))) +
  geom_point(alpha = 0.2, size=0.5) +
  geom_smooth(color = "black",linewidth = .5)+
  scale_x_continuous(
    label = axis_man$chr,
    breaks = axis_man$center,
    expand = c(0.05,0)) +
  scale_color_manual(values = rep(
    c("grey70", "grey90"),
    unique(length(axis_man$chr)))) +
  facet_grid(~chr, scales="free_x") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  labs(
    x = "NOT_fall - RUP",
    y = "Fst") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0,'lines'),
        axis.title.y =  element_text(angle = 90, size = 12, vjust = 0.5),
        axis.text.x = element_blank(),
        strip.text.x = element_blank())

#Plot with NOT summer vs other populations

NOT_POP <- d/e/f/g

ggsave("../12.1_fst_manhattan/NOT_fall_LAG.jpeg", plot = NOT_POP)

