rm(list = ls())

####################################
###########Library##################
####################################

library(tidyverse)
library(magrittr)
library(data.table)
library(gplots)


####################################
######  Admixture analysis  ########
####################################

ld_all <- fread("data/chr34/all_maf0.01pctind0.8_maxdepth10_CM078251.1_by_50.ld")

mat_ld_all<-matrix(ncol=4)

colnames(mat_ld_all)<-c("pos1","pos2","var","r2")

head(ld_all)

colnames(ld_all)<-c("pos1","pos2","var","r2")

ld_all <- ld_all[-which(ld_all$pos1==ld_all$pos2),] # remove pairwise SNPs in which the same SNP is involved

variable<-levels(ld_all$var)

mat_ld_all<-rbind(mat_ld_all,ld_all)

mat_ld_all<-mat_ld_all[2:dim(mat_ld_all)[1],]

head(mat_ld_all)

##############################################################################
# ld aa karyotype only:

ld_aa <- fread("data/chr34/karyoaa_CM078251_by_50.ld")

mat_ld_aa <-matrix(ncol=4)

colnames(mat_ld_aa)<-c("pos1","pos2","var","r2")

colnames(ld_aa)<-c("pos1","pos2","var","r2")

ld_aa <- ld_aa[-which(ld_aa$pos1==ld_aa$pos2),] # remove pairwise SNPs in which the same SNP is involved

variable_aa<-levels(ld_aa$var)

mat_ld_aa<-rbind(mat_ld_aa,ld_aa)

mat_ld_aa <-mat_ld_aa[2:dim(mat_ld_aa)[1],]

targ_win <- c(34800000,35619999)
targ_win_plus <- c((34800000-1000000)/1000,(35619999+1000000)/1000)


# Filter both datasets to the target window
mat_ld_all_filt <- mat_ld_all %>% 
  filter(pos1 > targ_win_plus[1] & pos1 < targ_win_plus[2] & 
           pos2 > targ_win_plus[1] & pos2 < targ_win_plus[2])

mat_ld_aa_filt <- mat_ld_aa %>% 
  filter(pos1 > targ_win_plus[1] & pos1 < targ_win_plus[2] & 
           pos2 > targ_win_plus[1] & pos2 < targ_win_plus[2])

# Create upper triangle data (mat_ld_all) - where pos2 > pos1
upper_triangle <- mat_ld_all_filt %>%
  filter(pos2 > pos1) %>%
  mutate(dataset = "All")

# Create lower triangle data (mat_ld_summer) - SWAP pos1 and pos2
lower_triangle <- mat_ld_aa_filt %>%
  mutate(pos1_new = pos2,
         pos2_new = pos1) %>%
  select(pos1 = pos1_new, pos2 = pos2_new, var, r2) %>%
  mutate(dataset = "Summer")

# Combine the two triangles
combined_data <- bind_rows(upper_triangle, lower_triangle)

# Plot
ggplot(combined_data, aes(x=pos1/1000, y=pos2/1000)) +
  geom_tile(aes(fill=r2)) +
  geom_rect(inherit.aes = FALSE, 
            aes(xmin = 34.800, xmax = 35.619,
                ymin = 34.800, ymax = 35.619),
            color = "red",
            fill = "transparent",
            linewidth = 1) +
  scale_fill_viridis_c() +
  xlab("Chromosome 36 (Mbp)") + 
  ylab("Chromosome 36 (Mbp)") +
  theme_bw(base_size = 20) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  coord_fixed(ratio = 1)

ggsave(filename = "figures/05_karyoaa_all_zoom_LDmat_50000snp_combined.png", 
       width = 15, height = 12, unit = 'cm', dpi = 900)


##################################################################
# karyo bb only

ld_bb <- fread("data/chr34/karyobb_CM078251.1_by_50")

mat_ld_bb<-matrix(ncol=4)

colnames(mat_ld_bb)<-c("pos1","pos2","var","r2")

colnames(ld_bb)<-c("pos1","pos2","var","r2")

ld_bb <- ld_bb[-which(ld_bb$pos1==ld_bb$pos2),] # remove pairwise SNPs in which the same SNP is involved

variable_sum<-levels(ld_bb$var)

mat_ld_bb <-rbind(mat_ld_bb,ld_bb)

mat_ld_bb <-mat_ld_bb[2:dim(mat_ld_bb)[1],]

targ_win <- c(34800000,35619999)
targ_win_plus <- c((34800000-1000000)/1000,(35619999+1000000)/1000)


# Filter both datasets to the target window
mat_ld_all_filt <- mat_ld_all %>% 
  filter(pos1 > targ_win_plus[1] & pos1 < targ_win_plus[2] & 
           pos2 > targ_win_plus[1] & pos2 < targ_win_plus[2])

mat_ld_bb_filt <- mat_ld_bb %>% 
  filter(pos1 > targ_win_plus[1] & pos1 < targ_win_plus[2] & 
           pos2 > targ_win_plus[1] & pos2 < targ_win_plus[2])

# Create upper triangle data (mat_ld_all) - where pos2 > pos1
upper_triangle <- mat_ld_all_filt %>%
  filter(pos2 > pos1) %>%
  mutate(dataset = "All")

# Create lower triangle data (mat_ld_summer) - SWAP pos1 and pos2
lower_triangle2 <- mat_ld_bb_filt %>%
  mutate(pos1_new = pos2,
         pos2_new = pos1) %>%
  select(pos1 = pos1_new, pos2 = pos2_new, var, r2) %>%
  mutate(dataset = "bb")

# Combine the two triangles
combined_data_bb <- bind_rows(upper_triangle, lower_triangle2)

# Plot
ggplot(combined_data_bb, aes(x=pos1/1000, y=pos2/1000)) +
  geom_tile(aes(fill=r2)) +
  geom_rect(inherit.aes = FALSE,
            aes(xmin = 34800/1000, xmax = 35619/1000,
                ymin = 34800/1000, ymax = 35619/1000),
            color = "red",
            fill = "transparent",
            linewidth = 1) +
  scale_fill_viridis_c("LD") +
  xlab("Chromosome 34 (Mbp)") + 
  ylab("Chromosome 34 (Mbp)") +
  theme_bw(base_size = 20) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  coord_fixed(ratio = 1)

ggsave(filename = "figures/05_karyo_bb_all_zoom_LDmat_50000snp_combined.png", 
       width = 15, height = 12, unit = 'cm', dpi = 900)




