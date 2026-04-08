rm(list = ls())

#########################
####### LIBRARY #########
#########################

library(tidyverse)
library(magrittr)
library(rubias)
library(paletteer)
library(openxlsx)

tmp <-  read.table("data/gtscore_analysis/polyGenResults_singleSNP_rubias.txt",header=T, 
                   stringsAsFactors=FALSE, 
                   colClasses = c("character"))

#Are there different spelling of sample type?
tmp %>% group_by(sample_type) %>% count()

tmp %>% filter(is.na(sample_type)) %>% dplyr::select(indiv)

tmp$indiv <- str_replace(tmp$indiv, 
                          "([A-Za-z]*s)_(\\d+)_", 
                          function(x) {
                            prefix <- str_extract(x, "[A-Za-z]*s(?=_)")
                            num <- str_extract(x, "(?<=[A-Za-z]s_)\\d+")
                            paste0(prefix, "_", sprintf("%03d", as.numeric(num)), "_")
                          })

tmp %<>% mutate(sample_type = if_else(is.na(sample_type),"reference",sample_type))

#-----------------------------------------#
####### ASSESS BASELINE ACCURACY #########
#-----------------------------------------#

#First let's divide our dataset into source, hold-out and mixed samples
#Format source pop file

tmp_cisco_source <- tmp %>% filter(sample_type == "reference")

source <- read.table("data/ngsadmix/admixture_assign_k3_allpop.txt",header=T)

source %<>% dplyr::select(individual_id,genetic_cluster)

colnames(source) <- c("indiv","pop_clus")

source %<>% mutate(indiv = str_replace(indiv,"-","_"))

tmp_cisco_source %<>% 
  left_join(source, by = c("indiv")) %>% 
  dplyr::select(-collection,repunit) %>% 
  dplyr::rename(collection = pop_clus) %>% 
  mutate(repunit = collection)

cisco_source <- tmp_cisco_source[,c(1,2,824,3,4:823)]

gtdev_samples <- read.table("data/id_gtseq_dev.tsv")

colnames(gtdev_samples) <- "indiv"

gtdev_samples %<>% mutate(indiv = str_replace(indiv,"-","_"))

##################################
#Make training dataset

train <- cisco_source %>% filter(indiv %in% gtdev_samples$indiv) 

train %>% group_by(collection) %>% dplyr::count()

#filter out V1 and NA in collection
train %<>% filter(!is.na(collection))

train %>% group_by(collection) %>% dplyr::count()

##############################
#Make holdout dataset
holdout <- cisco_source %>% filter(!indiv %in% gtdev_samples$indiv & !is.na(collection)) 

holdout %>% group_by(collection) %>% dplyr::count()

#remove populations from holdout for which we did not have samples for in the
#training dataset

holdout %<>% filter(collection %in% train$collection)

holdout %<>% distinct()

#-------------------------------------------------------#
####### TEST ASSIGNMENT POWER USING HOLDOUT SET #########
#-------------------------------------------------------#


assign_holdout <- infer_mixture(reference = train, 
                         mixture = holdout, 
                         gen_start_col = 5,
                         reps = 20000,
                         burn_in = 5000,
                         pb_iter = 10000,
                         method = "MCMC")


assign_holdout_map <- assign_holdout$indiv_posteriors %>%
  dplyr::group_by(indiv) %>%
  top_n(1, PofZ) %>%
  ungroup()

assign_holdout_map %<>% left_join(holdout %>% dplyr::select(indiv,collection), by="indiv")

assign_holdout_map %<>% dplyr::rename("collection" = "collection.x",
                              "true_pop"= "collection.y")

assign_holdout_map %<>% dplyr::mutate(assignment_res = if_else(true_pop == collection,"1","0"))

supp_table <- assign_holdout_map %>% 
  group_by(true_pop) %>% 
  summarise(
    accuracy = sum(assignment_res == "1") / n(),
    correct = sum(assignment_res == "1"),
    total = n()
  )

#####
#Put results in supp table

wb <- loadWorkbook("res/supplementary_table.xlsx")

addWorksheet(wb, "GTseq_holdout_assignment")

writeData(
  wb,
  sheet = "GTseq_holdout_assignment",
  supp_table 
)

saveWorkbook(wb, "res/supplementary_table.xlsx", overwrite = TRUE)


######################################################
#Estimate dispersal between north south with samples taken within rivers
#during upstream migration
metadata_source <- read.csv("data/JB_samples_metadata.csv", header = T)

colnames(metadata_source) <- c("indiv","pop")

metadata_source %<>% mutate(indiv = str_replace(indiv,"-","_"))

unknown_source <- cisco_source %>% filter(is.na(collection)) 

unknown_source %>% group_by(collection) %>% dplyr::count()

unknown_source %<>% left_join(metadata_source, by = "indiv") 

unknown_source$collection <- unknown_source$pop

unknown_source %<>% dplyr::select(-pop) %>% distinct()

assign_dispersal <- infer_mixture(reference = train, 
                                  mixture = unknown_source, 
                                  gen_start_col = 5,
                                  reps = 20000,
                                  burn_in = 5000,
                                  pb_iter = 10000,
                                  method = "MCMC")


assign_disp_map <- assign_dispersal$indiv_posteriors %>%
  dplyr::group_by(indiv) %>%
  top_n(1, PofZ) %>%
  ungroup()

assign_disp_map %<>% left_join(unknown_source %>% dplyr::select(indiv,collection), by="indiv")

assign_disp_map %<>% dplyr::rename("collection" = "collection.x",
                                      "true_pop"= "collection.y")

assign_disp_map %<>% dplyr::mutate(assignment_res = if_else(true_pop == collection,"1","0"))

assign_disp_map  %>% 
  group_by(true_pop) %>% 
  dplyr::count(collection)

#-----------------------------------------#
#Get cluster assignment for all source individuals


unknown_gt <- assign_disp_map %>% 
  select(indiv,collection)

holdout_gt <- assign_holdout_map %>% 
  select(indiv, collection)

train_gt <- cisco_source %>% filter(indiv %in% gtdev_samples$indiv) %>% select(indiv,collection)

gt_pca <- rbind(train_gt,holdout_gt,unknown_gt)

write.table(gt_pca,"data/baseline_assignment_gtseq.txt",quote = F,row.names = F)

#-----------------------------------------#
####### INFER MIXTURE PROPORTION #########
#-----------------------------------------#

#Make reference dataset with training and holdout samples

ref_all <- cisco_source %>% filter(!is.na(collection)) 

#Subset mixed-stock samples

tmp_cisco_mixed <- tmp %>% filter(sample_type == "mixture")

cisco_mixed <- tmp_cisco_mixed %>% dplyr::select(starts_with("Chr"))

cisco_mixed<- cbind(tmp_cisco_mixed[,1:4],cisco_mixed)

mix_est <- infer_mixture(reference = ref_all , 
                         mixture = cisco_mixed, 
                         gen_start_col = 5,
                         reps = 20000,
                         burn_in = 5000,
                         pb_iter = 10000,
                         method = "MCMC")


map_rows <- mix_est$indiv_posteriors %>%
  dplyr::group_by(indiv) %>%
  top_n(1, PofZ) %>%
  ungroup()

normo <- tibble(z_score = rnorm(1e06))

thresh <- quantile(normo$z_score, 0.01)

ggplot(map_rows, aes(x = z_score)) +
  geom_density(colour = "blue") +
  geom_density(data = normo, colour = "black")+
  geom_vline(xintercept = -2.32,linetype = 2,color = "red")+
  theme_bw()

ggsave("figures/06_migrating_samples_zscore_distribution.png",height=10,width =10, unit = "cm", dpi = 600)

#Visual inspection of z-scores showed that the distribution of
#posterior assignment probabilities in mixed-stock samples were
#more skewed

length(which(map_rows$z_score < quantile(normo$z_score, 0.01)))

#42 samples have a z-score which is smaller than the 0.99 confidence
#interval

#-------------------------------------------------#
#--> Merge assignment results with geographic info
#-------------------------------------------------#

geo_location <- read.csv(file = "data/migrating_samples_metadata.csv")

geo_location %<>% dplyr::select(community,river.name,latitude,longitude,fishes_id) %>% 
  mutate(id = str_replace_all(fishes_id,"-","_")) %>% 
  dplyr::select(-fishes_id)

map_rows %<>% left_join(geo_location, by = c("indiv" = "id"))

map_rows %<>% dplyr::select(-missing_loci)

map_rows$community <- factor(map_rows$community , 
                             levels=c("Whapmagoostui", "Chisasibi", "Wemindji", "Eastmain","Waskaganish"))


map_rows %>% ggplot(aes(x = reorder(mixture_collection,latitude), y = z_score))+
  geom_boxplot()+
  geom_jitter(alpha = 0.1)+
  geom_hline(yintercept = unname(thresh),linetype = 2,color = "red")+
  facet_grid(~ community,scale = "free_x")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))

ggsave("figures/06_zscore_per_sampling_site.png",height=10,width =10, unit = "cm", dpi = 600)

map_rows %>% group_by(collection,community) %>% filter(z_score < unname(thresh)) %>% dplyr::count()

map_rows %>% mutate(assigned = if_else(PofZ > 0.99, "assigned","non-assigned")) %>% 
  group_by(assigned) %>% 
  dplyr::count()

#2 non-assigned based on PofZ threshold of 0.99

mean(map_rows$PofZ)

map_rows_filt <- map_rows %>% filter(z_score > unname(thresh)) 


###################

write.table(map_rows,"data/migrating_individuals_assignment_rubias.txt")

###################

supp_table_assign_res <- map_rows %>% mutate(repunit = case_when(repunit == "V1" ~ "Fall-run North",
                                        repunit == "V2" ~ "Fall-run South",
                                        repunit == "V3" ~ "Kâcikâsikumekw"),
                    collection = case_when(collection == "V1" ~ "Fall-run North",
                                        collection == "V2" ~ "Fall-run South",
                                        collection == "V3" ~ "Kâcikâsikumekw"),)

#####
#Put results in supp table

wb <- loadWorkbook("res/supplementary_table.xlsx")

addWorksheet(wb, "Mixed_stock_pop")

writeData(
  wb,
  sheet = "Mixed_stock_pop",
    supp_table_assign_res
)

saveWorkbook(wb, "res/supplementary_table.xlsx", overwrite = TRUE)


