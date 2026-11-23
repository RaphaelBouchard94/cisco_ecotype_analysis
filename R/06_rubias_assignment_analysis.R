#------------------------------------------------------------------------------#
# 06_rubias_assignment_analysis.R
# Baseline validation (training/holdout) and genetic assignment of cisco
# harvested in Cree subsistence fisheries using rubias.
#------------------------------------------------------------------------------#

rm(list = ls())

library(tidyverse)
library(magrittr)
library(rubias)
library(openxlsx)

#--------------------#
#---- PARAMETERS ----#
#--------------------#

reps     <- 20000
burn_in  <- 5000
seeds    <- c(101, 202, 303)   # three independent MCMC runs
pofz_min <- 0.99
z_thresh <- qnorm(0.01)        # one-sided 1% threshold on rubias z-scores

# Reporting group of each lcWGS sampling group, used to label NGSadmix clusters.
# Codes follow the convention of the downstream scripts.
group_of_pop <- c(LAG = "V1", SAB = "V1",
                  EAS = "V2", RUP = "V2", NOT_fall = "V2",
                  NOT_summer = "V3")
group_labels <- c(V1 = "Fall-run North", V2 = "Fall-run South", V3 = "Kâcikâsikumekw")

#-------------------#
#---- FUNCTIONS ----#
#-------------------#

# Harmonize sample IDs across files: hyphens to underscores and zero-padded
# sample number after the community code (e.g. Was_7_ -> Was_007_)
clean_id <- function(x) {
  x <- str_replace_all(x, "-", "_")
  str_replace(x, "([A-Za-z]*s)_(\\d+)_", function(m) {
    g <- str_match(m, "([A-Za-z]*s)_(\\d+)_")
    sprintf("%s_%03d_", g[, 2], as.integer(g[, 3]))
  })
}

run_rubias <- function(reference, mixture, seed) {
  set.seed(seed)
  infer_mixture(reference = reference, mixture = mixture, gen_start_col = 5,
                reps = reps, burn_in = burn_in, method = "MCMC")
}

# Most probable collection for each individual
best_assignment <- function(indiv_posteriors) {
  indiv_posteriors %>%
    group_by(indiv) %>%
    slice_max(PofZ, n = 1, with_ties = FALSE) %>%
    ungroup()
}

# Write a sheet to the supplementary workbook, replacing it if present
write_supp_sheet <- function(data, sheet, file = "res/supplementary_table.xlsx") {
  wb <- if (file.exists(file)) loadWorkbook(file) else createWorkbook()
  if (sheet %in% names(wb)) removeWorksheet(wb, sheet)
  addWorksheet(wb, sheet)
  writeData(wb, sheet, data)
  saveWorkbook(wb, file, overwrite = TRUE)
}

#--------------#
#---- DATA ----#
#--------------#

tmp <- read.table("data/gtscore_analysis/polyGenResults_singleSNP_rubias.txt", header = TRUE,
                  stringsAsFactors = FALSE, colClasses = "character")

tmp %>% count(sample_type)

tmp %<>% mutate(indiv = clean_id(indiv),
                sample_type = if_else(is.na(sample_type), "reference", sample_type))

# NGSadmix assignment (K = 3) of the lcWGS baseline
admix <- read.table("data/ngsadmix/admixture_assign_k3_allpop.txt", header = TRUE,
                    stringsAsFactors = FALSE) %>%
  mutate(indiv = clean_id(individual_id))

# Each NGSadmix cluster takes the reporting group of the majority of its individuals
table(admix$pop, admix$genetic_cluster)

cluster_to_group <- admix %>%
  filter(!is.na(pop)) %>%
  count(genetic_cluster, group = unname(group_of_pop[pop])) %>%
  group_by(genetic_cluster) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()

print(cluster_to_group)
stopifnot(setequal(cluster_to_group$group, names(group_labels)),
          !anyDuplicated(cluster_to_group$group))

admix %<>% left_join(cluster_to_group %>% dplyr::select(genetic_cluster, group),
                     by = "genetic_cluster")

#------------------------------#
#---- BASELINE (REFERENCE) ----#
#------------------------------#

cisco_source <- tmp %>%
  filter(sample_type == "reference") %>%
  dplyr::select(-collection, -repunit) %>%
  left_join(admix %>% dplyr::select(indiv, group), by = "indiv") %>%
  mutate(collection = group, repunit = group) %>%
  dplyr::select(-group) %>%
  relocate(sample_type, repunit, collection, indiv)

gtdev_samples <- read.table("data/id_gtseq_dev.tsv", col.names = "indiv") %>%
  mutate(indiv = clean_id(indiv))

message("lcWGS individuals without GT-seq genotype: ",
        length(setdiff(admix$indiv, cisco_source$indiv)))

# Training set: individuals used to design the GT-seq panel
train <- cisco_source %>%
  filter(indiv %in% gtdev_samples$indiv, !is.na(collection)) %>%
  distinct(indiv, .keep_all = TRUE)

# Holdout set: remaining baseline individuals from groups present in the training set
holdout <- cisco_source %>%
  filter(!indiv %in% gtdev_samples$indiv, !is.na(collection),
         collection %in% train$collection) %>%
  distinct(indiv, .keep_all = TRUE)

train %>% count(collection)
holdout %>% count(collection)

#--------------------------------------------------#
#---- ASSIGNMENT POWER USING THE HOLDOUT SET ----#
#--------------------------------------------------#

# Holdout individuals are analysed as a single mixture of unknown origin;
# their NGSadmix cluster is used only to score the assignments.
holdout_mix <- holdout %>%
  mutate(sample_type = "mixture", repunit = NA_character_, collection = "holdout")

assign_holdout <- run_rubias(train, holdout_mix, seeds[1])

assign_holdout_map <- best_assignment(assign_holdout$indiv_posteriors) %>%
  left_join(holdout %>% dplyr::select(indiv, true_pop = collection), by = "indiv") %>%
  mutate(assignment_res = as.integer(collection == true_pop))

summarise_holdout <- function(d) {
  d %>%
    summarise(accuracy = mean(assignment_res),
              correct = sum(assignment_res),
              total = n(),
              n_PofZ_above_0.99 = sum(PofZ > pofz_min),
              min_PofZ = min(PofZ),
              .groups = "drop")
}

supp_table <- bind_rows(
  assign_holdout_map %>% group_by(true_pop) %>% summarise_holdout() %>%
    mutate(true_pop = unname(group_labels[true_pop])),
  assign_holdout_map %>% summarise_holdout() %>% mutate(true_pop = "All"))

print(supp_table)
write_supp_sheet(supp_table, "GTseq_holdout_assignment")

#-----------------------------------------------------------------#
#---- BASELINE INDIVIDUALS SAMPLED IN RIVERS WITHOUT lcWGS DATA ----#
#-----------------------------------------------------------------#

metadata_source <- read.csv("data/JB_samples_metadata.csv", header = TRUE) %>%
  rename(indiv = 1, pop = 2) %>%
  mutate(indiv = clean_id(indiv)) %>%
  distinct(indiv, .keep_all = TRUE)

unknown_source <- cisco_source %>%
  filter(is.na(collection)) %>%
  dplyr::select(-collection) %>%
  left_join(metadata_source, by = "indiv") %>%
  mutate(sample_type = "mixture", repunit = NA_character_,
         collection = coalesce(pop, "unknown")) %>%
  dplyr::select(-pop) %>%
  relocate(sample_type, repunit, collection, indiv) %>%
  distinct(indiv, .keep_all = TRUE)

unknown_gt <- NULL
if (nrow(unknown_source) > 0) {
  assign_dispersal <- run_rubias(train, unknown_source, seeds[1])
  unknown_gt <- best_assignment(assign_dispersal$indiv_posteriors) %>%
    dplyr::select(indiv, collection)
  best_assignment(assign_dispersal$indiv_posteriors) %>%
    count(sampled_in = mixture_collection, assigned_to = collection) %>%
    print(n = Inf)
}

# Cluster of every baseline individual (used for the GT-seq PCA)
gt_pca <- bind_rows(train %>% dplyr::select(indiv, collection),
                    assign_holdout_map %>% dplyr::select(indiv, collection),
                    unknown_gt)

write.table(gt_pca, "data/baseline_assignment_gtseq.txt", quote = FALSE, row.names = FALSE)

#-----------------------------------------------#
#---- ASSIGNMENT OF HARVESTED INDIVIDUALS ----#
#-----------------------------------------------#

# Reference: training and holdout individuals
ref_all <- cisco_source %>% filter(!is.na(collection))

cisco_mixed <- tmp %>%
  filter(sample_type == "mixture") %>%
  dplyr::select(sample_type, repunit, collection, indiv, starts_with("Chr"))

stopifnot(identical(names(cisco_mixed)[-(1:4)], names(ref_all)[-(1:4)]))

# Three independent runs
mix_runs  <- map(seeds, ~ run_rubias(ref_all, cisco_mixed, .x))
best_runs <- map(mix_runs, ~ best_assignment(.x$indiv_posteriors))

run_agreement <- imap_dfr(best_runs, ~ .x %>%
                            dplyr::select(indiv, collection, PofZ) %>%
                            mutate(run = .y)) %>%
  group_by(indiv) %>%
  summarise(same_assignment = n_distinct(collection) == 1,
            pofz_range = max(PofZ) - min(PofZ),
            .groups = "drop")

message(sprintf("Individuals assigned to the same group in all runs: %d of %d (maximum PofZ range: %.2g)",
                sum(run_agreement$same_assignment), nrow(run_agreement),
                max(run_agreement$pofz_range)))

map_rows <- best_runs[[1]]

#---- z-scores (potential unsampled populations) ----#

ggplot(map_rows, aes(x = z_score)) +
  geom_density(colour = "blue") +
  stat_function(fun = dnorm, colour = "black") +
  geom_vline(xintercept = z_thresh, linetype = 2, color = "red") +
  theme_bw()

ggsave("figures/06_migrating_samples_zscore_distribution.png", height = 10, width = 10, units = "cm", dpi = 600)

#---- Merge with sampling information ----#

geo_location <- read.csv("data/migrating_samples_metadata.csv") %>%
  transmute(id = clean_id(fishes_id), community, river.name, latitude, longitude) %>%
  distinct(id, .keep_all = TRUE)

map_rows %<>%
  left_join(geo_location, by = c("indiv" = "id")) %>%
  dplyr::select(-missing_loci) %>%
  mutate(z_outlier = z_score < z_thresh,
         community = factor(community, levels = c("Whapmagoostui", "Chisasibi", "Wemindji",
                                                  "Eastmain", "Waskaganish")))

map_rows %>%
  ggplot(aes(x = reorder(mixture_collection, latitude), y = z_score)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.1) +
  geom_hline(yintercept = z_thresh, linetype = 2, color = "red") +
  facet_grid(~ community, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

ggsave("figures/06_zscore_per_sampling_site.png", height = 10, width = 10, units = "cm", dpi = 600)

map_rows %>% filter(z_outlier) %>% count(collection, community)

message(sprintf("Harvested individuals assigned: %d", nrow(map_rows)))
message(sprintf("PofZ <= %.2f: %d", pofz_min, sum(map_rows$PofZ <= pofz_min)))
message(sprintf("z-score < %.2f: %d (%.1f expected under the 1%% threshold)",
                z_thresh, sum(map_rows$z_outlier, na.rm = TRUE), 0.01 * nrow(map_rows)))

#---- Save ----#

write.table(map_rows, "data/migrating_individuals_assignment_rubias.txt")

supp_table_assign_res <- map_rows %>%
  mutate(repunit = unname(group_labels[repunit]),
         collection = unname(group_labels[collection]))

write_supp_sheet(supp_table_assign_res, "Mixed_stock_pop")
