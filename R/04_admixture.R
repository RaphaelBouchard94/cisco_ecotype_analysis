#########################
####### LIBRARY #########
#########################
library(tidyverse)
library(magrittr)
library(paletteer)
library(patchwork)

####################################
######  Admixture analysis  ########
####################################

# Load metadata
id <- read.table("data/id_part3.tsv")
metadata <- read.csv("data/data_with_gillrakers.csv", header = TRUE)
id %<>% left_join(metadata, by = c("V1" = "FISHES_id"))
id %<>% mutate(pop = if_else(is.na(pop), "RUP", pop))
id$pop <- factor(id$pop, levels = c("RUP", "NOT_fall", "NOT_summer"))

# Import ngsadmix proportion of admixture
files <- list.files(path = "data/admixture", 
                    pattern = "qopt",
                    full.names = TRUE)

# Define colors for dominant ancestry in each population
pop_colors <- c(
  "NOT_summer" = "#E69F00",
  "NOT_fall" = "skyblue1",
  "RUP" = "skyblue3"
)

# Loop over K to plot admixture plot
kfig <- vector('list', length(files) + 1)

for (i in 2:(length(files) + 1)) {
  df <- read.table(paste("data/admixture/all_maf0.01_pctind0.8_maxdepth10_singletons_", i, "_1.qopt", sep = ""))
  
  df$id <- id$V1
  
  # Join with metadata FIRST
  df <- df %>% 
    left_join(id, by = c("id" = "V1"))
  
  # Then reshape
  dfclean <- df %>% 
    gather('ancestry_cluster', 'prob', starts_with("V")) %>%
    group_by(id) %>% 
    mutate(likely_assignment = ancestry_cluster[which.max(prob)],
           assignment_prob = max(prob)) %>% 
    ungroup()
  
  # Identify dominant ancestry cluster for each population
  dominant_ancestry <- dfclean %>%
    group_by(pop, ancestry_cluster) %>%
    summarize(mean_prob = mean(prob), .groups = "drop") %>%
    group_by(pop) %>%
    slice_max(mean_prob, n = 1) %>%
    select(pop, dominant_cluster = ancestry_cluster)
  
  # Create color mapping
  # Start with a palette that has enough colors
  base_palette <- paletteer_d("ggsci::category20_d3", n = i)
  
  # Get all unique ancestry clusters
  all_clusters <- unique(dfclean$ancestry_cluster)
  
  # Initialize color vector
  cluster_colors <- setNames(rep(NA, length(all_clusters)), all_clusters)
  
  # Assign colors to dominant clusters first
  for (pop_name in names(pop_colors)) {
    dom_cluster <- dominant_ancestry %>%
      filter(pop == pop_name) %>%
      pull(dominant_cluster)
    
    if (length(dom_cluster) > 0) {
      cluster_colors[dom_cluster] <- pop_colors[pop_name]
    }
  }
  
  # Assign remaining colors to non-dominant clusters
  remaining_clusters <- all_clusters[is.na(cluster_colors)]
  remaining_colors <- base_palette[!base_palette %in% cluster_colors]
  cluster_colors[remaining_clusters] <- remaining_colors[1:length(remaining_clusters)]
  
  # Order individuals within each population by their assignment probability
  dfclean <- dfclean %>%
    arrange(pop, likely_assignment, desc(assignment_prob)) %>%
    mutate(id = forcats::fct_inorder(factor(id)))
  
  kfig[[i]] <-
    ggplot(dfclean, aes(id, prob, fill = ancestry_cluster)) +
    geom_col(width = 1) +
    facet_grid(~ pop, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = cluster_colors) +
    labs(y = paste("K =", i), x = "") +
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 18),
          axis.ticks.x = element_blank(),
          strip.text = if (i > 2) element_blank() else element_text(size = 18),
          legend.position = "none")
}

wrap_plots(kfig[2:3], ncol = 1)

ggsave("figures/04_admixture.png",height=15,width =27, unit = "cm", dpi = 600)






