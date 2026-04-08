library(tidyverse)
library(ggplot2)
library(magrittr)

# Set working directory to where your .cov files are located
setwd("data/chr34/sliding_pca/")

# Get all .cov files for chromosome CM078251.1 (chr34)
cov_files <- list.files(pattern = "window_CM078251\\.1_.*\\.cov$")

# Function to extract window information from filename
extract_window_info <- function(filename) {
  # Extract chromosome, start, and end positions
  # Pattern: window_CHR_START-END_1000_snps.cov
  pattern <- "window_(CM078251\\.1)_(\\d+)-(\\d+)_500_snps\\.cov"
  matches <- str_match(filename, pattern)
  
  data.frame(
    filename = filename,
    chr = matches[2],
    start = as.numeric(matches[3]),
    end = as.numeric(matches[4]),
    midpoint = (as.numeric(matches[3]) + as.numeric(matches[4])) / 2
  )
}

# Function to perform PCA on a single covariance matrix
perform_pca <- function(cov_matrix) {
  # Eigendecomposition of covariance matrix
  eigen_result <- eigen(cov_matrix)
  
  # PC scores (eigenvectors scaled by sqrt of eigenvalues)
  pc_scores <- eigen_result$vectors
  
  # Variance explained
  var_explained <- eigen_result$values / sum(eigen_result$values) * 100
  
  list(
    pc_scores = pc_scores,
    var_explained = var_explained
  )
}

# Process all files
results_list <- list()

for (i in seq_along(cov_files)) {
  file <- cov_files[i]
  
  # Read covariance matrix
  cov_matrix <- as.matrix(read.table(file))
  
  # Perform PCA
  pca_result <- perform_pca(cov_matrix)
  
  # Extract window information
  window_info <- extract_window_info(file)
  
  # Create data frame with PC1 scores for each individual
  pc_data <- data.frame(
    individual = 1:nrow(cov_matrix),
    PC1 = pca_result$pc_scores[, 1],
    PC2 = pca_result$pc_scores[, 2],  # in case you want PC2 later
    var_PC1 = pca_result$var_explained[1],
    var_PC2 = pca_result$var_explained[2],
    chr = window_info$chr,
    start = window_info$start,
    end = window_info$end,
    midpoint = window_info$midpoint,
    window = paste0(window_info$start, "-", window_info$end)
  )
  
  results_list[[i]] <- pc_data
  
  # Progress indicator
  if (i %% 10 == 0) cat("Processed", i, "of", length(cov_files), "files\n")
}

# Combine all results
all_pca_results <- bind_rows(results_list)

# Sort by genomic position
all_pca_results <- all_pca_results %>%
  arrange(start, individual)

# Save the combined results
write_csv(all_pca_results, "chr34_windowed_pca_results.csv")

##############################################################3
#Add metadata to file

id <- read.table("../../id_part3.tsv")

all_pca_results$fishes_id <- rep(id$V1,nrow(all_pca_results)/nrow(id))

#pop column
metadata <- read.table("../karyo_info.tsv",header=T)

all_pca_results %<>% left_join(metadata, by = c("fishes_id" = "id"))

targ_win <- c(34800000,35619999)
targ_win_plus <- c(34800000-1000000,35619999+1000000)


order_pca_res <- all_pca_results %>%
  arrange(geno_clus, fishes_id) %>%
  mutate(fishes_id = factor(fishes_id, levels = unique(fishes_id)))

# PLOTTING: PC1 values along chromosome for each individual
# Option 1: Line plot for each individual
p1 <- ggplot(order_pca_res  %>% filter(start > targ_win_plus[1] & end < targ_win_plus[2]), aes(x = midpoint, y = PC1, color = fishes_id))+
  geom_line(alpha = 0.7) +
  annotate("rect", 
           xmin = targ_win[1], xmax = targ_win[2], 
           ymin = -0.3, ymax = 0.25,
           fill = "indianred", 
           alpha = 0.2) +
  scale_color_viridis_d(direction = -1)+
  labs(
    y = "PC1 Score",
    x = ""
  ) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = "none")

p1

setwd("../../../")

ggsave("figures/05_windowpca_chr34_outlier.png", p1,
       width = 10, height = 10, unit = "cm", dpi = 900)




