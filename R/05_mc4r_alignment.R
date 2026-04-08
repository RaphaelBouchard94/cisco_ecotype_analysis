
library(ggmsa)
library(Biostrings)
library(tidyverse)
library(msa)
library(ape)
library(phangorn)
library(seqinr)
library(Peptides)
library(zoo)


# 1. Read alignment
aa_data <- read.alignment("data/mc4r/aligned_salmo_human_cisco.fa", format = "fasta")

names(mc4r_aa_aln) <- aa_data$nam


sequences <- sapply(aa_data$seq, paste, collapse = "")
mc4r_aa_aln <- AAStringSet(sequences)
names(mc4r_aa_aln) <- c("Homo sapiens",
                        "Kâcikâsikumekw",
                        "Nûtamesânîw",
                        "Salmo salar",
                        "Salvelinus namaycush",
                        "Oncorhynchus gorbuscha",
                        "Oncorhynchus tshawytscha")


# 5. Define region of interest
aln_pos <- 284 # Change this to your position of interest
window <- 10

# 6. Create alignment visualization
p <- ggmsa(
  mc4r_aa_aln,
  start = aln_pos - window,
  end   = aln_pos + window,
  color = "Hydrophobicity",  # Options: "Chemistry", "Hydrophobicity", "Clustal", "Taylor"
  char_width = 0.5,
  seq_name = TRUE,consensus_views = T
) +
  theme_bw(base_size = 18)+
  theme(axis.text.y = element_text(face = "italic",size = 14))

print(p)

# 7. Save the plot
ggsave("figures/05_mc4r_alignment.png", p, width = 8, height = 6)



library(ggplot2)
library(zoo)

# 1. Define the standard Kyte-Doolittle scale
kd_scale <- c(
  'A' = 1.8,  'R' = -4.5, 'N' = -3.5, 'D' = -3.5, 'C' = 2.5, 
  'Q' = -3.5, 'E' = -3.5, 'G' = -0.4, 'H' = -3.2, 'I' = 4.5, 
  'L' = 3.8,  'K' = -3.9, 'M' = 1.9,  'F' = 2.8,  'P' = -1.6, 
  'S' = -0.8, 'T' = -0.7, 'W' = -0.9, 'Y' = -1.3, 'V' = 4.2
)

# 2. Function to map residues to scores and calculate rollmean
calc_local_hydro <- function(sequence, window) {
  # Split sequence into individual letters
  residues <- strsplit(as.character(sequence), "")[[1]]
  # Map residues to KD scores (replace non-standard with 0)
  scores <- kd_scale[residues]
  scores[is.na(scores)] <- 0
  
  # Calculate sliding window average
  rollmean(scores, k = window, fill = NA, align = "center")
}

# 3. Apply to your mc4r_aa_aln (AAStringSet)
# We use window = 5 to see the local impact of a single mutation
hydro_list <- lapply(as.character(mc4r_aa_aln), calc_local_hydro, window = 2)

# 4. Prepare data for ggplot
plot_df <- data.frame(
  Position = rep(1:width(mc4r_aa_aln)[1], length(hydro_list)),
  Species = rep(names(mc4r_aa_aln), each = width(mc4r_aa_aln)[1]),
  Hydropathy = unlist(hydro_list)
)

# 5. Zoom into your region of interest (aln_pos 284)
subset_df <- plot_df[plot_df$Position >= (284 - 15) & plot_df$Position <= (284 + 15), ]

# 6. Plot
ggplot(subset_df, aes(x = Position, y = Hydropathy, color = Species)) +
  geom_line() +
  theme_minimal() +
  geom_vline(xintercept = 284, linetype = "dotted") +
  labs(y = "Hydropathy (Kyte-Doolittle)", x = "Alignment Position")



# 1. Prepare your highlighted species list
highlight_species <- c("Kuukamek", "Nutamesanio")

# 2. Create the plot
ggplot(subset_df, aes(x = Position, y = Hydropathy, group = Species)) +
  # Draw all OTHER species in grey
  geom_line(data = subset(subset_df, !(Species %in% highlight_species)), 
            color = "grey30", size = 0.8, alpha = 0.5) +
  # Draw Kuukamek and Nutamesanio with specific colors and thicker lines
  geom_line(data = subset(subset_df, Species %in% highlight_species), 
            aes(color = Species), size = 1) + 
  # Set your specific colors
  scale_color_manual(values = c("Kuukamek" = "#E69F00",    # Strong Orange/Red
                                "Nutamesanio" = "skyblue1")) + # Strong Green
  geom_vline(xintercept = 284, linetype = "dotted", color = "black") +
  theme_bw(base_size = 16) +
  labs(y = "Hydropathy (Kyte-Doolittle)", 
       x = "Alignment Position") +
  theme(legend.position = "none",
        panel.grid = element_blank())

ggsave("figures/04_hydropathy.png",width = 10, height = 10, unit = "cm", dpi = 600)


subset_df %>% filter(Position == 285)


