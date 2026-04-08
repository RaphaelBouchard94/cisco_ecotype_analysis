#!/bin/bash
#SBATCH -J "15_beagle_chr34_outlier_windo"
#SBATCH -o log_%j
#SBATCH -c 4
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=none
#SBATCH --time=1-00:00
#SBATCH --mem=16G


module load angsd

NB_CPU=4


# Run ANGSD to compute Beagle genotype likelihoods for specific region
angsd -P "${NB_CPU}" \
	-doMaf 1  -GL 2 -doGlf 2 -doMajorMinor 5 \
	-anc 02_info/genome.fasta \
	-r CM078251.1:35300000-35519999 \
	-sites 04_ngsParalog/canonical_sites_CM078251.1 \
	-b "02_info/bam.filelist" \
	-out "15_pca_chr34/all_ind_CM078251.1_35300000-35519999"



