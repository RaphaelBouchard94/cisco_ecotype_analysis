#!/bin/bash
#SBATCH -J "15_pcaoutlier"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=10G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#maybe edit
NB_CPU=1 #change accordingly in SLURM header

module load pcangsd

#prepare variables - avoid to modify
source 01_scripts/01_config.sh

echo "analyse covariance matrix on all individuals"

pcangsd -b 15_pca_chr34/all_ind_CM078251.1_35300000-35519999.beagle.gz \
        -e 4 \
        -t $NB_CPU \
	-o 15_pca_chr34/outlier_last_section

