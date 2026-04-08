#!/bin/bash
#SBATCH -J "04_pca"
#SBATCH -o log_%j
#SBATCH -c 5
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=100G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#maybe edit
NB_CPU=5 #change accordingly in SLURM header

module load pcangsd

#prepare variables - avoid to modify
source 01_scripts/01_config.sh

echo "analyse covariance matrix on all individuals"

pcangsd -b 04_ngsParalog/all_maf0.01pctind0.8_maxdepth10_ALL_CHR_canonical.beagle.gz \
	-e 4 \
	-t $NB_CPU \
	-o 06_pca/ALL_CHR_singletons
