#!/bin/bash
#SBATCH -J "04_pca_by_window"
#SBATCH -o log_%j
#SBATCH -c 2
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=2-00:00
#SBATCH --mem=30G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

NB_CPU=2

for i in $(cat 16_sliding_pca/list_window.txt)

do 
	python2 /project/lbernatchez/programs/pcangsd/pcangsd.py -threads $NB_CPU \
        -beagle 16_sliding_pca/beagle_by_window/${i} -o 16_sliding_pca/pca_by_window/${i}
done
