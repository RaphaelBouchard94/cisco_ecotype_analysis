#!/bin/bash
#SBATCH -J "15_pca_per_window"
#SBATCH -o log_%j
#SBATCH -c 4 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=none
#SBATCH --time=7-00:00
#SBATCH --mem=5G

###this script will split the beagle into windows and calculate the covariance matrix of each window
#maybe edit
NB_CPU=4 #change accordingly in SLURM header
window_size=1000 # for small window along a chromosome for instance

#prepare variables - avoid to modify
source 01_scripts/01_config.sh

####split beagle into non overlapping windows
mkdir -p "14_pca_chr34/"$window_size
mkdir -p "14_pca_chr34/"$window_size"/beagle_by_window"
mkdir -p "14_pca_chr34/"$window_size"/cov_by_window"


#if your reference contigs names are "Chr1" use this line
python 01_scripts/utility_scripts/beagle_sliding_window.py 04_ngsParalog/all_maf0.01pctind0.8_maxdepth10_CM078251.1_canonical.beagle.gz $window_size 14_pca_chr34/"$window_size"/beagle_by_window/

#run pcangsd on each window
#this is the input file for the pca
ls -1 10_pca_by_window/"$window_size"/beagle_by_window/ |
    sort -u |
    while read i
    do
        echo "GL for $i"
INPUT=14_pca_chr34/"$window_size"/beagle_by_window/$i

echo "analyse covariance matrix on all individuals"
python3 $PCA_ANGSD_PATH/pcangsd.py -threads $NB_CPU -beagle $INPUT -o 14_pca_chr34/"$window_size"/cov_by_window/$i
done

