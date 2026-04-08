#!/bin/bash
#SBATCH -J "14_make_window_ngsLD"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=X
#SBATCH --time=1-00:00
#SBATCH --mem=100G


# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#module
module load angsd
ulimit -S -n 2048



source 01_scripts/01_config.sh

i="CM078251.1"

#run Eric's optimized script on compressed file
#gzip $LDOUT_FILE

#cat 12_ngsLD/header.gz 12_ngsLD/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i".ld.gz > 12_ngsLD/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i"_header.ld.gz

python3 01_scripts/utility_scripts/ld_by_blocks_optimized_gzinput.py 14_ngsld/karyobb_CM078251.1.ld 50 14_ngsld/karyobb_CM078251_by_50.ld


