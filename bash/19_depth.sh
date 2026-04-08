#!/bin/bash
#SBATCH -J "mosdepth"
#SBATCH -o log_%j
#SBATCH -c 4
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=none
#SBATCH --time=01-00:00
#SBATCH --mem=5G

module load mosdepth
mkdir -p 19_depth_chr36

while read id
do
    bam=$(grep "$id" 02_info/bam.filelist | cut -f2)
    echo "Processing: $id -> $bam"  # Debug line
    mosdepth 19_depth_chr36/"$id" "$bam" \
        --fast-mode \
        --threads 4 \
        --b 14_ngsld/chr_36_outlier.bed
done < 02_info/was.id
