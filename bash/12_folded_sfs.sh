#!/bin/bash
#SBATCH -J "12_folded_sfs"
#SBATCH -o log_%j
#SBATCH -c 10
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=10G


POP="$1"

realSFS 10_saf_maf_by_pop/${POP}_genome_wide_saf.saf.idx -fold 1 > 12_theta/${POP}_genome_wide.sfs

