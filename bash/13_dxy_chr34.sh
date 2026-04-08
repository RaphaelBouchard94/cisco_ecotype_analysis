#!/bin/bash
#SBATCH -J "13_dxy_chr34"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=7-00:00
#SBATCH --mem=50G

module load angsd

REGION=$1

if [[ -z "$REGION" ]]; then
	echo "Error missing parameters" >&2
	echo "Usage $0 <region>" >&2
	exit 1
fi

winsfs -v --output 13_dxy/notsummer_rup_"$REGION".sfs \
	-t 10 \
	10_saf_maf_by_pop/per_chr/notsummer_"$REGION".saf.idx \
	10_saf_maf_by_pop/per_chr/rup_"$REGION".saf.idx


# 2. Run winsfs with window size (-w) and step size (-s)
