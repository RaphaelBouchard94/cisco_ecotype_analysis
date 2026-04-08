#!/bin/bash
#SBATCH -J "08_ngsamova"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=30G


ngsAMOVA --bcf-src 1 \
	--in-vcf 08_ngsamova/all_chromosomes.renamed.bcf \
	-doMajorMinor 1 \
	-doEM 1 \
	-doJGTM 1 \
	-doDist 1 \
	-doAMOVA 1 \
	--print-amova 1 \
	--metadata 08_ngsamova/metadata_reordered.tsv \
	--formula 'Individuals ~ Regions/Populations' \
	-out amova_all_chr

