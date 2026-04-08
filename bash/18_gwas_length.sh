#!/bin/bash
#SBATCH -J "04_pca"
#SBATCH -o log_%j
#SBATCH -c 20
#SBATCH -p large
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=7-00:00
#SBATCH --mem=100G

module load angsd

angsd -b 02_info/bam.filelist \
  -ref 02_info/genome.fasta \
  -rf 02_info/regions.txt \
  -out 18_gwas/gwas_length \
  -doAsso 2 \
  -doMaf 1 \
  -doPost 1 \
  -SNP_pval 1e-6 \
  -GL 2 \
  -doMajorMinor 3 \
  -sites 04_ngsParalog/canonical_sites_ALL_CHR \
  -yQuant 18_gwas/pheno.yquant \
  -Pvalue 1 \
  -P 20  
