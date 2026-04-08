#!/bin/bash
#SBATCH -J "04_pca"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=x
#SBATCH --time=1-00:00
#SBATCH --mem=25G


cd $SLURM_SUBMIT_DIR

# LOAD MODULES
module load angsd

ulimit -S -n 2048

angsd -bam 16_het_chr34/aa_geno_bam.filelist  -dohaplocall 1 -doCounts 1 -r CM078251.1:35059641-35060636 -out 17_make_fasta/aa

angsd -bam 16_het_chr34/aa_geno_bam.filelist  -dohaplocall 1 -doCounts 1 -r CM078251.1:35059641-35062636 -out 17_make_fasta/aa_promo

angsd -bam 16_het_chr34/bb_geno_bam.filelist  -dohaplocall 1 -doCounts 1 -r CM078251.1:35059641-35060636 -out 17_make_fasta/bb

angsd -bam 16_het_chr34/bb_geno_bam.filelist  -dohaplocall 1 -doCounts 1 -r CM078251.1:35059641-35062636 -out 17_make_fasta/bb_promo
