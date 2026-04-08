#!/bin/bash

#parallel -a 02_info/regions.txt -j 10 srun -p medium -c 10 --mem=50G --time=12:00:00 -J saf_{} -o saf_{}_%j.log /bin/sh 01_scripts/09_calculate_maf_per_chr.sh notsummer {} &
#parallel -a 02_info/regions.txt -j 10 srun -p medium -c 10 --mem=50G --time=12:00:00 -J saf_{} -o saf_{}_%j.log /bin/sh 01_scripts/09_calculate_maf_per_chr.sh notfall {} &
#parallel -a 02_info/regions.txt -j 10 srun -p medium -c 10 --mem=50G --time=12:00:00 -J saf_{} -o saf_{}_%j.log /bin/sh 01_scripts/09_calculate_maf_per_chr.sh rup {} &

# Configuration
NB_CPU=10  # Change accordingly in SLURM header
POP="${1}"
CHR="${2}"  # New parameter for chromosome

# Validate inputs
if [[ -z "${POP}" ]] || [[ -z "${CHR}" ]]; then
    echo "Error: Missing required parameters" >&2
    echo "Usage: $0 <population> <chromosome>" >&2
    exit 1
fi

# Important: Move to directory where job was submitted
cd "${SLURM_SUBMIT_DIR}" || exit 1

# Load required module
module load angsd/0.931

# Prepare variables
source 01_scripts/01_config.sh || exit 1

# Create output directory if needed
mkdir -p 10_saf_maf_by_pop/per_chr

echo "Computing SAF for population ${POP}, chromosome ${CHR}"


# Run ANGSD to compute MAF per population with consistent major allele
angsd -P "${NB_CPU}" \
    -GL 2 \
    -doMajorMinor 3 \
    -doMaf 1 \
    -anc 02_info/genome.fasta \
    -r ${CHR} \
    -sites 04_ngsParalog/canonical_sites_${CHR} \
    -b "09_pop_bamfile/${POP}_bam.filelist" \
    -out "10_saf_maf_by_pop/per_chr/${POP}_${CHR}"


# Check if ANGSD succeeded
if [[ $? -eq 0 ]]; then
    echo "Successfully completed SAF calculation for ${POP} ${CHR}"
else
    echo "Error: ANGSD failed for ${POP} ${CHR}" >&2
    exit 1
fi
