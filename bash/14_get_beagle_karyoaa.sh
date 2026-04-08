#!/bin/bash
#srun -p medium -c 10 --mem=50G --time=12:00:00 -J beagle_notsummer_ch34 -o beagle_notsummer_chr34.log /bin/sh 01_scripts/14_get_beagle_karyoaa.sh karyo_aa CM078251.1 &

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

echo "Computing SAF for population ${POP}, chromosome ${CHR}"

# Run ANGSD to compute Beagle genotype likelihoods
angsd -P "${NB_CPU}" \
    -doMaf 1  -GL 2 -doGlf 2 -doMajorMinor 5 \
    -anc 02_info/genome.fasta \
    -r ${CHR}:33800000-36619999\
    -sites 04_ngsParalog/canonical_sites_${CHR} \
    -b "16_het_chr34/aa_geno_bam.filelist" \
    -out "14_ngsld/${POP}_${CHR}"


# Check if ANGSD succeeded
if [[ $? -eq 0 ]]; then
    echo "Successfully completed SAF calculation for ${POP} ${CHR}"
else
    echo "Error: ANGSD failed for ${POP} ${CHR}" >&2
    exit 1
fi

