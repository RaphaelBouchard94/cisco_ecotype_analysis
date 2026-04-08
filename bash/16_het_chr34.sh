#!/bin/bash
#SBATCH -J "16_calculate_het_chr34"
#SBATCH -o log_%j
#SBATCH -c 2
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=none
#SBATCH --time=1-00:00
#SBATCH --mem=10G

# Configuration
NB_CPU=2
CHR="${1}"

# Validate inputs
if [[ -z "${CHR}" ]]; then
    echo "Error: Missing required parameter" >&2
    echo "Usage: $0 <chromosome>" >&2
    exit 1
fi

cd "${SLURM_SUBMIT_DIR}" || exit 1
module load angsd/0.931
source 01_scripts/01_config.sh || exit 1

POPULATIONS=("aa_geno" "ab_geno" "bb_geno")

for POP in "${POPULATIONS[@]}"; do
    echo "Computing HWE for population: ${POP}"
    
    angsd -P "${NB_CPU}" \
        -doHWE 1 \
        -doMajorMinor 5 \
	-anc 02_info/genome.fasta \
        -doMaf 1 \
        -GL 2 \
        -r ${CHR} \
        -sites 04_ngsParalog/canonical_sites_${CHR} \
        -b "16_het_chr34/${POP}_bam.filelist" \
        -out "16_het_chr34/${POP}_${CHR}"
done
