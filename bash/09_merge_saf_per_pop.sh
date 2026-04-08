#!/bin/bash
#SBATCH -J "09_merge_saf"
#SBATCH -o log_%j
#SBATCH -c 4
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=16G

# Configuration
POP="${1}"

# Validate inputs
if [[ -z "${POP}" ]]; then
    echo "Error: Missing required parameter" >&2
    echo "Usage: $0 <population>" >&2
    exit 1
fi

module load angsd/0.931

source 01_scripts/01_config.sh || exit 1

# Create output directory
mkdir -p 10_saf_maf_by_pop

echo "Merging SAF files for population ${POP}"

# Merge all chromosome SAF files for this population
realSFS cat \
    10_saf_maf_by_pop/per_chr/${POP}_*.saf.idx \
    -outnames 10_saf_maf_by_pop/${POP}_genome_wide_saf

if [[ $? -eq 0 ]]; then
    echo "Successfully merged SAF files for ${POP}"
    # List output files
    ls -lh 10_saf_maf_by_pop/${POP}_genome_wide_saf.*
else
    echo "Error: SAF merge failed for ${POP}" >&2
    exit 1
fi
