#!/bin/bash
#SBATCH -J "10_fst"
#SBATCH -o log_%j
#SBATCH -c 10
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=8G

# Configuration
NB_CPU=10
POP1="${1}"
POP2="${2}"

# Validate inputs
if [[ -z "${POP1}" ]] || [[ -z "${POP2}" ]]; then
    echo "Error: Missing required parameters" >&2
    echo "Usage: $0 <pop1> <pop2>" >&2
    exit 1
fi

cd "${SLURM_SUBMIT_DIR}" || exit 1
module load angsd/0.931

source 01_scripts/01_config.sh || exit 1

# Create output directory
mkdir -p 11_fst

echo "Computing genome-wide 2D-SFS for ${POP1} vs ${POP2}"

# Calculate 2D-SFS
realSFS \
    10_saf_maf_by_pop/${POP1}_genome_wide_saf.saf.idx \
    10_saf_maf_by_pop/${POP2}_genome_wide_saf.saf.idx \
    -P ${NB_CPU} \
    > 11_fst/${POP1}_${POP2}_genome_wide.ml

if [[ $? -ne 0 ]]; then
    echo "Error: 2D-SFS calculation failed" >&2
    exit 1
fi

echo "Computing genome-wide FST for ${POP1} vs ${POP2}"

# Calculate FST index
realSFS fst index \
    10_saf_maf_by_pop/${POP1}_genome_wide_saf.saf.idx \
    10_saf_maf_by_pop/${POP2}_genome_wide_saf.saf.idx \
    -sfs 11_fst/${POP1}_${POP2}_genome_wide.ml \
    -fstout 11_fst/${POP1}_${POP2}_genome_wide \
    -P ${NB_CPU}

if [[ $? -ne 0 ]]; then
    echo "Error: FST index calculation failed" >&2
    exit 1
fi

echo "Calculating global FST estimate"

# Get global FST estimate
realSFS fst stats \
    11_fst/${POP1}_${POP2}_genome_wide.fst.idx \
    > 11_fst/${POP1}_${POP2}_genome_wide.global.fst

echo "Calculating windowed FST estimates"

# Get windowed FST 
realSFS fst stats2 \
    11_fst/${POP1}_${POP2}_genome_wide.fst.idx \
    -win 25000 -step 5000 \
    > 11_fst/${POP1}_${POP2}_genome_wide.25kb_windows.fst

if [[ $? -eq 0 ]]; then
    echo "Successfully computed FST for ${POP1} vs ${POP2}"
    echo ""
    echo "=== Global FST ==="
    cat 11_fst/${POP1}_${POP2}_genome_wide.global.fst
    echo ""
    echo "Output files:"
    ls -lh 11_fst/${POP1}_${POP2}_genome_wide.*
else
    echo "Error: FST calculation failed" >&2
    exit 1
fi
