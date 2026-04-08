#!/bin/bash
# To run:
# parallel -a ./07_ngsadmix/k.txt -j 5 srun -p medium -c 20 --mem=200G --time=7-00:00 -J k_{} -o ngsadmixk_%j.log 01_scripts/07_ngsadmix.sh {} &

# Configuration
NB_CPU=20  # Change accordingly in SLURM header
K="${1}"

# Validate input
if [[ -z "${K}" ]]; then
    echo "Error: K value not provided" >&2
    echo "Usage: $0 <K_value>" >&2
    exit 1
fi

# Important: Move to directory where job was submitted
cd "${SLURM_SUBMIT_DIR}" || exit 1

# Prepare variables - avoid modifying
source 01_scripts/01_config.sh || exit 1
source ~/.bashrc

# Load required module
module load angsd

# Run NGSadmix iterations
for i in {1..50}; do
    echo "Running iteration ${i} for K=${K}"
    
    NGSadmix \
        -P "${NB_CPU}" \
        -seed 2358 \
        -likes 04_ngsParalog/all_maf0.01pctind0.8_maxdepth10_ALL_CHR_canonical.beagle.gz \
        -K "${K}" \
        -tol 0.0000001 \
        -maxiter 20000 \
        -outfiles "07_ngsadmix/all_maf${MIN_MAF}_pctind${PERCENT_IND}_maxdepth${MAX_DEPTH_FACTOR}_singletons_${K}_${i}"
    
    # Check if NGSadmix succeeded
    if [[ $? -ne 0 ]]; then
        echo "Warning: NGSadmix failed for iteration ${i}" >&2
    fi
done

echo "Completed all 50 iterations for K=${K}"

