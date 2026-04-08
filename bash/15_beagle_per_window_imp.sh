#!/bin/bash
#SBATCH -J "15_pca_per_window"
#SBATCH -o log_%j.out
#SBATCH -e log_%j.err
#SBATCH -c 4 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=none
#SBATCH --time=7-00:00:00
#SBATCH --mem=5G

################################################################################
# PCA Analysis Per Window
# This script splits a beagle file into windows and calculates the covariance 
# matrix for each window using pcangsd
################################################################################

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# Configuration
readonly NB_CPU=4
readonly WINDOW_SIZE=500
readonly CHROMOSOME="CM078251.1"

# Paths
readonly CONFIG_SCRIPT="01_scripts/01_config.sh"
readonly PYTHON_SPLIT_SCRIPT="01_scripts/utility_scripts/beagle_sliding_window.py"
readonly INPUT_BEAGLE="04_ngsParalog/all_maf0.01pctind0.8_maxdepth10_${CHROMOSOME}_canonical.beagle.gz"
readonly OUTPUT_BASE="15_window_pca_chr34/${WINDOW_SIZE}"
readonly BEAGLE_WINDOW_DIR="${OUTPUT_BASE}/beagle_by_window"
readonly COV_WINDOW_DIR="${OUTPUT_BASE}/cov_by_window"

# Function: Print error message and exit
error_exit() {
    echo "ERROR: $1" >&2
    exit 1
}

# Function: Check if file/directory exists
check_exists() {
    if [[ ! -e "$1" ]]; then
        error_exit "$1 does not exist"
    fi
}

# Validate prerequisites
echo "=== Validating prerequisites ==="
check_exists "$CONFIG_SCRIPT"
check_exists "$PYTHON_SPLIT_SCRIPT"
check_exists "$INPUT_BEAGLE"

# Source configuration
echo "=== Loading configuration ==="
source "$CONFIG_SCRIPT" || error_exit "Failed to source config script"

# Load pcangsd
module load pcangsd

# Create output directories
echo "=== Creating output directories ==="
mkdir -p "$BEAGLE_WINDOW_DIR" "$COV_WINDOW_DIR"

# Split beagle into non-overlapping windows
echo "=== Splitting beagle file into ${WINDOW_SIZE}bp windows ==="
python "$PYTHON_SPLIT_SCRIPT" \
    "$INPUT_BEAGLE" \
    "$WINDOW_SIZE" \
    "$BEAGLE_WINDOW_DIR/" || error_exit "Failed to split beagle file"

# Count windows created
WINDOW_COUNT=$(find "$BEAGLE_WINDOW_DIR" -type f | wc -l)
echo "=== Created $WINDOW_COUNT windows ==="

# Run pcangsd on each window
echo "=== Running pcangsd on each window ==="
PROCESSED=0
FAILED=0

find "$BEAGLE_WINDOW_DIR" -type f -name "*.beagle.gz" | sort | while read -r window_file; do
    window_name=$(basename "$window_file")
    
    echo "Processing window: $window_name"
    
    # Run pcangsd
    if pcangsd --threads "$NB_CPU" -b "$window_file"  -o "$COV_WINDOW_DIR/${window_name%.beagle.gz}"; then
        ((PROCESSED++)) || true
        echo "  ✓ Successfully processed $window_name"
    else
        ((FAILED++)) || true
        echo "  ✗ Failed to process $window_name" >&2
    fi
done

# Summary
echo "=== Analysis Complete ==="
echo "Total windows: $WINDOW_COUNT"
echo "Successfully processed: $PROCESSED"
echo "Failed: $FAILED"

if [[ $FAILED -gt 0 ]]; then
    echo "WARNING: Some windows failed to process" >&2
    exit 1
fi

echo "=== All done! ==="
