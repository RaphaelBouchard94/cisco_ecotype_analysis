#!/bin/bash
#SBATCH -J "13_dxy"
#SBATCH -o 13_dxy_%j
#SBATCH -c 10
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=none
#SBATCH --time=7-00:00
#SBATCH --mem=100G


# Simple bash-only pipeline to calculate dxy in sliding windows
# No Python dependencies except for the final dxy calculation
#

set -e  # Exit on error

module load angsd/0.931

# ============================================================================
# PARAMETERS - MODIFY THESE FOR YOUR DATA
# ============================================================================

# Reference genome
REF="02_info/genome.fasta"

# BAM file lists
POP1_BAM="09_pop_bamfile/notsummer_bam.filelist"
POP2_BAM="09_pop_bamfile/rup_bam.filelist"

# Number of individuals (not chromosomes)
NPOP1=55
NPOP2=59

# Region to analyze
CHROM="CM078251.1"
START=33800000
END=36619999

# Window parameters
WINDOW=20000
STEP=10000

# Output prefix
PREFIX="chr36_outlier"

# Number of threads for ANGSD
THREADS=10

# ============================================================================
# PIPELINE=
# ============================================================================

echo "=========================================="
echo "dxy calculation pipeline (SIMPLE BASH)"
echo "=========================================="
echo "Region: ${CHROM}:${START}-${END}"
echo "Window size: ${WINDOW} bp"
echo "Step size: ${STEP} bp"
echo "Pop1 individuals: ${NPOP1}"
echo "Pop2 individuals: ${NPOP2}"
echo "=========================================="

# Step 1: Create windows BED file using bash
echo ""
echo "[Step 1] Creating windows BED file..."

> 13_dxy/${PREFIX}_windows.bed  # Create empty file

window_start=${START}
window_count=0

while [ ${window_start} -lt ${END} ]; do
    window_end=$((window_start + WINDOW))
    
    # Don't let window extend past region end
    if [ ${window_end} -gt ${END} ]; then
        window_end=${END}
    fi
    
    # Don't create tiny windows at the end
    window_size=$((window_end - window_start))
    if [ ${window_size} -lt ${STEP} ]; then
        break
    fi
    
    echo -e "${CHROM}\t${window_start}\t${window_end}" >> 13_dxy/${PREFIX}_windows.bed
    window_count=$((window_count + 1))
    window_start=$((window_start + STEP))
done

echo "Created ${window_count} windows"

# Step 2: Calculate SAF and 2D-SFS for each window
echo ""
echo "[Step 2] Calculating SAF and 2D-SFS for each window..."

# Initialize the window SFS output file
> 13_dxy/${PREFIX}_window_sfs.txt

processed=0
total_windows=${window_count}

while IFS=$'\t' read -r chrom start end; do
    processed=$((processed + 1))
    region="${chrom}:${start}-${end}"
    window_id="${chrom}_${start}_${end}"
    
    echo "  Processing window ${processed}/${total_windows}: ${region}"
    
    # Calculate SAF for notsummer
    angsd -bam ${POP1_BAM} \
        -anc ${REF} \
        -r ${region} \
        -dosaf 1 \
        -GL 2 \
        -P ${THREADS} \
        -out 13_dxy/${PREFIX}_notsummer_${window_id} \
        2> ${PREFIX}_notsummer_${window_id}.log
    
    # Calculate SAF for rup
    angsd -bam ${POP2_BAM} \
        -anc ${REF} \
        -r ${region} \
        -dosaf 1 \
        -GL 2 \
	-P ${THREADS} \
        -out 13_dxy/${PREFIX}_rup_${window_id} \
        2> ${PREFIX}_rup_${window_id}.log
    
    # Calculate 2D-SFS for this window
    realSFS 13_dxy/${PREFIX}_notsummer_${window_id}.saf.idx 13_dxy/${PREFIX}_rup_${window_id}.saf.idx -P ${THREADS} > 13_dxy/${PREFIX}_2dsfs_${window_id}.txt

    # Flatten the 2D-SFS and append to window SFS file
    # Convert newlines to spaces, remove extra spaces
    sfs_line=$(cat 13_dxy/${PREFIX}_2dsfs_${window_id}.txt | tr '\n' ' ' | sed 's/  */ /g' | sed 's/ $//')
    
    # Write: chrom start end sfs_values
    echo "${chrom} ${start} ${end} ${sfs_line}" >> 13_dxy/${PREFIX}_window_sfs.txt
    
done < 13_dxy/${PREFIX}_windows.bed

echo ""
echo "[Step 2] Complete! Processed ${processed} windows."

# Step 3: Calculate dxy using the Python script
echo ""
echo "[Step 3] Calculating dxy for all windows..."

python3 01_scripts/calculatedxy.py \
    -w 13_dxy/${PREFIX}_window_sfs.txt \
    -m ${NPOP1} \
    -n ${NPOP2} \
    -o 13_dxy/${PREFIX}_dxy_results.txt

