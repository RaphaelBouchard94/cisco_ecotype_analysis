#!/bin/bash

#parallel -a ./02_info/regions.txt -j 10 srun -p medium -c 1 --mem=20G --time=7-00:00 -J 08_dovcf_{} -o  08_dovcf_{}_%j.log /bin/sh 01_scripts/08.0_generate_vcf.sh {} &


# VARIABLES
REGIONS="$1" # positional argument 1, refers to 1st arg provided when calling the script (03_saf_maf_gl_all_parallel_LL.sh HERE)
NB_CPU=1

# LOAD MODULES
module load angsd/0.937

ulimit -S -n 2048

#============================================
# SOURCE CONFIG AND CALCULATE PARAMETERS
#============================================
source 01_scripts/01_config.sh

# Calculate individual and depth filters
N_IND=$(wc -l < 02_info/bam.filelist)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)" | bc -l)
MIN_IND=${MIN_IND_FLOAT%.*}
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" | bc -l)

# Log parameters
echo "================================================"
echo "Processing region: $REGIONS"
echo "Number of individuals: $N_IND"
echo "Minimum individuals (${PERCENT_IND}%): $MIN_IND"
echo "Maximum depth: $MAX_DEPTH"
echo "Minimum MAF: $MIN_MAF"
echo "================================================"

#============================================
# RUN ANGSD
#============================================

OUTPUT_PREFIX="08_ngsamova/all_maf${MIN_MAF}_pctind${PERCENT_IND}_maxdepth${MAX_DEPTH_FACTOR}_${REGIONS}"

echo "Processing region: ${REGIONS}"
echo "Output prefix: ${OUTPUT_PREFIX}"
echo "CPUs: ${NB_CPU}"

#Generate vcf file
angsd -P $NB_CPU \
	-nQueueSize 50 \
	-doBcf 1 \
	-doMaf 1  \
	-GL 2 \
	-doPost 1 \
	-doMajorMinor 1 \
	-sites 04_ngsParalog/canonical_sites_${REGIONS} \
	-anc 02_info/genome.fasta \
	-b 02_info/bam.filelist \
	-r ${REGIONS} \
	-out ${OUTPUT_PREFIX} 

# Check if ANGSD completed successfully
if [ $? -eq 0 ]; then
    echo "ANGSD completed successfully for $REGIONS"
    # Compress and index VCF if created
    if [ -f "${OUTPUT_PREFIX}.vcf" ]; then
        bgzip -f ${OUTPUT_PREFIX}.vcf
        tabix -p vcf ${OUTPUT_PREFIX}.vcf.gz
        echo "VCF compressed and indexed for $REGIONS"
    fi
else
    echo "ERROR: ANGSD failed for $REGIONS"
    exit 1
fi

echo "Processing complete for $REGIONS"
