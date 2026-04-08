#!/bin/bash
#SBATCH -J "08_ngsamova"
#SBATCH -o log_%j
#SBATCH -c 4
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=8G


set -e 

INPUT_DIR="08_ngsamova"
OUTPUT_DIR="08_ngsamova"
THREADS=4

echo "Step 1: Converting VCF to BCF for each chromosome..."

# Array to store BCF files
BCF_FILES=()

# Loop through each chromosome VCF
for vcf in ${INPUT_DIR}/*.vcf.gz; do
    
    filename=$(basename "$vcf" .vcf.gz)
    bcf="${OUTPUT_DIR}/${filename}.bcf"
    
    echo "  Processing: $filename" 
    # Convert VCF to BCF
    bcftools view -O b --threads $THREADS -o "$bcf" "$vcf"
    
    # Index the BCF
    bcftools index -f "$bcf"
    
    # Add to array
    BCF_FILES+=("$bcf")
done

echo ""
echo "Step 2: Merging all chromosomes..."

# Merge all BCF files
bcftools concat -O b --threads $THREADS -o "${OUTPUT_DIR}/all_chromosomes.bcf" "${BCF_FILES[@]}"

echo ""
echo "Step 3: Indexing merged file..."

# Index the merged BCF
bcftools index "${OUTPUT_DIR}/all_chromosomes.bcf"

echo ""
echo "Done! Merged BCF: ${OUTPUT_DIR}/all_chromosomes.bcf"

# Show summary
echo ""
echo "Summary:"
bcftools stats "${OUTPUT_DIR}/all_chromosomes.bcf" | grep "^SN" | head -10

