#!/bin/bash

## To run:
## parallel -a ./02_info/regions.txt -j 40 srun -c 1 --mem 20G -p small --time 1-00:00 -J 04_ngsParalog_{} -o 04_ngsParalog_{}_%j.log /bin/sh ./01_scripts/04_ngsParalog.sh {} &

module load samtools/1.15

REGION="$1"

SNP_LIST="02_info/sites_all_maf0.01pctind0.8_maxdepth10_"$REGION""
BED_FILE="02_info/sites_all_maf0.01pctind0.8_maxdepth10_"$REGION".bed"
BAM_LIST="02_info/bam.filelist"

### Convert site file to bed
awk '{print $1"\t"$2-1"\t"$2}' $SNP_LIST > $BED_FILE

### Calculate number of individuals to match ANGSD
source 01_scripts/01_config.sh
N_IND=$(wc -l $BAM_LIST | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 

### mpileup and ngsParalog without intermidate files

samtools mpileup -b $BAM_LIST -l $BED_FILE -r $REGION -q 0 -Q 0 --ff UNMAP,DUP |
01_scripts/ngsParalog/ngsParalog calcLR \
    -infile - \
    -outfile 04_ngsParalog/ngsParalog_$REGION \
    -minQ 20 -minind $MIN_IND -mincov $MIN_DEPTH -allow_overwrite 1


