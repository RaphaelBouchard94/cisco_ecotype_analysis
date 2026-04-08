#!/bin/bash

#parallel -a ./02_info/regions.txt -j 30 srun -p small -c 1 --mem=5G --time=1-00:00 -J 04_filter_deviant_{} -o  04_filter_deviant_{}_%j.log /bin/sh 01_scripts/04.1_filter_based_on_lr.sh {} &

REGION_NUM=$1
PVAL_THRESHOLD="0.05"

### Convert ngsparalog output in list of canonical and deviant SNPs based on p-value threshold
Rscript 01_scripts/Rscripts/convert_ngsparalog_to_sitelist.R 04_ngsParalog/ngsParalog_"$REGION_NUM" "$REGION_NUM" $PVAL_THRESHOLD
