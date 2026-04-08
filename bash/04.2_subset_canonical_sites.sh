#!/bin/bash

#parallel -a ./02_info/regions.txt -j 10 srun -p small -c 1 --mem=4G --time=1-00:00 -J 08_dovcf{} -o  08_dovcf_{}_%j.log /bin/sh 01_scripts/04.2_subset_canonical_sites.sh {} &


REGION=$1
SINGLETON="04_ngsParalog/canonical_${REGION}"
SITES="02_info/sites_all_maf0.01pctind0.8_maxdepth10_${REGION}"
OUT="04_ngsParalog/canonical_sites_${REGION}"

python3 01_scripts/subset_canonical_sites.py $SINGLETON $SITES $OUT


