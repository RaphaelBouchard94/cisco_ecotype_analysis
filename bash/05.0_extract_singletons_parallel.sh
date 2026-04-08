#!/bin/bash

#To call:
#parallel -a 02_info/regions.txt  -j 30 srun -p small -c 1 --mem=5G --time=1:00:00 -J extract_wanted_snps_{} -o extract_wanted_snps_{}_%j.log /bin/sh 01_scripts/05.0_extract_singletons_parallel.sh {} &

REGION="$1"

python 01_scripts/Scripts/beagle_extract_wanted_snps.py 03_saf_maf_gl_all/all_maf0.01pctind0.8_maxdepth10_"${REGION}".beagle.gz 04_ngsParalog/canonical_"${REGION}" 04_ngsParalog/all_maf0.01pctind0.8_maxdepth10_"${REGION}"_canonical.beagle
