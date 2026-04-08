#!/bin/bash

cd 03_saf_maf_gl_all/

for i in $(ls *mafs); do cut -f1-4 $i | tail -n+2 > ../02_info/sites_${i%.mafs}; done

cd ../02_info/

for i in $(ls sites*) ; do angsd sites index $i ; done
