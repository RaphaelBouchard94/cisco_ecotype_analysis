#!/bin/bash
cd $SLURM_SUBMIT_DIR

#srun -p medium -c 5 --mem=5G --time=7-00:00 -J 05_ngsLD_chr34 -o 05_ngsLD_per_chr34.log /bin/sh 01_scripts/14_ngsLD_parallel.sh CM078251.1 &


#To run
#parallel -a ./02_info/regions.txt -j 4 srun -p medium -c 5 --mem=20G --time=7-00:00 -J 05_ngsLD_{} -o 05_ngsLD_per_chr_{}_%j.log /bin/sh 01_scripts/05_run_ngsLD_parallel.sh {} &

#Important ngsLD parameters:

#--probs: specification of whether the input is genotype probabilities (likelihoods or posteriors)?
#--n_ind INT: sample size (number of individuals).
#--n_sites INT: total number of sites.
#--max_kb_dist DOUBLE: maximum distance between SNPs (in Kb) to calculate LD. Set to 0(zero) to disable filter. [100]
#--max_snp_dist INT: maximum distance between SNPs (in number of SNPs) to calculate LD. Set to 0 (zero) to disable filter. [0]
#--rnd_sample= sample x% of comparison between snp
#--n_threads INT: number of threads to use. [1]
#--out FILE: output file name. [stdout]
#
REGION="$1"
NGSLD=./01_scripts/ngsLD
NIND=$(wc -l 02_info/bam.filelist)
NCPU=5

#Prepare pos file

#awk '{print $1"\t"$2}' 04_ngsParalog/canonical_"$REGION" > 05_ngsLD/canonical_"$REGION".pos

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
N_IND=$(wc -l 02_info/bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*}
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

NSITES=$(wc -l 14_ngsld/canonical_"$REGION".pos)

$NGSLD/ngsLD --geno 04_ngsParalog/all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$REGION"_canonical.beagle.gz \
	--probs \
	--pos 14_ngsld/canonical_"$REGION".pos \
	--n_ind "$NIND" \
	--n_sites "$NSITES" \
	--min_maf 0.01 \
	--max_kb_dist 0 \
	--rnd_sample 0.5 \
	--n_threads "$NCPU" \
	--out 14_ngsld/all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$REGION".ld 
