#!/bin/bash
#SBATCH -J "het_by_indv"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p ibis_small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=1-00:00
#SBATCH --mem=10G

#cat 02_info/bam.filelist.valeria | parallel -j40 srun -c 1 --mem 10G -p ibis_small -o log_%j --time 1-00:00 ./01_scripts/het_by_indv.sh {}

###this script will work on each bamfiles and calculate individual saf and SFS for heterozygosity
#maybe edit
NB_CPU=1 #change accordingly in SLURM header

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
module load angsd
ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
mkdir het_by_indv

BAM=$1
ID=$(grep $BAM 02_info/bam.filelist.valeria | grep -o [A-Z][A-Z][A-Z]s_[0-9][0-9][0-9]-[0-9][0-9])

# Do SAF for one sample

echo $ID
		
	angsd -P $NB_CPU -underFlowProtect 1 \
	-dosaf 1 -GL 2 -doMajorMinor 5 -doCounts 1 \
	-anc ../../ancestral_genome_maskdev_genes1kb.fasta \
    -rf 02_info/regions.txt \
    -uniqueOnly 1 -only_proper_pairs 1 \
	-remove_bads 1 -minMapQ 30 -minQ 20 \
    -minInd 1  -setMinDepthInd 1 \
	-i $BAM -out het_by_indv/"$ID"
	
	#the output if a saf
	
	#now we use winSFS to calculate a 1dimension sfs
    
    ## XAVIER : I now use winsfs. It's faster while still using all sites

	echo "estimate winsfs for indv $ID"
    
    winsfs shuffle -v --output het_by_indv/"$ID".saf.shuf \
    het_by_indv/"$ID".saf.idx
    
    winsfs -v het_by_indv/"$ID".saf.shuf > het_by_indv/"$ID".dsfs
	
    rm het_by_indv/"$ID".saf.shuf
    
	#Remove header
    tail -n 1 het_by_indv/"$ID".dsfs > het_by_indv/"$ID".dsfs.tmp
    mv het_by_indv/"$ID".dsfs.tmp het_by_indv/"$ID".dsfs
    
    ### Extract heterozygosity 
    HET=$(awk '{ sum += $2; total += $0 } END { print sum / total }' het_by_indv/"$ID".dsfs)
    echo -e $ID"\t"$HET >> het_by_indv/summary.het