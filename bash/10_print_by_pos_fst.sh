#!/bin/bash
#SBATCH -J "10_print_by_pos_fst"
#SBATCH -o log_%j
#SBATCH -c 10
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.bouchard.3@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=8G

# Configuration
NB_CPU=10
POP1="${1}"
POP2="${2}"

# Validate inputs
if [[ -z "${POP1}" ]] || [[ -z "${POP2}" ]]; then
    echo "Error: Missing required parameters" >&2
    echo "Usage: $0 <pop1> <pop2>" >&2
    exit 1
fi

cd "${SLURM_SUBMIT_DIR}" || exit 1
module load angsd/0.931

source 01_scripts/01_config.sh || exit 1

realSFS fst print 11_fst/"$POP1"_"$POP2"_genome_wide.fst.idx -P $NB_CPU > 11_fst/"$POP1"_"$POP2".bypos.sfs

