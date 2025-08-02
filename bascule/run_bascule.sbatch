#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem 20gb
#SBATCH --time=01:00:00
#SBATCH --job-name=bascule
#SBATCH --output=logs/bascule_%A_%a.out
#SBATCH --error=logs/bascule_%A_%a.err
#SBATCH --array=1-3

module load R/4.4.1

# change them accordingly
user="cdslab"
spns=(SPN04)
spn=${spns[$SLURM_ARRAY_TASK_ID-1]}
purities=(0.9 0.6 0.3)
purity=${purities[$SLURM_ARRAY_TASK_ID-1]}
# change with your own absolute path
path="/orfeo/scratch/cdslab/ggandolfi/Github/ProCESS-examples/"

sg cdslab << EOF
umask 007
Rscript ${path}/bascule/bascule.R --spn_id "${spn}" --coverage 50 --purity "$purity" --cna_caller ascat --vcf_caller mutect2
EOF
