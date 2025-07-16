#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --mem=30GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=out/parsing_chr%a.log
#SBATCH --job-name=parsing
#SBATCH --array=1-24

# Load R module
#source /u/cdslab/gsantacatterina/miniconda3/etc/profile.d/conda.sh
conda activate process


SPN='SPN04'
COV='100'
PUR='0.6'

# Map SLURM_ARRAY_TASK_ID to chromosome names
if [[ ${SLURM_ARRAY_TASK_ID} -le 22 ]]; then
  CHROMOSOME=${SLURM_ARRAY_TASK_ID}
elif [[ ${SLURM_ARRAY_TASK_ID} -eq 23 ]]; then
  CHROMOSOME="X"
elif [[ ${SLURM_ARRAY_TASK_ID} -eq 24 ]]; then
  CHROMOSOME="Y"
fi

echo "Processing chromosome ${CHROMOSOME}"

# Run R script with the current chromosome
Rscript preprocess.R --spn_id $SPN --chr ${CHROMOSOME} --purity $PUR --coverage $COV
