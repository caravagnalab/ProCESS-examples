#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20gb
#SBATCH --time=01:00:00
#SBATCH --job-name=assign_process
#SBATCH --output=logs/assign_process_%A_%a.out
#SBATCH --error=logs/assign_process_%A_%a.err
#SBATCH --array=1-6
# Adjust based on spn number

module load R/4.4.1

# Define your inputs
spns=(SPN01 SPN02 SPN03 SPN04 SPN06 SPN07)

# Compute indices
task_id=$((SLURM_ARRAY_TASK_ID - 1))
n_spns=${#spns[@]}

# Determine indices
spn_idx=$((task_id))

# Get actual values
spn=${spns[$spn_idx]}

# Base path
path="/orfeo/scratch/area/lvaleriani/races/ProCESS-examples"

echo $spn

Rscript ${path}/assign_signature/assign_process.R --spn_id "${spn}"
