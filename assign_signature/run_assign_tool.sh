#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20gb
#SBATCH --time=01:00:00
#SBATCH --job-name=assign
#SBATCH --output=logs/assign_%A_%a.out
#SBATCH --error=logs/assign_%A_%a.err
#SBATCH --array=1-9
# Adjust based on total combinations: spns × purities × coverages

module load R/4.4.1

# Define your inputs
#spns=(SPN01 SPN02 SPN03 SPN04 SPN05 SPN06 SPN07)
spns=(SPN07)
purities=(0.9 0.6 0.3)
coverages=(50 100 150)  

# Compute indices
task_id=$((SLURM_ARRAY_TASK_ID - 1))
n_spns=${#spns[@]}
n_purities=${#purities[@]}
n_coverages=${#coverages[@]}

# Total combinations per SPN
combinations_per_spn=$((n_purities * n_coverages))

# Determine indices
spn_idx=$((task_id / combinations_per_spn))
rem=$((task_id % combinations_per_spn))
purity_idx=$((rem / n_coverages))
coverage_idx=$((rem % n_coverages))

# Get actual values
spn=${spns[$spn_idx]}
purity=${purities[$purity_idx]}
coverage=${coverages[$coverage_idx]}

# Base path
path="/orfeo/scratch/area/lvaleriani/races/ProCESS-examples"

echo $spn
echo $purity
echo $coverage

#Rscript ${path}/assign_signature/assign_tool.R --spn_id "${spn}" --coverage "$coverage" --purity "$purity" --signature "SigProfiler"
#Rscript ${path}/assign_signature/assign_tool.R --spn_id "${spn}" --coverage "$coverage" --purity "$purity" --signature "BASCULE"

Rscript ${path}/assign_signature/assign_tool_mobster.R --spn_id "${spn}" --coverage "$coverage" --purity "$purity" --signature "SigProfiler"
Rscript ${path}/assign_signature/assign_tool_mobster.R --spn_id "${spn}" --coverage "$coverage" --purity "$purity" --signature "BASCULE"


