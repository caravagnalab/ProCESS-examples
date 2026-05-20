#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --job-name=seq_mut_copy
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-1

set -euo pipefail

base_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
spn="SPN01"
dst_root="${base_dir}/zenodo"

# only purities vary
purities=(0.6 0.9)

purity=${purities[$SLURM_ARRAY_TASK_ID]}

src="${base_dir}/${spn}/sequencing/tumour/purity_${purity}/data/"
dst="${dst_root}/${spn}/sequencing/tumour/purity_${purity}/data/"

mkdir -p "${dst}"

# copy ONLY merged mutation files
rsync -av \
    --include="mutations/seq_results_muts_merged_coverage_*.rds" \
    --exclude="mutations/seq_results_muts_SPN01_*.rds" \
    --exclude="parameters/" \
    --exclude="resources/" \
    "${src}" \
    "${dst}"

echo "Copied merged mutation files for purity ${purity}"
