#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --job-name=sarek_normal_copy
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-0

set -euo pipefail

base_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
dst_root="${base_dir}/zenodo"

# SPN list
spns=(
    SPN01
)

# Select SPN from array index
spn=${spns[$SLURM_ARRAY_TASK_ID]}

src="${base_dir}/${spn}/sarek/normal"
dst="${dst_root}/${spn}/sarek/normal"

mkdir -p "${dst}"

echo "Copying sarek for ${spn}"

rsync -av \
    --exclude="multiqc/" \
    --exclude="preprocessing/" \
    --exclude="reports/" \
    --exclude="pipeline_info/" \
    --exclude="csv/" \
    "${src}" \
    "${dst}"

echo "Done: ${spn}"
