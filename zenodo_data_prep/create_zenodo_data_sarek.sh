#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --job-name=sarek_copy
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-0

set -euo pipefail

base_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
spn="SPN01"
dst_root="${base_dir}/zenodo"

# Parameter grids
coverages=(50)
purities=(0.3)

# Convert array index → parameter pair
n_purities=${#purities[@]}

coverage_idx=$(( SLURM_ARRAY_TASK_ID / n_purities ))
purity_idx=$(( SLURM_ARRAY_TASK_ID % n_purities ))

coverage=${coverages[$coverage_idx]}
purity=${purities[$purity_idx]}

src="${base_dir}/${spn}/sarek/${coverage}x_${purity}p"
dst="${dst_root}/${spn}/sarek/${coverage}x_${purity}p"

mkdir -p "${dst}"

if [[ -d "${src}/variant_calling" ]]; then

    echo "Copying variant_calling for ${coverage}x_${purity}p"

    rsync -av \
        --exclude="multiqc/" \
	--exclude="preprocessing/" \
	--exclude="reports/" \
	--exclude="pipeline_info/" \
	--exclude="csv/" \
        "${src}/variant_calling/" \
        "${dst}/variant_calling/"

    echo "Done: ${coverage}x_${purity}p"

else
    echo "Missing ${src}/variant_calling"
fi
