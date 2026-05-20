#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --job-name=tumourevo_copy
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-7 #coverages*purities*vcf_caller*cna_caller

set -euo pipefail

base_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
spn="SPN01"
dst_root="${base_dir}/zenodo"

# Parameter grids
coverages=(150)
purities=(0.6 0.9)
cna_callers=(ascat sequenza)
vcf_callers=(mutect2 strelka)

n_purities=${#purities[@]}
n_cna=${#cna_callers[@]}
n_vcf=${#vcf_callers[@]}

i=$SLURM_ARRAY_TASK_ID

cov_idx=$(( i / (n_purities * n_cna * n_vcf) ))
rem=$(( i % (n_purities * n_cna * n_vcf) ))

pur_idx=$(( rem / (n_cna * n_vcf) ))
rem=$(( rem % (n_cna * n_vcf) ))

cna_idx=$(( rem / n_vcf ))
vcf_idx=$(( rem % n_vcf ))

coverage=${coverages[$cov_idx]}
purity=${purities[$pur_idx]}
cna=${cna_callers[$cna_idx]}
vcf=${vcf_callers[$vcf_idx]}

src="${base_dir}/${spn}/tumourevo/${coverage}x_${purity}p_${vcf}_${cna}/"
dst="${dst_root}/${spn}/tumourevo/${coverage}x_${purity}p_${vcf}_${cna}/"

mkdir -p "${dst}"

# COPY EVERYTHING EXCEPT SELECTED FOLDERS
rsync -av \
  --exclude="lifter/" \
  --exclude="bcftools/" \
  --exclude="pipeline_info/" \
  "${src}" \
  "${dst}"

echo "Copied ${coverage}x_${purity}p_${vcf}_${cna} (filtered)"
