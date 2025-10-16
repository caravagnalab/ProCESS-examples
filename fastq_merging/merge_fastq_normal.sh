#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=80GB
#SBATCH --time=05:00:00
#SBATCH --job-name=merge_fastq_normal
#SBATCH --output=merge_fastq_normal_%A_%a.out
#SBATCH --error=merge_fastq_normal_%A_%a.err
#SBATCH --array=1-2  # 2 SPNs ×  2 reads = 6 combos


spns=(SPN01)
reads=(R1 R2)

task_id=$((SLURM_ARRAY_TASK_ID - 1))
n_spns=${#spns[@]}
n_reads=${#reads[@]}

combinations_per_spn=$((n_reads))

spn_idx=$((task_id / combinations_per_spn))
rem=$((task_id % combinations_per_spn))
read_idx=$((rem % n_reads))

spn=${spns[$spn_idx]}
read=${reads[$read_idx]}

echo "Running on:"
echo "  SPN: $spn"
echo "  Read: $read"

fastq_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/sequencing/normal/purity_1/FASTQ"
output_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/sequencing/normal/purity_1/MERGED_FASTQ"

input_files="$fastq_dir/n*_normal_sample.${read}.fastq.gz"
output_file="$output_dir/normal_sample.${read}.fastq.gz"

echo "mkdir -p ${output_dir}"
mkdir -p ${output_dir}

echo "pigz -dc -p ${SLURM_CPUS_PER_TASK:-1} ${input_files} | pigz -p ${SLURM_CPUS_PER_TASK:-1} > ${output_file}"
time pigz -dc -p ${SLURM_CPUS_PER_TASK:-1} ${input_files} | pigz -p ${SLURM_CPUS_PER_TASK:-1} > ${output_file}
