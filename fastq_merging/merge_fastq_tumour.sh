#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200GB
#SBATCH --time=05:00:00
#SBATCH --job-name=merge_fastq_tumour
#SBATCH --output=merge_fastq_tumour_%A_%a.out
#SBATCH --error=merge_fastq_tumour_%A_%a.err
#SBATCH --array=1-2  # 2 SPNs × 3 purities × 2 reads = 12 combos


spns=(SPN01)
purities=(0.3)
reads=(R1 R2)

task_id=$((SLURM_ARRAY_TASK_ID - 1))
n_spns=${#spns[@]}
n_purities=${#purities[@]}
n_reads=${#reads[@]}

combinations_per_spn=$((n_purities * n_reads))

spn_idx=$((task_id / combinations_per_spn))
rem=$((task_id % combinations_per_spn))
purity_idx=$((rem / n_reads))
read_idx=$((rem % n_reads))

spn=${spns[$spn_idx]}
purity=${purities[$purity_idx]}
read=${reads[$read_idx]}

echo "Running on:"
echo "  SPN: $spn"
echo "  Purity: $purity"
echo "  Read: $read"

fastq_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/sequencing/tumour/purity_${purity}/FASTQ"
output_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/sequencing/tumour/purity_${purity}/MERGED_FASTQ"


samples=$(for i in $(ls ${fastq_dir}/*${read}*); do basename $i| cut -f 2,3 -d "_" | cut -f 1,2 -d "." ; done | sort -u)
echo ${samples[@]}

for s in $samples; do
  echo "Merging sample $s (${read})..."
  input_files="$fastq_dir/t*${s}.${read}.fastq.gz"
  output_file="$output_dir/${s}.${read}.fastq.gz"

  echo "mkdir -p ${output_dir}"
  mkdir -p ${output_dir}

  echo "pigz -dc -p ${SLURM_CPUS_PER_TASK:-1} ${input_files} | pigz -p ${SLURM_CPUS_PER_TASK:-1} > ${output_file}"
  time pigz -dc -p ${SLURM_CPUS_PER_TASK:-1} ${input_files} | pigz -p ${SLURM_CPUS_PER_TASK:-1} > ${output_file}

done

echo "✅ Done for ${spn} purity=${purity} read=${read}"
