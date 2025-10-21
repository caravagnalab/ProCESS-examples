#!/bin/bash
#SBATCH --partition=GENOA

sample_id="$1"
read="$2"
sample_type="$3"
purity="$4"
spn="$5"

echo "Running on:"
echo "  SPN: $spn"
echo "  Read: $read"
echo "  Type: $sample_type"
echo "  Purity: $purity"

fastq_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/sequencing/${sample_type}/purity_${purity}/FASTQ"
output_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/ENA_Submission"

if [[ "$sample_type" == "normal" ]]; then
    input_files="$fastq_dir/n*_normal_sample.${read}.fastq.gz"
    output_file="$output_dir/${spn}_normal.${read}.fastq.gz"
elif [[ "$sample_type" == "tumour" ]]; then
    input_files="$fastq_dir/t*_${sample_id}.${read}.fastq.gz"
    output_file="$output_dir/${sample_id}_${purity}.${read}.fastq.gz"
else
    echo "Error: Unknown type '$type'. Exiting."
    exit 1
fi

if [ ! -f ${output_file} ]; then
    echo "pigz -dc -p ${SLURM_CPUS_PER_TASK:-1} ${input_files} | pigz -p ${SLURM_CPUS_PER_TASK:-1} > ${output_file}"
    #time pigz -dc -p ${SLURM_CPUS_PER_TASK:-1} ${input_files} | pigz -p ${SLURM_CPUS_PER_TASK:-1} > ${output_file}
else 
    echo "File already exists"
fi
