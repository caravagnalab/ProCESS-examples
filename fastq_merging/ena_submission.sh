#!/bin/bash

# -----------------------------
# Default values
# -----------------------------
sample_id="$1"
sample_type="$2"

# -----------------------------
# Submit Job 1 based on type
# -----------------------------
job_ids_merging=()
if [[ "$sample_type" == "normal" ]]; then
    purity=1
    JOB1_R1_ID=$(sbatch --parsable --job-name=merging_fastq_${sample_id}_${sample_type}_R1 --nodes=1 \
	    --cpus-per-task=12 --mem=50GB --time=04:00:00 \
	    --output=out/merging_fastq_${sample_id}_${sample_type}_R1.out --error=out/merging_fastq_${sample_id}_${sample_type}_R1.err merge_fastq_normal.sh "$sample_id" "R1" "$sample_type" "$purity")
    JOB1_R2_ID=$(sbatch --parsable --job-name=merging_fastq_${sample_id}_${sample_type}_R2 --nodes=1 \
	    --cpus-per-task=12 --mem=50GB --time=04:00:00 \
	    --output=out/merging_fastq_${sample_id}_${sample_type}_R2.out --error=out/merging_fastq_${sample_id}_${sample_type}_R2.err merge_fastq_normal.sh "$sample_id" "R2" "$sample_type" "$purity")
    echo "Merging normal sample for $sample_id in R1: $JOB1_R1_ID"
    echo "Merging normal sample for $sample_id in R2: $JOB1_R2_ID"
elif [[ "$sample_type" == "tumour" ]]; then
    purities=(0.3 0.6 0.9)
    for purity in ${purities[@]}; do
      JOB1_R1_ID=$(sbatch --parsable --job-name=merging_fastq_${sample_id}_${sample_type}_R1 --nodes=1 \
	      --cpus-per-task=24 --mem=20GB --time=04:00:00 \
              --output=out/merging_fastq_${sample_id}_${sample_type}_${purity}_R1.out --error=out/merging_fastq_${sample_id}_${sample_type}_${purity}_R1.err merge_fastq_normal.sh \
              "$sample_id" "R1" "$sample_type" "$purity")
      JOB1_R2_ID=$(sbatch --parsable --job-name=merging_fastq_${sample_id}_${sample_type}_R2 --nodes=1 \
	      --cpus-per-task=24 --mem=20GB --time=04:00:00 \
              --output=out/merging_fastq_${sample_id}_${sample_type}_${purity}_R2.out --error=out/merging_fastq_${sample_id}_${sample_type}_${purity}_R2.err merge_fastq_normal.sh \
              "$sample_id" "R2" "$sample_type" "$purity")
      echo "Merging sample for $sample_id in R1 for purity $purity: $JOB1_R1_ID"
      echo "Merging sample for $sample_id in R2 for purity $purity: $JOB1_R2_ID"
      job_ids_merging+=("$JOB1_R1_ID","$JOB1_R2_ID")
    done
else
    echo "Error: Unknown type '$type'. Exiting."
    exit 1
fi
# -----------------------------
# Submit Job 2 with dependency on Job 1
# -----------------------------
jobs_to_wait=$(IFS=, ; echo "${job_ids_merging[*]}")
echo $jobs_to_wait
JOB2_ID=$(sbatch --parsable --dependency=afterok:$jobs_to_wait --nodes=1 \
               --cpus-per-task=1 \
               --mem=2GB \
               --job-name=generating_manifest_${sample_id} \
               --output=out/generating_manifest_${sample_id}.out \
               --error=out/generating_manifest_${sample_id}.err \
               --time=00:10:00 \
               write_manifest.sh "$sample_id" "$sample_type")
echo "Job 2 submitted with ID: $JOB2_ID (will start after Job 1 finishes successfully)"

## --------------------------------------------------------
## Submit Job 3 with dependency on Job 2 - write checksums
## --------------------------------------------------------
#
JOB3_ID=$(sbatch --parsable --dependency=afterok:$jobs_to_wait --nodes=1 \
              --cpus-per-task=4 \
              --mem=40GB \
              --job-name=write_checksum_${sample_id} \
              --output=out/write_checksum_${sample_id}.out \
              --error=out/write_checksum_${sample_id}.err \
              --time=03:00:00 \
              write_checksum.sh "$sample_id" "$sample_type")
echo "Job 2 submitted with ID: $JOB2_ID (will start after Job 1 finishes successfully)"
