#!/bin/bash
#SBATCH --partition=GENOA


sample_id="$1"
sample_type="$2"


echo "Running on:"
echo "  SPN: $spn"

if [[ "$sample_type" == "normal" ]]; then
    sample_prefix="${sample_id}_normal"
elif [[ "$sample_type" == "tumour" ]]; then
    sample_prefix="${sample_id}"
else
    echo "Error: Unknown type '$type'. Exiting."
    exit 1
fi

ena_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/ENA_Submission"

md5sum $ena_dir/$sample_prefix*gz > $ena_dir/checksum_${sample_prefix}.tsv
