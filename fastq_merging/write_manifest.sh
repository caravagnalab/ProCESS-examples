#!/bin/bash
#SBATCH --partition=GENOA

sample_id="$1"
sample_type="$2"
purity="$3"

echo "Running on:"
echo "  SPN: $spn"

if [[ "$sample_type" == "normal" ]]; then
    sample_prefix="${sample_id}_normal"
elif [[ "$sample_type" == "tumour" ]]; then
    sample_prefix="${sample_id}_${purity}"
else
    echo "Error: Unknown type '$type'. Exiting."
    exit 1
fi

ena_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/ENA_Submission/"
manifest_outdir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/ENA_Submission/Manifests/"
sub_file="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/ENA_Submission/samples-2025-10-20T12_26_18.csv"
python3 make_manifest.py --input ${ena_dir} --output ${manifest_outdir}/manifest_${sample_prefix}.json --sample "${sample_prefix}" --ena_submission_file "${sub_file}"
