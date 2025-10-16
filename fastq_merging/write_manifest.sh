#!/bin/bash
#SBATCH --partition=GENOA

sample_id="$1"
type="$2"


echo "Running on:"
echo "  SPN: $spn"

if [[ "$type" == "normal" ]]; then
    sample_prefix="${sample_id}_normal"
elif [[ "$type" == "tumour" ]]; then
    sample_prefix="${sample_id}"
else
    echo "Error: Unknown type '$type'. Exiting."
    exit 1
fi

ena_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/ENA_Submission/"
manifest_outdir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/ENA_Submission/Manifests/"
python3 make_manifest.py --input ${ena_dir} --output ${manifest_outdir}/manifest_${sample_prefix}.json --sample "${sample_prefix}"
