#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --job-name=zenodo_pack
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-0

set -euo pipefail

base_root="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/zenodo"

# List of SPNs
spns=(
    SPN01
)

# Select SPN from array index
spn=${spns[$SLURM_ARRAY_TASK_ID]}

base="${base_root}/${spn}"

# pigz threads
THREADS=${SLURM_CPUS_PER_TASK:-4}

echo "Processing ${spn}"
echo "Using ${THREADS} threads"

cd "${base}"
#
## -----------------------------
## tumourevo
## -----------------------------
#if [[ -d "tumourevo" ]]; then
#
#    echo "Compressing tumourevo"
#
#    tar -I "pigz -p ${THREADS}" \
#        -cf tumourevo.tar.gz \
#        tumourevo/
#
#    echo "Created tumourevo.tar.gz"
#
#else
#    echo "Missing tumourevo/"
#fi
#
## -----------------------------
## sarek
## -----------------------------
#if [[ -d "sarek" ]]; then
#
#    echo "Compressing sarek"
#
#    tar -I "pigz -p ${THREADS}" \
#        -cf sarek.tar.gz \
#        sarek/
#
#    echo "Created sarek.tar.gz"
#
#else
#    echo "Missing sarek/"
#fi

# -----------------------------
# sequencing
# -----------------------------
if [[ -d "sequencing" ]]; then

    echo "Compressing sequencing"

    tar -I "pigz -p ${THREADS}" \
        -cf sequencing.tar.gz \
        sequencing/

    echo "Created sequencing.tar.gz"

else
    echo "Missing sequencing/"
fi

# -----------------------------
# checksums
# -----------------------------
#sha256sum *.tar.gz > CHECKSUMS.txt

#echo "Finished ${spn}"
