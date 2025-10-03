#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-756
#SBATCH --job-name=subclonal
#SBATCH --mem=150GB
#SBATCH --partition THIN
#SBATCH --output=./logs/log-%A-%a

module load R/4.4.1

Rscript /orfeo/cephfs/scratch/cdslab/ebusca00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/generate_tables.R FALSE $SLURM_ARRAY_TASK_ID

