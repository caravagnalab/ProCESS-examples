#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH -A cdslab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-63
#SBATCH --job-name=tables_subclonal
#SBATCH --mem=100GB
#SBATCH --partition EPYC
#SBATCH --output=./logs/log-%A-%a

module load R/4.4.1
export R_LIBS_USER=/u/cdslab/erivar00/new_r_env
Rscript /orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/generate_tables.R FALSE $SLURM_ARRAY_TASK_ID

