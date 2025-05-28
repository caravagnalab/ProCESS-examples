#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 50gb
#SBATCH --time=1:00:00
#SBATCH --output=ProCESS.out
#SBATCH --error=ProCESS.err

module load singularity
image="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/process_1.0.2.sif"

base="/orfeo/cephfs/scratch/cdslab/ahaghighi/packages/ProCESS-examples/SCOUT/SPN06"

singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/simulate_tissue.R
singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/simulate_mutation.R

