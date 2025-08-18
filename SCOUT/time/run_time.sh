#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 50gb
#SBATCH --time=10:00:00
#SBATCH --output=ProCESS.out
#SBATCH --error=ProCESS.err

module load singularity
image="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/process_1.1.0.sif"
base="/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/SCOUT"

for spn in 1 2 3 4 6 7
do
  echo $spn
  singularity exec --bind /orfeo:/orfeo --no-home $image Rscript ${base}/time/SPN0${spn}_simulate_tissue.R
  singularity exec --bind /orfeo:/orfeo --no-home $image Rscript ${base}/time/SPN0${spn}_simulate_mutation.R
done
