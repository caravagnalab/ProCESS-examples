#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-11
#SBATCH --mem 50gb
#SBATCH --time=2:00:00
#SBATCH --output=map_muts_%A_%a.out
#SBATCH --error=map_muts_%A_%a.err

# to change
# path to a temporary directory 
dir="/orfeo/fast/area/lvaleriani/tmp_nf"
mkdir -p $dir
# path to ProCESS-examples repo
base="/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/report"
# spn name
spn=SPN03

cov_values=(50 100 150 200)
pur_values=(0.3 0.6 0.9)

module load singularity/3.10.4
export TMPDIR=$dir
export SINGULARITY_TMPDIR=$dir
export TMP=$dir
image="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/process_1.1.0.sif"

input="/orfeo/cephfs/scratch/cdslab/shared/SCOUT"

cov_idx=$(( SLURM_ARRAY_TASK_ID / 3 ))
pur_idx=$(( SLURM_ARRAY_TASK_ID % 3 ))

cov=${cov_values[$cov_idx]}
pur=${pur_values[$pur_idx]}

echo "Running job with COV=$cov and PUR=$pur"

singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/map_cna_to_muts.R -i $input -s $spn -c $cov -p $pur 
