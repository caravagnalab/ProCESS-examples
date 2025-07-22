#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 100gb
#SBATCH --time=96:00:00
#SBATCH --output=seq.out
#SBATCH --error=seq.err

module load singularity/3.10.4

# change them accordingly
user="cdslab"
spn="SPN07"

# change with your own absolute path
path="/orfeo/scratch/cdslab/antonelloa/ProCESS-examples/build_cohorts"

# keep them as they are
partition=THIN
phylo="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/process/phylo_forest.sff"
tmp="/orfeo/cephfs/fast/cdslab/${USER}/tmp_files"
image="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/process_1.0.2.sif"
out="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/sequencing"
sarek_output_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/sarek"
tumourevo_output_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/tumourevo"

$path/benchmark_build_cohort.py -P $partition -A $user -s $tmp -I $image $spn $phylo $out -C $path/orfeo.config -SD $sarek_output_dir -TD $tumourevo_output_dir
