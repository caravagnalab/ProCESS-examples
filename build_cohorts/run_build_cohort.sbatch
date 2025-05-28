#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 20gb
#SBATCH --time=96:00:00
#SBATCH --output=seq.out
#SBATCH --error=seq.err

module load singularity/3.10.4

user="cdslab"
spn="SPN06"

path="/orfeo/cephfs/scratch/cdslab/ahaghighi/packages/ProCESS-examples/build_cohorts"

partition=EPYC
phylo="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/process/phylo_forest.sff"
tmp="/orfeo/cephfs/fast/cdslab/${USER}/tmp_files"
image="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/process_1.0.2.sif"
out="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/sequencing"
sarek_output_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/sarek"
tumourevo_output_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/tumourevo"

$path/benchmark_build_cohort.py -P $partition -A $user -s $tmp -I $image $spn $phylo $out -C $path/orfeo.config -SD $sarek_output_dir -TD $tumourevo_output_dir
