#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=validation
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=96:00:00
#SBATCH --output=validation.out 
#SBATCH --error=validation.err
#SBATCH -A cdslab

module load java
module load singularity
module load R/4.4.1

nf=/orfeo/scratch/area/lvaleriani/
base="/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/nf-validation"
results=$base/results
config=/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/build_cohorts/orfeo.config
profile=epyc

$nf/nextflow run $base/main.nf --input $base/samplesheet.csv --outdir $results -profile singularity,$profile -c $config
