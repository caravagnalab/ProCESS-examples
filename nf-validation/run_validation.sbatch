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


config=orfeo.config
profile=$(sinfo -h -o "%P %a %D %t" | grep -w 'EPYC\|GENOA\|THIN' |awk '$2 == "up" && $4 ~ /idle|mix/ {print tolower($1), $3}' | awk '{sum[$1] += $2} END {for (p in sum) print p, sum[p]}' | sort -k2 -nr | head -n1 | cut -f1 -d " ")
#/orfeo/cephfs/scratch/cdslab/shared/SCOUT/nextflow run nf-core/sarek -r 3.5.1 --genome GATK.GRCh38 --input $input     --step variant_calling --tools cnvkit,freebayes,strelka,ascat,mutect2 --joint_mutect2 true     --outdir $output_dir_combination -profile singularity,${profile} -c $config

nextflow run main.nf --input samplesheet.csv --outdir /orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/validate_nextflow/results/ -profile $profile -c $config
