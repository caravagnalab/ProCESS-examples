#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20gb
#SBATCH --time=01:00:00
#SBATCH --job-name=bcftools_concat
#SBATCH --output=logs/bcftools_%A_%a.out
#SBATCH --error=logs/bcftools_%A_%a.err
#SBATCH --array=1-27
# Adjust based on total combinations: spns × purities × coverages

module load bcftools

# Define your inputs
spns=(SPN03 SPN06 SPN07)
tool=(ascat sequenza)
purities=(0.3 0.6 0.9)
coverages=(50 100 150)  

# Compute indices
task_id=$((SLURM_ARRAY_TASK_ID - 1))
n_spns=${#spns[@]}
n_purities=${#purities[@]}
n_coverages=${#coverages[@]}

# Total combinations per SPN
combinations_per_spn=$((n_purities * n_coverages))

# Determine indices
spn_idx=$((task_id / combinations_per_spn))
rem=$((task_id % combinations_per_spn))
purity_idx=$((rem / n_coverages))
coverage_idx=$((rem % n_coverages))

# Get actual values
spn=${spns[$spn_idx]}
purity=${purities[$purity_idx]}
coverage=${coverages[$coverage_idx]}
sarek_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/sarek/${coverage}x_${purity}p/variant_calling/strelka/"
tumourevo_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/sequencing/tumourevo/"

for i in $(ls ${sarek_dir} | grep vs | cut -f 1-2 -d "_"); do 
    indel_vcf=$sarek_dir/$i\_vs_normal_sample/$i\_vs_normal_sample.strelka.somatic_indels.vcf.gz
    snv_vcf=$sarek_dir/$i\_vs_normal_sample/$i\_vs_normal_sample.strelka.somatic_snvs.vcf.gz
    final_vcf=$sarek_dir/$i\_vs_normal_sample/$i\_vs_normal_sample.strelka.somatic.vcf.gz
    final_tbi=$sarek_dir/$i\_vs_normal_sample/$i\_vs_normal_sample.strelka.somatic.vcf.gz.tbi
    echo "bcftools concat $indel_vcf $snv_vcf -a -Oz -o $final_vcf"
    bcftools concat $indel_vcf $snv_vcf -a -Oz -o $final_vcf
    echo "bcftools index -t $final_vcf -o $final_tbi"
    bcftools index -t $final_vcf -o $final_tbi
    
    for t in ${tool[@]}; do
        tumourevo_csv=${tumourevo_dir}/tumourevo_${coverage}x_${purity}p_strelka_${t}.csv
        sed -i 's/_snvs//g' ${tumourevo_csv}
    done
done





