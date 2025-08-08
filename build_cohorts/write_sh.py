#!/usr/bin/python3

import os
import sys
import math
import glob
import time
import subprocess
import argparse


## This part is currently run sequentially


write_tumour_type_file="""#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=tumour_type
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=01:00
#SBATCH --output=tumour_type_%J.out 
#SBATCH --error=tumour_type_%J.err
#SBATCH -A {ACCOUNT}

echo ${DIR}/../references/SCOUT_tumour_types.txt
tumour_type=$(grep ${SPN} ${DIR}/../references/SCOUT_tumour_types.txt | cut -f 3)
echo $tumour_type > "${BASEDIR}/tumour_type.txt"
"""

sarek_file_launcher="""#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=sarek_mapping_{JOB_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=96:00:00
#SBATCH --output=sarek_mapping_{JOB_NAME}_%J.out 
#SBATCH --error=sarek_mapping_{JOB_NAME}_%J.err
#SBATCH -A {ACCOUNT}

module load java
module load singularity

input_dir={INPUT_DIR}
input="${input_dir}/sarek_{JOB_NAME}.csv"

output_base_dir={SAREK_OUT}
output_dir_combination="${output_base_dir}/{JOB_NAME}"

config={CONFIG}
profile=$(sinfo -h -o "%P %a %D %t" | grep -w 'EPYC\|GENOA\|THIN' |awk '$2 == "up" && $4 ~ /idle|mix/ {print tolower($1), $3}' | awk '{sum[$1] += $2} END {for (p in sum) print p, sum[p]}' | sort -k2 -nr | head -n1 | cut -f1 -d " ")
/orfeo/cephfs/scratch/cdslab/shared/SCOUT/nextflow run nf-core/sarek -r 3.5.1 --input $input \
    --outdir $output_dir_combination -profile singularity,${profile} -c $config
    
chgrp -R cdslab $output_dir_combination
chmod -R g+wrx $output_dir_combination
"""

sarek_file_normal_launcher="""#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=sarek_mapping_vc_{JOB_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=96:00:00
#SBATCH --output=sarek_mapping_vc_{JOB_NAME}_%J.out 
#SBATCH --error=sarek_mapping_vc_{JOB_NAME}_%J.err
#SBATCH -A {ACCOUNT}

module load java
module load singularity

input_dir={INPUT_DIR}
input="${input_dir}/sarek_{JOB_NAME}.csv"

output_base_dir={SAREK_OUT}
output_dir_combination="${output_base_dir}/{JOB_NAME}"

config={CONFIG}
profile=$(sinfo -h -o "%P %a %D %t" | grep -w 'EPYC\|GENOA\|THIN' |awk '$2 == "up" && $4 ~ /idle|mix/ {print tolower($1), $3}' | awk '{sum[$1] += $2} END {for (p in sum) print p, sum[p]}' | sort -k2 -nr | head -n1 | cut -f1 -d " ")
/orfeo/cephfs/scratch/cdslab/shared/SCOUT/nextflow run nf-core/sarek -r 3.5.1 --input $input \
    --outdir $output_dir_combination --tools haplotypecaller,freebayes -profile singularity,${profile} -c $config

chgrp -R cdslab $output_dir_combination
chmod -R g+wrx $output_dir_combination
"""

sarek_variant_calling_launcher="""#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=sarek_VC_{JOB_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=96:00:00
#SBATCH --output=sarek_VC_{JOB_NAME}_%J.out 
#SBATCH --error=sarek_VC_{JOB_NAME}_%J.err
#SBATCH -A {ACCOUNT}

module load java
module load singularity

input_dir={INPUT_DIR}
input="${input_dir}/sarek_variant_calling_{JOB_NAME}.csv"

output_base_dir={SAREK_OUT}
output_dir_combination="${output_base_dir}/{JOB_NAME}"

config={CONFIG}
profile=$(sinfo -h -o "%P %a %D %t" | grep -w 'EPYC\|GENOA\|THIN' |awk '$2 == "up" && $4 ~ /idle|mix/ {print tolower($1), $3}' | awk '{sum[$1] += $2} END {for (p in sum) print p, sum[p]}' | sort -k2 -nr | head -n1 | cut -f1 -d " ")

/orfeo/cephfs/scratch/cdslab/shared/SCOUT/nextflow run nf-core/sarek -r 3.5.1 --genome GATK.GRCh38 --input $input \
    --step variant_calling --tools cnvkit,freebayes,strelka,haplotypecaller,ascat,mutect2 --joint_mutect2 true \
    --outdir $output_dir_combination -profile singularity,${profile} -c $config
    
chgrp -R cdslab $output_dir_combination
chmod -R g+wrx $output_dir_combination
"""


sequenza_launcher="""#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=sequenza_{JOB_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --output=sequenza_{JOB_NAME}_%J.out 
#SBATCH --error=sequenza_{JOB_NAME}_%J.err
#SBATCH -A {ACCOUNT}

module load java
module load singularity

input_dir={INPUT_DIR}
input="${input_dir}/sarek_variant_calling_{JOB_NAME}.csv"

output_base_dir={SAREK_OUT}
output_dir_combination="${output_base_dir}/{JOB_NAME}"

config={CONFIG}
base={PROCESS_DIR}
profile=$(sinfo -h -o "%P %a %D %t" | grep -w 'EPYC\|GENOA\|THIN' |awk '$2 == "up" && $4 ~ /idle|mix/ {print tolower($1), $3}' | awk '{sum[$1] += $2} END {for (p in sum) print p, sum[p]}' | sort -k2 -nr | head -n1 | cut -f1 -d " ")

/orfeo/cephfs/scratch/cdslab/shared/SCOUT/nextflow run $base/main.nf -profile singularity,${profile} --input $input --outdir $output_dir_combination -c $config --tool sequenza

chgrp -R cdslab $output_dir_combination
chmod -R g+wrx $output_dir_combination
"""


battenberg_launcher="""#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=battenberg_{JOB_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --output=battenberg_{JOB_NAME}_%J.out 
#SBATCH --error=battenberg_{JOB_NAME}_%J.err
#SBATCH -A {ACCOUNT}

module load java/11.0.22
module load singularity

input_dir={INPUT_DIR}
input="${input_dir}/sarek_variant_calling_{JOB_NAME}.csv"

output_base_dir={SAREK_OUT}
output_dir_combination="${output_base_dir}/{JOB_NAME}"

base={PROCESS_DIR}
config=$base/batt.config
profile=$(sinfo -h -o "%P %a %D %t" | grep -w 'EPYC\|GENOA\|THIN' |awk '$2 == "up" && $4 ~ /idle|mix/ {print tolower($1), $3}' | awk '{sum[$1] += $2} END {for (p in sum) print p, sum[p]}' | sort -k2 -nr | head -n1 | cut -f1 -d " ")

/orfeo/cephfs/scratch/cdslab/shared/SCOUT/nextflow run $base/main.nf -profile singularity,${profile} --input $input --outdir $output_dir_combination -c $config -plugins nf-schema@2.0.0 --tool battenberg

chgrp -R cdslab $output_dir_combination
chmod -R g+wrx $output_dir_combination
"""


tumourevo_launcher="""#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=tumourevo_{JOB_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=96:00:00
#SBATCH --output=tumourevo_{JOB_NAME}_%J.out 
#SBATCH --error=tumourevo_{JOB_NAME}_%J.err
#SBATCH -A {ACCOUNT}

module load java
module load singularity

input_dir={INPUT_DIR}
input="${input_dir}/tumourevo_{JOB_NAME}.csv"

output_base_dir={TUMOUREVO_OUT}
output_dir_combination="${output_base_dir}/{JOB_NAME}"

config={CONFIG}
profile=$(sinfo -h -o "%P %a %D %t" | grep -w 'EPYC\|GENOA\|THIN' |awk '$2 == "up" && $4 ~ /idle|mix/ {print tolower($1), $3}' | awk '{sum[$1] += $2} END {for (p in sum) print p, sum[p]}' | sort -k2 -nr | head -n1 | cut -f1 -d " ")

/orfeo/cephfs/scratch/cdslab/shared/SCOUT/nextflow run /orfeo/cephfs/scratch/cdslab/shared/SCOUT/tumourevo/main.nf --input $input \
    --tools mobster,viber,pyclone-vi,sparsesignatures,sigprofiler \
    --genome GRCh38 \
    --fasta /orfeo/LTS/CDSLab/LT_storage/ref_genomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
    --download_cache_vep false \
    --vep_cache /orfeo/LTS/CDSLab/LT_storage/ref_genomes/VEP \
    --vep_genome GRCh38 \
    --vep_cache_version 110 \
    --vep_species Homo_sapiens \
    --filter false \
    --outdir $output_dir_combination -profile singularity,${profile} -c $config
    
chgrp -R cdslab $output_dir_combination
chmod -R g+wrx $output_dir_combination
"""

def get_lot_prefix(seq_type):
    if seq_type=='normal':
        return 'n'
    if seq_type=='normal_with_preneoplastic':
        return 's'
    if seq_type=='tumour':
        return 't'
    raise TypeError("Only \"tumour\", \"normal\", and "
                    + "\"normal_with_preneoplastic\" are supported")

def get_sample_names_from_FASTQ(fastq_dir):
    suffix = '.R1.fastq.gz'
    fastq_files = glob.glob(f'{fastq_dir}/t*_*{suffix}')

    sample_names = set()
    for fastq_file in fastq_files:
        fastq_basename = os.path.basename(fastq_file)
        prefix_up_to = fastq_basename.find('_')
        
        sample_names.add(fastq_basename[prefix_up_to+1:-len(suffix)])
    
    return sorted(list(sample_names))
  
  
    
def write_tumourevo_lines(tumourevo_file, SPN, sample_name, combination, coverage, purity, sarek_output_dir, cancer_type):
    variant_caller = combination[0]
    cna_caller = combination[1]
    
    base_path = f'{sarek_output_dir}/{coverage}x_{purity}p'
    path = f'{base_path}/variant_calling'
    
    if variant_caller == 'mutect2':
        rel_path = f'{path}/{variant_caller}/{SPN}'
        name = f'{SPN}.mutect2.filtered.vcf.gz'
    elif variant_caller == 'strelka':
        rel_path = f'{path}/{variant_caller}/{sample_name}_vs_normal_sample'
        name= f'{sample_name}_vs_normal_sample.strelka.somatic_snvs.vcf.gz'
    elif variant_caller == 'freebayes':
        rel_path = f'{path}/{variant_caller}/{sample_name}_vs_normal_sample'
        name = f'{sample_name}_vs_normal_sample.freebayes.vcf.gz'
    
    if cna_caller == 'ascat':
        path_cn = f'{path}/ascat/{sample_name}_vs_normal_sample'
        segment = f'{sample_name}_vs_normal_sample.segments.txt'
        purity = f'{sample_name}_vs_normal_sample.purityploidy.txt'
        cn_caller = 'ASCAT'
        
    elif cna_caller == 'sequenza':
        path_cn = f'{path}/sequenza/{sample_name}_vs_normal_sample'
        segment = f'{sample_name}_vs_normal_sample_segments.txt'
        purity = f'{sample_name}_vs_normal_sample_confints_CP.txt'
        cn_caller = 'sequenza'
    
    if  variant_caller == 'mutect2':
        tumourevo_file.write(f'\nSCOUT,{SPN},{SPN}_{sample_name},{SPN}_normal_sample,{rel_path}/{name},{rel_path}/{name}.tbi,{path_cn}/{segment},{path_cn}/{purity},{cn_caller},{cancer_type}')
    else:
        cram_tumour = f'{base_path}/preprocessing/recalibrated/{sample_name}/{sample_name}.recal.cram'
        tumourevo_file.write(f'\nSCOUT,{SPN},{SPN}_{sample_name},{SPN}_normal_sample,{rel_path}/{name},{rel_path}/{name}.tbi,{path_cn}/{segment},{path_cn}/{purity},{cn_caller},{cancer_type},{cram_tumour},{cram_tumour}.crai')



if (__name__ == '__main__'):
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     description=('Produces the cohorts of a SPN'))
    parser.add_argument('SPN', type=str, help='The SPN name (e.g., SPN01)')
    parser.add_argument('phylogenetic_forest', type=str,
                        help = ('A ProCESS phylogenetic forest'))
    parser.add_argument('output_dir', type=str,
                        help = ('The output directory'))
    parser.add_argument('-P', '--partition', type=str, required=True,
                        help="The cluster partition")
    parser.add_argument('-A', '--account', type=str, required=True,
                        help="The cluster account")
    parser.add_argument('-s', '--node_scratch_directory', type=str,
                        default='/local_scratch',
                        help="The nodes' scratch directory")
    parser.add_argument('-j', '--parallel_jobs', type=int, default=40,
                        help="The number of parallel jobs")
    parser.add_argument('-x', '--exclude', type=str, default="genoa011,genoa008",
                        help=("A list of nodes to exclude from the "
                              + "computation"))
    parser.add_argument('-F', '--force_completed_jobs', action='store_true',
                        help=("A Boolean flag to force rerun of "
                              + "already completed job."))
    parser.add_argument('-S', '--scratch_per_node', type=float, default=300,
                        help=("The scratch space available in each "
                              + "node (in GB)."))
    parser.add_argument('-M', '--mem_per_node', type=float, default=512,
                        help="The memory of each node in GB")
    parser.add_argument('-I', '--image_path', type=str, default="",
                        help="Path to singularity image")
    parser.add_argument('-C', '--config', type=str, default="",
                        help="Path to nextflow config file")
    parser.add_argument('-SD', '--sarek_output_dir', type=str, default="",
                       help="Path to sarek launching dir")
    parser.add_argument('-TD', '--tumourevo_output_dir', type=str, default="",
                       help="Path to tumourevo result path")
    
    cohorts = { 'normal': {
                    'max_coverage': 30,
                    'purities': list([1])
                    },
               'tumour': {
                    'max_coverage': 200,
                    'purities': [0.3,0.6,0.9]
                    }
                }

    num_of_lots_T = 40
    num_of_lots_N = 6
    
    cohort_coverages = list([50, 100, 150, 200])
    args = parser.parse_args()

    if args.account is None:
        process = subprocess.Popen(['whoami'],
                                stdout=subprocess.PIPE)
        account = process.communicate()
    else:
        account = args.account

    gender_filename = os.path.join(os.path.dirname(args.phylogenetic_forest),
                                   "subject_gender.txt")

    curr_dir = os.getcwd()

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    
    sarek_dir = os.path.join(args.output_dir, 'sarek')
    tumourevo_dir = os.path.join(args.output_dir, 'tumourevo')
    
    tumour_type_file=os.path.join(os.path.dirname(args.phylogenetic_forest),"tumour_type.txt")

    #if not os.path.exists(tumour_type_file):

    # with open('ProCESS_tumour_type.sh', 'w') as outstream:
    #     outstream.write(write_tumour_type_file)

    # cmd = ['sbatch', '--account={}'.format(account),
    #       '--partition={}'.format(args.partition),
    #       ('--export=DIR={},SPN={},BASEDIR={}').format(curr_dir,args.SPN,
    #                                             os.path.dirname(args.phylogenetic_forest)),
    #       './ProCESS_tumour_type.sh']

    # subprocess.run(cmd)

    config_file = args.config
    if not os.path.exists(sarek_dir):
        os.mkdir(sarek_dir)
        
    if not os.path.exists(tumourevo_dir):
        os.mkdir(tumourevo_dir)  
        
    for seq_type, cohorts_data in cohorts.items():
        if seq_type == 'normal':
            num_of_lots = num_of_lots_N
        else:
            num_of_lots = num_of_lots_T
            
        zeros = math.ceil(math.log10(num_of_lots))
        
        lot_coverage = cohorts_data['max_coverage']/num_of_lots
        lot_prefix = get_lot_prefix(seq_type)
        type_output_dir = f'{args.output_dir}/{seq_type}'


        for purity in cohorts_data['purities']:

            if seq_type == 'normal':
                job_id=seq_type
                sarek_file_launcher_orig = sarek_file_normal_launcher
                sarek_file_normal_launcher = sarek_file_normal_launcher.replace('{ACCOUNT}', str(account))
                sarek_file_normal_launcher = sarek_file_normal_launcher.replace('{JOB_NAME}', str(job_id))
                sarek_file_normal_launcher = sarek_file_normal_launcher.replace('{INPUT_DIR}', str(sarek_dir))
                sarek_file_normal_launcher = sarek_file_normal_launcher.replace('{CONFIG}', str(config_file))
                sarek_file_normal_launcher = sarek_file_normal_launcher.replace('{SAREK_OUT}', str(args.sarek_output_dir))

                with open(f'{sarek_dir}/sarek_mapping_vc_normal.sh', 'w') as outstream:
                    outstream.write(sarek_file_normal_launcher)
                sarek_file_normal_launcher = sarek_file_launcher_orig

                    
            else:
                tumour_fastq_dir = os.path.join(f'{args.output_dir}', f'tumour/purity_{purity}/FASTQ')
                sample_names = get_sample_names_from_FASTQ(tumour_fastq_dir)
                for cohort_cov in cohort_coverages:

                    #sarek mapping sh file
                    sarek_file_launcher_orig = sarek_file_launcher    
                    job_id=f'{cohort_cov}x_{purity}p'

                    sarek_file_launcher = sarek_file_launcher.replace('{ACCOUNT}', str(account))
                    sarek_file_launcher = sarek_file_launcher.replace('{JOB_NAME}', str(job_id))
                    sarek_file_launcher = sarek_file_launcher.replace('{INPUT_DIR}', str(sarek_dir))
                    sarek_file_launcher = sarek_file_launcher.replace('{CONFIG}', str(config_file))
                    sarek_file_launcher = sarek_file_launcher.replace('{SAREK_OUT}', str(args.sarek_output_dir)) 

                    with open(f'{sarek_dir}/sarek_mapping_{cohort_cov}x_{purity}p.sh', 'w') as outstream:
                        outstream.write(sarek_file_launcher)
                    sarek_file_launcher = sarek_file_launcher_orig
                    
                    #sarek VC sh file
                    sarek_variant_calling_launcher_orig = sarek_variant_calling_launcher
                    job_id=f'{cohort_cov}x_{purity}p'
                    
                    sarek_variant_calling_launcher = sarek_variant_calling_launcher.replace('{ACCOUNT}', str(account))
                    sarek_variant_calling_launcher = sarek_variant_calling_launcher.replace('{JOB_NAME}', str(job_id))
                    sarek_variant_calling_launcher = sarek_variant_calling_launcher.replace('{INPUT_DIR}', str(sarek_dir))
                    sarek_variant_calling_launcher = sarek_variant_calling_launcher.replace('{CONFIG}', str(config_file))
                    sarek_variant_calling_launcher = sarek_variant_calling_launcher.replace('{SAREK_OUT}', str(args.sarek_output_dir))                    
                    with open(f'{sarek_dir}/sarek_variant_calling_{cohort_cov}x_{purity}p.sh', 'w') as outstream:
                        outstream.write(sarek_variant_calling_launcher)
                    sarek_variant_calling_launcher = sarek_variant_calling_launcher_orig
                    
                    # sequenza sh file
                    sequenza_launcher_orig = sequenza_launcher
                    job_id=f'{cohort_cov}x_{purity}p'
                    process_path = os.path.join('/'.join(os.path.normpath(config_file).split(os.path.sep)[:-2]), 'sequenza')
                    
                    sequenza_launcher = sequenza_launcher.replace('{ACCOUNT}', str(account))
                    sequenza_launcher = sequenza_launcher.replace('{JOB_NAME}', str(job_id))
                    sequenza_launcher = sequenza_launcher.replace('{INPUT_DIR}', str(sarek_dir))
                    sequenza_launcher = sequenza_launcher.replace('{CONFIG}', str(config_file))
                    sequenza_launcher = sequenza_launcher.replace('{SAREK_OUT}', str(args.sarek_output_dir))
                    sequenza_launcher = sequenza_launcher.replace('{PROCESS_DIR}', str(process_path))

                    
                    with open(f'{sarek_dir}/sequenza_{cohort_cov}x_{purity}p.sh', 'w') as outstream:
                        outstream.write(sequenza_launcher)
                    sequenza_launcher = sequenza_launcher_orig
                    
                    # battenberg sh file
                    battenberg_launcher_orig = battenberg_launcher
                    job_id=f'{cohort_cov}x_{purity}p'
                    process_path = os.path.join('/'.join(os.path.normpath(config_file).split(os.path.sep)[:-2]), 'sequenza')
                    
                    battenberg_launcher = battenberg_launcher.replace('{ACCOUNT}', str(account))
                    battenberg_launcher = battenberg_launcher.replace('{JOB_NAME}', str(job_id))
                    battenberg_launcher = battenberg_launcher.replace('{INPUT_DIR}', str(sarek_dir))
                    battenberg_launcher = battenberg_launcher.replace('{CONFIG}', str(config_file))
                    battenberg_launcher = battenberg_launcher.replace('{SAREK_OUT}', str(args.sarek_output_dir))
                    battenberg_launcher = battenberg_launcher.replace('{PROCESS_DIR}', str(process_path))

                    
                    with open(f'{sarek_dir}/battenberg_{cohort_cov}x_{purity}p.sh', 'w') as outstream:
                        outstream.write(battenberg_launcher)
                    battenberg_launcher = battenberg_launcher_orig
                    
                    with open(tumour_type_file) as cancer_type_file:
                        cancer_type = cancer_type_file.read().strip()
                    #tumourevo sh file and csv file
                    variant_callers = ['freebayes', 'strelka', 'mutect2']
                    cn_caller = ['ascat', 'sequenza']
                    combinations = []
                    for vc in variant_callers:
                        for cnc in cn_caller:
                            combinations.append([vc, cnc])


                    for comb in combinations:
                        vc = comb[0]
                        cc = comb[1]
                        # with open(f'{tumourevo_dir}/tumourevo_{cohort_cov}x_{purity}p_{vc}_{cc}.csv', 'w') as tumourevo_file:
                        #     if comb[0] == 'mutect2':
                        #         tumourevo_file.write('dataset,patient,tumour_sample,normal_sample,vcf,tbi,cna_segments,cna_extra,cna_caller,cancer_type')
                        #     else:
                        #         tumourevo_file.write('dataset,patient,tumour_sample,normal_sample,vcf,tbi,cna_segments,cna_extra,cna_caller,cancer_type,tumour_alignment,tumour_alignment_index')
                        #     for sample_name in sample_names:
                        #         write_tumourevo_lines(tumourevo_file, args.SPN, sample_name, comb, cohort_cov, purity, args.sarek_output_dir,cancer_type)
                    
                        tumourevo_launcher_orig = tumourevo_launcher
                        job_id=f'{cohort_cov}x_{purity}p_{vc}_{cc}'
                        tumourevo_launcher = tumourevo_launcher.replace('{ACCOUNT}', str(account))
                        tumourevo_launcher = tumourevo_launcher.replace('{JOB_NAME}', str(job_id))
                        tumourevo_launcher = tumourevo_launcher.replace('{INPUT_DIR}', str(tumourevo_dir))
                        tumourevo_launcher = tumourevo_launcher.replace('{CONFIG}', str(config_file))
                        tumourevo_launcher = tumourevo_launcher.replace('{TUMOUREVO_OUT}', str(args.tumourevo_output_dir))


                        with open(f'{tumourevo_dir}/tumourevo_{cohort_cov}x_{purity}p_{vc}_{cc}.sh', 'w') as outstream:
                            outstream.write(tumourevo_launcher)
                        tumourevo_launcher = tumourevo_launcher_orig
