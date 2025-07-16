#!/usr/bin/python3

import os
import sys
import math
import glob
import time

import subprocess
import argparse
import fnmatch

germline_processing_shell_script="""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 50g
#SBATCH --time=4:00:00

outdir="/orfeo/scratch/cdslab/shared/SCOUT/${SPN}/validation/germline/vcf"
mkdir -p $outdir

tool="haplotypecaller"
vcf="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${SPN}/sarek/normal/variant_calling/${tool}/normal_sample/normal_sample.${tool}.filtered.vcf.gz"
bcftools view ${vcf} --regions chr${CHROMOSOME} -o ${outdir}/chr${CHROMOSOME}_normal_sample.${tool}.vcf.gz -Oz
        
tool="freebayes"
vcf="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${SPN}/sarek/normal/variant_calling/${tool}/normal_sample/normal_sample.${tool}.vcf.gz"
bcftools view ${vcf} --regions chr${CHROMOSOME} -o ${outdir}/chr${CHROMOSOME}_normal_sample.${tool}.vcf.gz -Oz

tool="strelka"
vcf="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${SPN}/sarek/${COVERAGE}x_${PURITY}p/variant_calling/${tool}/normal_sample/normal_sample.${tool}.variants.vcf.gz"
bcftools view ${vcf} --regions chr${CHROMOSOME} -o ${outdir}/chr${CHROMOSOME}_normal_sample.${tool}.vcf.gz -Oz
"""

germline_report_shell_script="""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem 300g
#SBATCH --time=10:00:00

module load R/4.4.1

cd ${DIRECTORY}/Germline/
Rscript ${DIRECTORY}/Germline/compare_caller.R -s ${SPN}
"""

if (__name__ == '__main__'):
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     description=('Produces the cohorts of a SPN'))
    parser.add_argument('SPN', type=str, help='The SPN name (e.g., SPN01)')
    parser.add_argument('coverage', type=int, help='The coverage')
    parser.add_argument('purity', type=str, help='The purity')
    parser.add_argument('directory', type=str,
                        help="The main directory of ProCESS-examples")
    parser.add_argument('-P', '--partition', type=str, required=True,
                        help="The cluster partition")
    parser.add_argument('-A', '--account', type=str, required=True,
                        help="The cluster account")
                        
    args = parser.parse_args()

    if args.account is None:
        process = subprocess.Popen(['whoami'],
                                stdout=subprocess.PIPE)
        account = process.communicate()
    else:
        account = args.account
   
    scout_dir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
    with open(os.path.join(scout_dir, args.SPN, "process", "subject_gender.txt")) as gender_file:
      gender = gender_file.read().strip()
      
    if gender == "XX":
      chromosomes = list(range(1, 23)) + ['X']
    elif gender == "XY":
      chromosomes = list(range(1, 23)) + ['X', 'Y']
    else:
      raise ValueError(f"Unexpected gender value: {gender}")
    print(chromosomes)
    
    #curr_dir = os.getcwd()
    base_dir = args.directory
    log_dir = '{}/out/'.format(base_dir)
    
    os.makedirs(f'/orfeo/scratch/cdslab/shared/SCOUT/{args.SPN}/validation/germline/vcf', exist_ok = True)
    os.makedirs(f'/orfeo/scratch/cdslab/shared/SCOUT/{args.SPN}/validation/germline/report', exist_ok = True)
    if len(os.listdir(f'/orfeo/scratch/cdslab/shared/SCOUT/{args.SPN}/validation/germline/vcf')) != 72 and len(os.listdir(f'/orfeo/scratch/cdslab/shared/SCOUT/{args.SPN}/validation/germline/report')) != 6:
    
        ## Run processing for germline
        with open('run_processing_germline.sh', 'w') as outstream:
            outstream.write(germline_processing_shell_script)
        chr_job_ids_germline = []
        
        for chr in chromosomes:
            cmd_germline = ['sbatch', '--parsable','--account={}'.format(account),
                '--partition={}'.format(args.partition),
                '--job-name=split_chr_{}_{}_{}_{}'.format(args.SPN, args.coverage,args.purity,chr),
                ('--export=SPN={},COVERAGE={},PURITY={},CHROMOSOME={},DIRECTORY={}').format(args.SPN,
                                                    args.coverage, args.purity,chr,args.directory),
                '--output={}/split_chr{}.log'.format(log_dir, chr),
                './run_processing_germline.sh'] 
            result_germline = subprocess.run(cmd_germline, stdout=subprocess.PIPE)
            job_id_germline = result_germline.stdout.decode().strip()
            chr_job_ids_germline.append(job_id_germline)
            
        ## Run germline report
        with open('run_reports_germline.sh', 'w') as outstream:
            outstream.write(germline_report_shell_script)
        dependency_str_germline = ':'.join(chr_job_ids_germline)
        cmd = ['sbatch', '--account={}'.format(account),
            '--partition={}'.format(args.partition),
            '--job-name=germline_report_{}'.format(args.SPN),
            '--dependency=afterok:{}'.format(dependency_str_germline),
            ('--export=SPN={},DIRECTORY={}').format(args.SPN,args.directory),
            '--output={}/germline_report_{}.log'.format(log_dir, args.SPN),
            './run_reports_germline.sh']
        subprocess.run(cmd)
    
    else:
        print('Germline validation already exist!')
