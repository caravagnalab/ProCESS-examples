#!/usr/bin/python3

import os
import sys
import math
import glob
import time
from pathlib import Path
import re


import subprocess
import argparse
import fnmatch

somatic_processing_shell_script="""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem=50GB

module load R/4.4.1 
echo "Processing chromosome ${CHROMOSOME}"

cd ${DIRECTORY}/Somatic/
Rscript ${DIRECTORY}/Somatic/preprocess.R --spn_id ${SPN} --purity ${PURITY} --coverage ${COVERAGE} --chr ${CHROMOSOME}
"""

somatic_prepare_report_shell_script="""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 50g
#SBATCH --time=4:00:00

module load R/4.4.1

cd ${DIRECTORY}/Somatic/
Rscript ${DIRECTORY}/Somatic/prepare_report.R --spn_id ${SPN} --purity ${PURITY} --coverage ${COVERAGE}
"""

somatic_prepare_multi_caller_report_shell_script="""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 50g
#SBATCH --time=4:00:00

module load R/4.4.1

cd ${DIRECTORY}/Somatic/
Rscript ${DIRECTORY}/Somatic/prepare_multi_caller_report.R --spn_id ${SPN} --purity ${PURITY} --coverage ${COVERAGE}
"""


cna_report_shell_script="""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 50g
#SBATCH --time=4:00:00

module load R/4.4.1

cd ${DIRECTORY}/CNA/
Rscript ${DIRECTORY}/CNA/Validate_CNA_calls.R --spn_id ${SPN} --purity ${PURITY} --coverage ${COVERAGE}
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
    parser.add_argument('-S', '--skip', type=str, default='',
                        help="Which step to skip, select among: cna, somatic, none")

    args = parser.parse_args()
    
    def get_sample_names(spn, base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"):
        search_path = Path(base_path) / spn / "process" / spn
        name_files = [f.name for f in search_path.glob("*.rff")]
        sample_names = [re.sub(r'.{4}$', '', name) for name in name_files]
        return sample_names


    def required_preprocessing(chr, search_path,expected_files):
        filename="chr"+str(chr)+".rds"
        founded_files = []
        for root, dirs, files in os.walk(search_path):
            if filename in files:
                founded_files.append(os.path.join(root, filename))
        if (len(founded_files)==expected_files):
            return False
        return True


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
    
    #somatic_callers = ["mutect2","strelka",""]

    if "somatic" not in args.skip.split(","):
        with open('run_processing.sh', 'w') as outstream:
            outstream.write(somatic_processing_shell_script)

        chr_job_ids = []
        outfile_preprocessing = outfile_preprocessing = scout_dir+"/"+args.SPN+"/validation/somatic/"+args.SPN+"/"+str(args.coverage)+"x_"+args.purity+"p"
        n_exp =len(get_sample_names(args.SPN,scout_dir))*4*2
        

        for chr in chromosomes:
            
            if (required_preprocessing(chr, outfile_preprocessing,n_exp)):
              
              print("Run preprocessing for chromosome "+str(chr))
              
              cmd = ['sbatch','--parsable',
                  '--account={}'.format(account),
                  '--partition={}'.format(args.partition),
                  '--job-name=parsing_{}_{}_{}_{}'.format(args.SPN, args.coverage,args.purity,chr),
                  ('--export=SPN={},COVERAGE={},PURITY={},CHROMOSOME={},DIRECTORY={}').format(args.SPN,
                                                      args.coverage, args.purity,chr,args.directory),
                  '--output={}/parsing_chr{}.log'.format(log_dir, chr),
                  './run_processing.sh']

              result = subprocess.run(cmd, stdout=subprocess.PIPE)
              job_id = result.stdout.decode().strip()
              chr_job_ids.append(job_id)
            else:
              print("No need to run preprocessing for chromosome "+str(chr))
        
        # Run final somatic report

        with open('run_single_caller_reports.sh', 'w') as outstream:
            outstream.write(somatic_prepare_report_shell_script) ## single caller
        dependency_str = ':'.join(chr_job_ids)
        print("Run report")
        cmd = ['sbatch', '--parsable',
            '--account={}'.format(account),
            '--partition={}'.format(args.partition),
            '--job-name=somatic_report_{}_{}_{}'.format(args.SPN, args.coverage,args.purity),
            '--dependency=afterok:{}'.format(dependency_str),
            ('--export=SPN={},COVERAGE={},PURITY={},DIRECTORY={}').format(args.SPN,
                                                args.coverage, args.purity,args.directory),
            '--output={}/somatic_report_{}_{}_{}.log'.format(log_dir, args.SPN,args.coverage,args.purity),
            './run_single_caller_reports.sh']
        result = subprocess.run(cmd, stdout=subprocess.PIPE)
        job_id = result.stdout.decode().strip()
        
        with open('run_multicaller_reports.sh', 'w') as outstream:
            outstream.write(somatic_prepare_multi_caller_report_shell_script)     
        cmd = ['sbatch', '--parsable',
            '--account={}'.format(account),
            '--partition={}'.format(args.partition),
            '--job-name=somatic_multicaller_report_{}_{}_{}'.format(args.SPN, args.coverage,args.purity),
            '--dependency=afterok:{}'.format(job_id),
            ('--export=SPN={},COVERAGE={},PURITY={},DIRECTORY={}').format(args.SPN,
                                                args.coverage, args.purity,args.directory),
            '--output={}/somatic_multicaller_report_{}_{}_{}.log'.format(log_dir, args.SPN,args.coverage,args.purity),
            './run_multicaller_reports.sh']
        result = subprocess.run(cmd, stdout=subprocess.PIPE)
        job_id = result.stdout.decode().strip()

    
    if "cna" not in args.skip.split(","):

        ## Run final CNA report

        base_scout_dir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
        cna_dir = os.path.join(base_scout_dir, args.SPN, "process", "cna_data")
        cna_files = [f for f in os.listdir(cna_dir) if fnmatch.fnmatch(f, '*_cna.rds')]

        with open('validate_cna.sh', 'w') as outstream:
            outstream.write(cna_report_shell_script)
      
        cmd = ['sbatch', '--account={}'.format(account),
            '--partition={}'.format(args.partition),
            '--job-name=cna_report_{}_{}_{}'.format(args.SPN, args.coverage,args.purity),
            ('--export=SPN={},COVERAGE={},PURITY={},DIRECTORY={}').format(args.SPN,
                                                args.coverage, args.purity,args.directory),
            '--output={}/cna_report_{}_{}_{}.log'.format(log_dir, args.SPN,args.coverage,args.purity),
            './validate_cna.sh']
        subprocess.run(cmd)
