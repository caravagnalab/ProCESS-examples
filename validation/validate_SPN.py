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

somatic_prepare_final_somatic_report="""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem=50GB

module load R/4.4.1 

cd ${DIRECTORY}/Somatic/

echo "Rscript ${DIRECTORY}/Somatic/prepare_final_report.R --spn_id ${SPN} --coverages ${COVERAGES} --purities ${PURITIES}"
Rscript ${DIRECTORY}/Somatic/prepare_final_report.R --spn_id ${SPN} --coverages ${COVERAGES} --purities ${PURITIES}
"""


#prepare_final_cna_report="""#!/bin/bash
##SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=1
##SBATCH --time=4:00:00
##SBATCH --mem=50GB
#
#module load R/4.4.1
#
#cd ${DIRECTORY}/Somatic/
#
#echo "Rscript ${DIRECTORY}/CNA/prepare_final_report.R --spn_id ${SPN} --coverages ${COVERAGES} --purities ${PURITIES}"
#Rscript ${DIRECTORY}/CNA/prepare_final_report.R --spn_id ${SPN} --coverages ${COVERAGES} --purities ${PURITIES}
#"""
#

if (__name__ == '__main__'):
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     description=('Produces the cohorts of a SPN'))
    parser.add_argument('SPN', type=str, help='The SPN name (e.g., SPN01)')
    parser.add_argument('coverages', type=str, help='The coverage')
    parser.add_argument('purities', type=str, help='The purity')
    parser.add_argument('directory', type=str)
    parser.add_argument('-P', '--partition', type=str, required=True,
                        help="The cluster partition")
    parser.add_argument('-A', '--account', type=str, required=True,
                        help="The cluster account")
    parser.add_argument('-S', '--skip', type=str, default='',
                        help="Which step to skip, select among: cna, somatic, none")

    args = parser.parse_args()

    if args.account is None:
        process = subprocess.Popen(['whoami'],
                                stdout=subprocess.PIPE)
        account = process.communicate()
    else:
        account = args.account
   
    scout_dir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT"

    base_dir = args.directory
    log_dir = '{}/out/'.format(base_dir)

    with open('run_final_somatic_report.sh', 'w') as outstream:
        outstream.write(somatic_prepare_final_somatic_report) ## single caller
    print("Run report")
    cmd = ['sbatch',
        '--account={}'.format(account),
        '--partition={}'.format(args.partition),
        '--job-name=somatic_final_{}'.format(args.SPN),
        ('--export=SPN={},COVERAGES="{}",PURITIES="{}",DIRECTORY={}').format(args.SPN,
                                                args.coverages,args.purities,
                                                args.directory),
        '--output={}/somatic_final_report_{}.log'.format(log_dir, args.SPN),
        './run_final_somatic_report.sh']
    subprocess.run(cmd)
