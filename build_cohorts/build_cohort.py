#!/usr/bin/python3

import os
import sys
import math
import glob
import time
import subprocess
import argparse

gender_shell_script="""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=20GB

module load R/4.3.3

echo "Rscript ProCESS_subject_gender.R ${PHYLO_FOREST}"


Rscript ProCESS_subject_gender.R ${PHYLO_FOREST}
"""

gender_R_script="""rm(list = ls())
library(ProCESS)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop(paste("Syntax error: ProCESS_subject_gender.R",
	         "<phylo_forest>"),
       call. = FALSE)
}

forest <- load_phylogenetic_forest(args[1])

dir <- dirname(args[1])

gender <- forest$get_germline_subject()$gender
if (gender == "male") {
    gender <- "XY"
} else if (gender == "female") { 
    gender <- "XY"
} else {
    stop(paste0("Unsupported germline subject gender ",
                "\\"",gender,"\\"."),
         call. = FALSE)
}

fileConn<-file(file.path(dir, "subject_gender.txt"))
writeLines(c(gender), fileConn)
close(fileConn)
"""


shell_script="""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --time=8:00:00
#SBATCH --mem={MEMORY}GB

module load R/4.3.3
module load samtools

echo "Rscript ProCESS_seq.R ${PHYLO_FOREST} ${SPN} ${LOT} ${NODE_SCRATCH} ${DEST} ${COVERAGE} ${TYPE} 4 ${SEED} ${PURITY}"


Rscript ProCESS_seq.R ${PHYLO_FOREST} ${SPN} ${LOT} ${NODE_SCRATCH} ${DEST} ${COVERAGE} ${TYPE} 4 ${SEED} ${PURITY}

rm -rf ${NODE_SCRATCH}/${SPN}_${LOT}
"""

R_script="""rm(list = ls())
library(ProCESS)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 10) {
  stop(paste("Syntax error: ProCESS_seq.R",
	         "<phylo_forest> <SPN> <lot_name>",
	         "<node_local_dir> <output_dir>",
	         "<coverage> <type> <num_of_cores>",
	         "<seed> <purity>"),
       call. = FALSE)
}

phylo_forest_filename <- args[1]
spn_name <- args[2]
lot_name <- args[3]
node_local_dir <- args[4]
output_dir <- args[5]
coverage <- as.double(args[6])
type <- args[7]
num_of_cores <- strtoi(args[8])
seed <- strtoi(args[9])
purity <- as.double(args[10])

if (type == "tumour") {
    seq_tumour <- TRUE
} else if (type == "normal") {
    seq_tumour <- FALSE
    with_preneoplastic <- FALSE
} else if (type == "normal_with_preneoplastic") {
    seq_tumour <- FALSE
    with_preneoplastic <- TRUE
} else {
    stop(paste("The paramter <type> must be among",
               "\\"tumour\\", \\"normal\\", and
               \\"normal_with_preneoplastic\\"."),
         call. = FALSE)
}

merge_sams <- function(output_local_dir, BAM_file,
                       sam_filename_prefix, chromosomes,
		       num_of_cores) {
    
    SAM_files <- ""
    for (i in 1:length(chromosomes)) {
        chr_SAM_file <- file.path(output_local_dir,
                                  paste0(sam_filename_prefix,
                                         chromosomes[i], ".sam"))

        SAM_files <- paste(SAM_files, chr_SAM_file)
    }
    
    cmd <- paste("samtools merge -fc -@", num_of_cores,
		 "-o", BAM_file, SAM_files)

    invisible(system(cmd, intern = TRUE))
}

delete_sams <- function(output_local_dir, sam_filename_prefix, chromosomes) {
    for (i in 1:length(chromosomes)) {
        chr_SAM_file <- file.path(output_local_dir,
                                  paste0(sam_filename_prefix,
                                         chromosomes[i], ".sam"))

        unlink(chr_SAM_file)
    }
}

if (!file.exists(node_local_dir)) {
  dir.create(node_local_dir)
}

output_local_dir <- file.path(node_local_dir,
                              paste0(spn_name, "_",
                                     lot_name))
if (file.exists(output_local_dir)) {
  unlink(output_local_dir, recursive=TRUE)
}
dir.create(output_local_dir)

if (!file.exists(output_dir)) {
  dir.create(output_dir)
}

bam_dir <- file.path(output_dir, "BAM")
if (!file.exists(bam_dir)) {
  dir.create(bam_dir)
}

fastq_dir <- file.path(output_dir, "FASTQ")
if (!file.exists(fastq_dir)) {
  dir.create(fastq_dir)
}

data_dir <- file.path(output_dir, "data")
if (!file.exists(data_dir)) {
  dir.create(data_dir)
}

set.seed(seed)

filename_prefix <- lot_name

sam_filename_prefix <- paste0(filename_prefix, "_chr_")

BAM_filename <- paste0(filename_prefix, ".bam")

BAM_file <- file.path(bam_dir, BAM_filename)
BAM_local_file <- file.path(output_local_dir, BAM_filename)

BAM_done_filename <- file.path(output_dir, paste0(lot_name, "_BAM.done"))

step <- 1

if (!file.exists(BAM_done_filename) || !file.exists(BAM_file)) {
    unlink(BAM_done_filename)

    cat("1. Reading phylogenetic forest...\\n")
    phylo_forest <- load_phylogenetic_forest(phylo_forest_filename)
    
    cat("done\\n2. Copying reference genome...")
    
    ref_path <- file.path(output_local_dir, "reference.fasta")
    
    invisible(file.copy(phylo_forest$get_reference_path(), ref_path))
    
    cat("done\\n3. Simulating reads...\\n")
    
    # Simulate sequencing ####
    #no_error_seq <- ErrorlessIlluminaSequencer()
    basic_seq <- BasicIlluminaSequencer(1e-3) ## only for testing purpose
    chromosomes <- phylo_forest$get_absolute_chromosome_positions()$chr
    if (seq_tumour) {
      seq_results <- parallel::mclapply(chromosomes, function(c) {
        simulate_seq(phylo_forest, reference_genome = ref_path,
    	             chromosomes = c,
                     coverage = coverage,
                     purity = purity, 
                     write_SAM = TRUE, read_size = 150,
                     sequencer = basic_seq,
                     insert_size_mean = 350,
                     insert_size_stddev = 10,
                     output_dir = output_local_dir,
                     update_SAM = TRUE,
                     filename_prefix = sam_filename_prefix,
                     template_name_prefix = paste0(lot_name,'r'),
                     with_normal_sample = FALSE)
      }, mc.cores = num_of_cores)
    } else {
      seq_results <- parallel::mclapply(chromosomes, function(c) {
        simulate_normal_seq(phylo_forest, reference_genome = ref_path,
                            chromosomes = c,
                            coverage = coverage,
                            write_SAM = TRUE, read_size = 150,
                            sequencer = basic_seq,
                            insert_size_mean = 350,
                            insert_size_stddev = 10,
                            filename_prefix = sam_filename_prefix,
                     	    template_name_prefix = paste0(lot_name,'r'),
                            output_dir = output_local_dir,
                            with_preneoplastic = with_preneoplastic,
                            update_SAM = TRUE)
      }, mc.cores = num_of_cores)
    }
    seq_results_final<- do.call("bind_rows", seq_results)
    saveRDS(seq_results_final,
            file.path(data_dir,
                      paste0("seq_results_", spn_name,
    			  "_", lot_name, ".rds")))
    
    cat("done\\n4. Building overall BAM file...")
    merge_sams(output_local_dir, BAM_local_file,
               sam_filename_prefix, chromosomes,
    		   num_of_cores)
    
    cat("done\\n5. Deleting SAM files...")
    delete_sams(output_local_dir, sam_filename_prefix, chromosomes)
    
    cat("done\\n6. Moving the BAM file to output directory...")

    cmd <- paste0("cp ", BAM_local_file, " ", bam_dir, "/")

    invisible(system(cmd, intern = TRUE))

    invisible(file.create(BAM_done_filename))
    cat("done\\n")

    remove_local_bam <- TRUE
    step <- 7
} else {
    BAM_local_file <- BAM_file

    cat("Found the lot BAM file\\n")
    
    remove_local_bam <- FALSE

    step <- 1
}

cat(paste0(step, ". Splitting BAM file by sample..."))
step <- step + 1

split_bam_by_samples <- function(output_local_dir, BAM_file, remove_local_bam) {
    cmd <- paste0("samtools split -f \\"",
                  file.path(output_local_dir,"%*_%!.bam"),
                  "\\" ", BAM_file, " -@ ", num_of_cores)
    invisible(system(cmd, intern = TRUE))

    if (remove_local_bam) {
        file.remove(BAM_file)
    }
}
invisible(split_bam_by_samples(output_local_dir, BAM_local_file, remove_local_bam))

cat(paste0("done\\n", step,
           ". Generating the FASTQs and deleting the BAMs..."))
step <- step + 1

BAM_files <- list.files(output_local_dir, pattern = "\\\\.bam$")

generate_fastq <- function(orig_file, fastq_dir) {
  base_orig_file <- tools::file_path_sans_ext(basename(orig_file))

  file_prefix <- file.path(fastq_dir, base_orig_file)
  R1 <- paste0(file_prefix, ".R1.fastq.gz")
  R2 <- paste0(file_prefix, ".R2.fastq.gz")
  unpaired <- paste0(file_prefix, ".unpaired.fastq.gz")
  singleton <- paste0(file_prefix, ".singleton.fastq.gz")

  cmd <- paste("samtools fastq -@ 20 -c 9 -N -1", R1, "-2", R2, "-0", unpaired, 
               "-s", singleton, orig_file)
  invisible(system(cmd, intern = TRUE))
}

result <- parallel::mclapply(BAM_files, function(c) {
    curr_BAM_file <- file.path(output_local_dir, c)
    if (BAM_file != curr_BAM_file) {
        generate_fastq(curr_BAM_file, output_local_dir)

        unlink(curr_BAM_file)
    }
}, mc.cores = num_of_cores)

cat(paste0("done\\n", step,
           ". Moving the FASTQ files to output directory..."))
step <- step + 1

cmd <- paste0("mv ", file.path(output_local_dir, "*.fastq.gz"),
              " ", fastq_dir, "/")
invisible(system(cmd, intern = TRUE))

cat(paste0("done\\n", step, ". Removing local files..."))
step <- step + 1

unlink(output_local_dir, recursive = TRUE)

done_filename <- file.path(output_dir, paste0(lot_name, "_final.done"))
invisible(file.create(done_filename))

cat("done\\n")
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


def get_completed_jobs(done_file_dir, lot_prefix):
    common_prefix = os.path.normpath(f"{done_file_dir}/{lot_prefix}")
    common_suffix = '_final.done'
    done_files = glob.glob(f"{common_prefix}*{common_suffix}")
    prefix_len = len(common_prefix)
    suffix_len = len(common_suffix)

    done_ids = list()
    for done_file in done_files:
        done_ids.append(int(done_file[prefix_len:-suffix_len]))
    return done_ids


def remove_old_done_files(output_dir, lot_prefix):
    done_files = glob.glob(f"{output_dir}/{lot_prefix}*.done")
    for done_file in done_files:
        os.unlink(done_file)

def get_sample_names_from_FASTQ(fastq_dir):
    suffix = '.R1.fastq.gz'
    fastq_files = glob.glob(f'{fastq_dir}/t*_*{suffix}')

    sample_names = set()
    for fastq_file in fastq_files:
        fastq_basename = os.path.basename(fastq_file)
        prefix_up_to = fastq_basename.find('_')
        
        sample_names.add(fastq_basename[prefix_up_to+1:-len(suffix)])
    
    return sorted(list(sample_names))

def write_sarek_sample_lines(sarek_file, SPN, seq_type, sample_name, num_of_lots, fastq_dir,
                             lot_padding_zeros, line_padding_zeros=2):
    if (seq_type == 'tumour'):
        status = 1
    elif (seq_type == 'normal'):
        status = 0
    else:
        raise RuntimeError(f'Unsupported sequence type "{seq_type}"')

    fastq_suffix = '.fastq.gz'

    line = 1
    for lot in range(num_of_lots):
        lot_name = f'{get_lot_prefix(seq_type)}{str(lot).zfill(lot_padding_zeros)}'
        fastq_base_name = f'{lot_name}_{sample_name}'
        line_name = f'L{str(line).zfill(line_padding_zeros)}'
        line += 1
        R1_filename = os.path.abspath(os.path.join(fastq_dir,
                                                    fastq_base_name+'.R1' + fastq_suffix))
        R2_filename = os.path.abspath(os.path.join(fastq_dir,
                                                    fastq_base_name+'.R2' + fastq_suffix))
        sarek_file.write(f'\n{SPN},{subject_gender},{status},{sample_name},'
                        + f'{line_name},{R1_filename},{R2_filename}')

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
    parser.add_argument('-x', '--exclude', type=str, default="",
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

    cohorts = { 'tumour': {
                    'max_coverage': 200, 
                    'purities': list([0.3, 0.6, 0.9])
                    },
                'normal': {
                    'max_coverage': 50, 
                    'purities': list([1])
                    }
                }

    num_of_lots = 40
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

    if not os.path.exists(gender_filename):
        with open('ProCESS_subject_gender.R', 'w') as outstream:
            outstream.write(gender_R_script)

        with open('ProCESS_subject_gender.sh', 'w') as outstream:
            outstream.write(gender_shell_script)

        cmd = ['sbatch', '--account={}'.format(account),
            '--partition={}'.format(args.partition),
            f'--export=PHYLO_FOREST={args.phylogenetic_forest}',
            './ProCESS_subject_gender.sh']

        subprocess.run(cmd)

    with open('ProCESS_seq.R', 'w') as outstream:
        outstream.write(R_script)

    space_per_lot = 3 * cohorts['tumour']['max_coverage'] * 5 / num_of_lots
    memory_per_lot = math.ceil(args.mem_per_node*space_per_lot/args.scratch_per_node)
    memory_per_lot = max(memory_per_lot, math.ceil(args.mem_per_node/5))

    shell_script = shell_script.replace('{MEMORY}', str(memory_per_lot))

    with open('ProCESS_seq.sh', 'w') as outstream:
        outstream.write(shell_script)

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    zeros = math.ceil(math.log10(num_of_lots))

    for seq_type, cohorts_data in cohorts.items():
        lot_coverage = cohorts_data['max_coverage']/num_of_lots

        lot_prefix = get_lot_prefix(seq_type)

        type_output_dir = f'{args.output_dir}/{seq_type}'

        if not os.path.exists(type_output_dir):
            os.mkdir(type_output_dir)

        for purity in cohorts_data['purities']:

            output_dir = f'{type_output_dir}/purity_{purity}'

            if not os.path.exists(output_dir):
                os.mkdir(output_dir)

            log_dir = '{}/log/'.format(output_dir)

            if not os.path.exists(log_dir):
                os.mkdir(log_dir)

            if args.force_completed_jobs:
                remove_old_done_files(output_dir, lot_prefix)
            
            lot_ids = set(range(num_of_lots))

            completed_ids = get_completed_jobs(output_dir, lot_prefix)
            submitted = list(completed_ids)
            lot_ids = list(lot_ids.difference(set(completed_ids)))

            while len(lot_ids) != 0:
                completed_ids = get_completed_jobs(output_dir, lot_prefix)
                
                to_be_submitted = (args.parallel_jobs
                                + len(completed_ids)
                                - len(submitted))

                for i in lot_ids[:to_be_submitted]:
                    lot_name = '{}{}'.format(lot_prefix, str(i).zfill(zeros))

                    sys.stdout.write('Submitting lot {}...'.format(lot_name))
                    sys.stdout.flush()

                    cmd = ['sbatch', '--account={}'.format(account),
                        '--partition={}'.format(args.partition),
                        '--job-name={}_{}'.format(args.SPN, lot_name),
                        ('--export=PHYLO_FOREST={},SPN={},LOT={},DEST={},'
                            + 'COVERAGE={},TYPE={},NODE_SCRATCH={},'
                            + 'SEED={},PURITY={}').format(args.phylogenetic_forest,
                                                args.SPN, lot_name,
                                                output_dir, lot_coverage,
                                                seq_type, args.node_scratch_directory,
                                                i, purity),
                        '--output={}/lot_{}.log'.format(log_dir, lot_name),
                        './ProCESS_seq.sh']
                    if args.exclude != "":
                        cmd.insert(-1,"--exclude={}".format(args.exclude))

                    subprocess.run(cmd)
                    sys.stdout.write('done\n')
                    sys.stdout.flush()

                if to_be_submitted>0:
                    submitted.extend(lot_ids[:to_be_submitted])
                    lot_ids = lot_ids[to_be_submitted:]
                time.sleep(60)  # wait 1 minute
            completed_ids = get_completed_jobs(output_dir, lot_prefix)
            while (len(completed_ids) != len(submitted)):
                time.sleep(60)
                completed_ids = get_completed_jobs(output_dir, lot_prefix)

    with open(gender_filename, "r") as gender_file:
        subject_gender = gender_file.read().strip('\n')
    
    sarek_dir = os.path.join(args.output_dir, 'sarek')
    if not os.path.exists(sarek_dir):
        os.mkdir(sarek_dir)

    normal_fastq_dir = os.path.join(f'{args.output_dir}', 'normal/purity_1/FASTQ')
    for purity in cohorts['tumour']['purities']:
        tumour_fastq_dir = os.path.join(f'{args.output_dir}', f'tumour/purity_{purity}/FASTQ')
        sample_names = get_sample_names_from_FASTQ(tumour_fastq_dir)

        lines = math.ceil(math.log10(num_of_lots*(len(sample_names)+1)))

        for cohort_cov in cohort_coverages:
            num_of_tumour_lots = math.ceil((cohort_cov*num_of_lots)/cohorts['tumour']['max_coverage'])
            with open(f'{sarek_dir}/sarek_{cohort_cov}x_{purity}p.csv', 'w') as sarek_file:
                sarek_file.write('patient,sex,status,sample,lane,fastq_1,fastq_2')
                for sample_name in sample_names:
                    write_sarek_sample_lines(sarek_file, args.SPN, 'tumour', sample_name,
                                             num_of_tumour_lots, tumour_fastq_dir, zeros, lines)
                write_sarek_sample_lines(sarek_file, args.SPN, 'normal', 'normal_sample',
                                         num_of_lots, normal_fastq_dir, zeros, lines)
