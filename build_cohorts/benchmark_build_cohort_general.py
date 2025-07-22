#!/usr/bin/python3

import os
import sys
import math
import glob
import time
import subprocess
import argparse

write_tumour_type_file="""#!/bin/bash
#SBATCH --partition={PARTITION}
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

## This part is currently run sequentially

merging_shell_script="""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --time=6:00:00
#SBATCH --mem=80GB

module load singularity
lots_list=${LOTS_LIST}

lots_list=($(echo $lots_list))

for i in ${lots_list[@]}
do
  echo "singularity exec --bind ${SINGULARITY_BIND} --no-home ${IMAGE} Rscript ${DIR}/ProCESS_merge_rds.R ${i} ${SPN} ${INPUT_DIR} ${PURITY} ${TYPE} ${MAX_COVERAGE} ${TOT_LOTS}"
  singularity exec --bind ${SINGULARITY_BIND} --no-home ${IMAGE} Rscript ${DIR}/ProCESS_merge_rds.R ${i} ${SPN} ${INPUT_DIR} ${PURITY} ${TYPE} ${MAX_COVERAGE} ${TOT_LOTS}
done
"""

merging_R_script=merging_R_script="""rm(list = ls())
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 7) {
  stop(paste("Syntax error: ProCESS_merge_rds.R",
             "<num_of_lots> <SPN> <input_dir> <purity> <type> <max_coverage> <tot_num_of_lots>"),
       call. = FALSE)
}

lot_end <- as.double(args[1])
spn <- args[2]
input_dir <- args[3]
purity <- args[4]
type <- args[5]
p_info <- ps::ps_handle()
start_time <- Sys.time()
initial_cpu <- ps::ps_cpu_times(p_info)
initial_mem <- ps::ps_memory_info(p_info)["rss"] / 1024^3
if (type=="tumour"){
  muts_dir <- paste0(input_dir,"/tumour/purity_",purity,"/data/mutations/")
  max_coverage <- as.double(args[6]) #200 ## this is hard-coded now
  num_of_lots <- as.double(args[7]) #40 ## this is hard-coded now
  coverage<-(max_coverage*lot_end)/num_of_lots
  data <- list()
 
  if (lot_end<=10){
    rds_files <- list.files(path = muts_dir, pattern = paste0("seq_results_muts_",spn,"_"),full.names = T)[1:lot_end]
    data <- lapply(rds_files,function(x){
      readRDS(x)  %>%
        dplyr::select(-ends_with(".VAF"))
    })
  } else{
    lot_start<-lot_end-10+1
    previous_lot <- lot_end-10
    previous_coverage <- (max_coverage*previous_lot)/num_of_lots
    print(previous_coverage)
    previous_lots <- list.files(path = muts_dir, pattern = paste0("seq_results_muts_merged_coverage_",previous_coverage),full.names = T)
    rds_files <- list.files(path = muts_dir, pattern = paste0("seq_results_muts_",spn,"_"),full.names = T)[lot_start:lot_end]
    
    rds_files_all <- c(rds_files,previous_lots)
    message(paste0("Merging the following rds: ",rds_files_all))
    data <- lapply(rds_files_all,function(x){
      readRDS(x)  %>%
        dplyr::select(-ends_with(".VAF"))
    })
  }
  
  ids <- grep(pattern = "coverage",x = colnames(data[[1]]),value = T) %>% strsplit("\\\.")
  print("Combining dataframes ...")
  combined_df <- bind_rows(data)
  sample_names <- sapply(ids, function(x) {
    parts <- unlist(strsplit(x, "\\\."))
    paste(parts[1], parts[2], sep = ".")
  })
  
  print("Summing up NV and DP...")
  
  columns  <- colnames(combined_df)[c(7:ncol(combined_df))]
  result <- combined_df %>%
    group_by(chr, chr_pos,ref,alt,classes,causes) %>% 
    summarize(across(all_of(columns), sum, .names = "{.col}"))
  
  print("Recalculate VAF...")
  for (s in sample_names){
    col_name_DP <- paste0(s,".coverage")
    col_name_NV <- paste0(s,".occurrences")
    col_name_VAF <- paste0(s,".VAF")
    result <- result %>% 
      mutate(!!col_name_VAF := .data[[col_name_NV]] / .data[[col_name_DP]])
    print(s)
  }
  print("Saving merged rds...")
  saveRDS(result, file = paste0(muts_dir,"seq_results_muts_merged_coverage_",coverage,"x", ".rds"))
  print("Done merging!")
} else if (type=="normal"){
  max_coverage <- as.double(args[6])
  num_of_lots <- as.double(args[7]) #40 ## this is hard-coded now
  coverage<-(max_coverage*lot_end)/num_of_lots
  muts_dir <- paste0(input_dir,"/normal/purity_1/data/mutations/")
  rds_files <- list.files(path = muts_dir, pattern = paste0("seq_results_muts_",spn,"_"),full.names = T)
  data <- lapply(rds_files,function(x){
    readRDS(x)  %>%
      dplyr::select(-ends_with(".VAF"))
  })
  print("Combining dataframes ...")
  combined_df <- bind_rows(data)

  print("Summing up NV and DP...")
  columns  <- colnames(combined_df)[c(7:ncol(combined_df))]

  result <- combined_df %>%
    group_by(chr, chr_pos,ref,alt,classes,causes) %>%
    summarize(across(all_of(columns), sum, .names = "{.col}"))
  
  print("Recalculate VAF...")
  s <- "normal_sample"
  col_name_DP <- paste0(s,".coverage")
  col_name_NV <- paste0(s,".occurrences")
  col_name_VAF <- paste0(s,".VAF")
  result <- result %>% 
    mutate(!!col_name_VAF := .data[[col_name_NV]] / .data[[col_name_DP]])
  saveRDS(result, file = paste0(muts_dir,"seq_results_muts_merged_coverage_",max_coverage,"x", ".rds"))
}

end_time <- Sys.time()
final_cpu <- ps::ps_cpu_times(p_info)
final_mem <- ps::ps_memory_info(p_info)["rss"] / 1024^3
elapsed_time <- end_time - start_time
elapsed_time <- as.numeric(elapsed_time, units = "mins")
cpu_used <- (final_cpu["user"] + final_cpu["system"]) - (initial_cpu["user"]+ initial_cpu["system"])
mem_used <- final_mem - initial_mem
resource_usage <- data.frame(
  elapsed_time_mins =  elapsed_time,
  cpu_time_secs = cpu_used,
  memory_used_MB = mem_used,
  coverage = coverage,
  purity = purity,
  type = type
)
rownames(resource_usage) <- NULL

if (type=="tumour"){
  data_dir_resources <- paste0(input_dir,"/tumour/purity_",purity,"/data/resources/")
  saveRDS(resource_usage, file = paste0(data_dir_resources,"seq_results_merging_coverage_",coverage,"x", ".rds"))
} else if (type=="normal"){
  data_dir_resources <- paste0(input_dir,"/normal/purity_1/data/resources/")
  saveRDS(resource_usage, file = paste0(data_dir_resources,"seq_results_merging_coverage_",max_coverage,"x", ".rds"))
}
"""

gender_shell_script="""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=40GB

module load singularity
singularity exec --bind ${SINGULARITY_BIND} --no-home ${IMAGE} Rscript ${DIR}/ProCESS_subject_gender.R ${PHYLO_FOREST}
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
    gender <- "XX"
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
#SBATCH --time=24:00:00
#SBATCH --mem={MEMORY}GB

module load singularity
singularity exec --bind ${SINGULARITY_BIND} --no-home ${IMAGE} Rscript ${DIR}/ProCESS_seq.R ${PHYLO_FOREST} ${SPN} ${LOT} ${NODE_SCRATCH} ${DEST} ${COVERAGE} ${TYPE} 4 ${SEED} ${PURITY}

rm -rf ${NODE_SCRATCH}/${SPN}_${LOT}
"""

R_script="""rm(list = ls())
library(ProCESS)
library(dplyr)
library(bench)

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
                       sam_filename_prefix, chromosomes, num_of_cores, resources_dir) {
    
    SAM_files <- ""
    for (i in 1:length(chromosomes)) {
        chr_SAM_file <- file.path(output_local_dir,
                                  paste0(sam_filename_prefix,
                                         chromosomes[i], ".sam"))

        SAM_files <- paste(SAM_files, chr_SAM_file)
    }
    
    name <- paste0(resources_dir, "/out_samtools_merge_", purity, "_", sam_filename_prefix, chromosomes[i])
    cmd <- paste("/usr/bin/time -o", name, "-p -v samtools merge -fc -@", num_of_cores,
		 "-o", BAM_file, SAM_files)
    print(cmd)
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

simulate_seq_resources <- function(c,tumour,ref_path,coverage,purity){
  p_info <- ps::ps_handle()
  start_time <- Sys.time()
  initial_cpu <- ps::ps_cpu_times(p_info)
  initial_mem <- ps::ps_memory_info(p_info)["rss"] / 1024^3
  if (tumour){
    seq_res <- simulate_seq(phylo_forest, reference_genome = ref_path,
                            chromosomes = c,
                            coverage = coverage,
                            purity = purity, 
                            write_SAM = TRUE, 
                            read_size = 150,
                            sequencer = basic_seq,
                            insert_size_mean = 350,
                            insert_size_stddev = 10,
                            output_dir = output_local_dir,
                            include_non_sequenced_mutations = TRUE,
                            update_SAM = TRUE,
                            filename_prefix = sam_filename_prefix,
                            template_name_prefix = paste0(lot_name,'r'),
                            with_normal_sample = FALSE)
  } else{
    cat("simulate_normal_seq")
    seq_res <- simulate_normal_seq(phylo_forest, reference_genome = ref_path,
                                   chromosomes = c,
                                   coverage = coverage,
                                   write_SAM = TRUE, 
                                   read_size = 150,
                                   sequencer = basic_seq,
                                   insert_size_mean = 350,
                                   include_non_sequenced_mutations = TRUE,
                                   insert_size_stddev = 10,
                                   filename_prefix = sam_filename_prefix,
                                   template_name_prefix = paste0(lot_name,'r'),
                                   output_dir = output_local_dir,
                                   with_preneoplastic = with_preneoplastic,
                                   update_SAM = TRUE)
  }

  end_time <- Sys.time()
  final_cpu <- ps::ps_cpu_times(p_info)
  final_mem <- ps::ps_memory_info(p_info)["rss"] / 1024^3
  elapsed_time <- end_time - start_time
  elapsed_time <- as.numeric(elapsed_time, units = "mins")
  cpu_used <- (final_cpu["user"] + final_cpu["system"]) - (initial_cpu["user"]+ initial_cpu["system"])
  mem_used <- final_mem - initial_mem
  resource_usage <- data.frame(
    elapsed_time_mins =  elapsed_time,
    cpu_time_secs = cpu_used,
    memory_used_MB = mem_used,
    chr = c,
    coverage = coverage,
    purity = purity,
    tumour = tumour
  )
  rownames(resource_usage) <- NULL
  out <- list(mutations=seq_res$mutations,parameters=seq_res$parameters,
              resource_usage = resource_usage)
  return(out)
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

resources_dir <- file.path(output_dir, "TIME")
if (!file.exists(resources_dir)) {
  dir.create(resources_dir)
}

data_dir_muts <- file.path(output_dir, "data/mutations")
data_dir_params <- file.path(output_dir, "data/parameters")
data_dir_resources <- file.path(output_dir, "data/resources")

if (!file.exists(data_dir_muts)) {
  dir.create(data_dir_muts,recursive = T)
}

if (!file.exists(data_dir_params)) {
  dir.create(data_dir_params,recursive = T)
}

if (!file.exists(data_dir_resources)) {
  dir.create(data_dir_resources,recursive = T)
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
    basic_seq <- BasicIlluminaSequencer(1e-3) ## only for testing purpose
    chromosomes <- phylo_forest$get_absolute_chromosome_positions()$chr
    cat("parallel::mclapply(",chromosomes,",simulate_seq_resources,mc.cores = 4,tumour = ",
    seq_tumour,",ref_path=",ref_path,",coverage = ",coverage,",purity=",purity,")")
    
    seq_results <- parallel::mclapply(chromosomes, simulate_seq_resources,
                          mc.cores = 4,
                          tumour = seq_tumour, ref_path= ref_path,
                          coverage = coverage, purity = purity)
                          
    seq_results_muts_final <- lapply(1:length(seq_results), function(i) {
        s <- seq_results[[i]]$mutations
        }) %>% do.call(bind_rows, .)


    seq_results_params_final <- lapply(1:length(chromosomes), function(i) {
        pp <- seq_results[[i]]$parameters
        dplyr::tibble(chr = chromosomes[i], parameters = list(pp))
        }) %>% do.call("bind_rows", .)
        
    seq_results_resources_final <- lapply(1:length(seq_results), function(i) {
        s <- seq_results[[i]]$resource_usage
        }) %>% do.call(bind_rows, .)

    saveRDS(seq_results_muts_final,
            file.path(data_dir_muts,
                      paste0("/seq_results_muts_", spn_name,
    			  "_", lot_name, ".rds")))
    saveRDS(seq_results_params_final,
            file.path(data_dir_params,
                      paste0("/seq_results_params_", spn_name,
                  "_", lot_name, ".rds")))
                  
    saveRDS(seq_results_resources_final,
            file.path(data_dir_resources,
                      paste0("/seq_results_resources_", spn_name,
                  "_", lot_name, ".rds")))
    

    cat("done\\n4. Building overall BAM file...")
    merge_sams(output_local_dir, BAM_local_file,
               sam_filename_prefix, chromosomes,
    		   num_of_cores, resources_dir)
    
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

split_bam_by_samples <- function(output_local_dir, BAM_file, remove_local_bam, num_of_cores, resources_dir) {
    name <- strsplit(BAM_file, '/') %>% unlist()
    name <- name[length(name)]
    cmd <- paste0("/usr/bin/time -o ", resources_dir, "/out_samtools_split_", purity,"_",
                  name, " -p -v samtools split -f \\"",
                  file.path(output_local_dir,"%*_%!.bam"),
                  "\\" ", BAM_file, " -@ ", num_of_cores)
    invisible(system(cmd, intern = TRUE))

    if (remove_local_bam) {
        file.remove(BAM_file)
    }
}
invisible(split_bam_by_samples(output_local_dir, BAM_local_file, remove_local_bam, num_of_cores, resources_dir))
## unlink(bam_dir, recursive = TRUE)

cat(paste0("done\\n", step,
           ". Generating the FASTQs and deleting the BAMs..."))
step <- step + 1

BAM_files <- list.files(output_local_dir, pattern = "\\\\.bam$")

generate_fastq <- function(orig_file, fastq_dir, resources_dir) {
  base_orig_file <- tools::file_path_sans_ext(basename(orig_file))

  file_prefix <- file.path(fastq_dir, base_orig_file)
  R1 <- paste0(file_prefix, ".R1.fastq.gz")
  R2 <- paste0(file_prefix, ".R2.fastq.gz")
  unpaired <- paste0(file_prefix, ".unpaired.fastq.gz")
  singleton <- paste0(file_prefix, ".singleton.fastq.gz")

  prefix <- strsplit(file_prefix,'/') %>% unlist()
  prefix <- prefix[length(prefix)]
  name <- paste0(resources_dir, "/out_fastq_", purity,"_",prefix) 
  cmd <- paste("/usr/bin/time -o", name, "-p -v samtools fastq -@ 20 -c 9 -N -1", R1, "-2", R2, "-0", unpaired, 
               "-s", singleton, orig_file)
  invisible(system(cmd, intern = TRUE))
}

result <- parallel::mclapply(BAM_files, function(c) {
    curr_BAM_file <- file.path(output_local_dir, c)
    if (BAM_file != curr_BAM_file) {
        generate_fastq(curr_BAM_file, output_local_dir, resources_dir)

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

sarek_file_launcher="""#!/bin/bash
#SBATCH --partition={PARTITION}
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
nextflow run nf-core/sarek -r 3.5.1 --input $input \
    --outdir $output_dir_combination -profile singularity -c $config
"""

sarek_file_normal_launcher="""#!/bin/bash
#SBATCH --partition={PARTITION}
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
nextflow run nf-core/sarek -r 3.5.1 --input $input \
    --outdir $output_dir_combination --tools haplotypecaller,freebayes -profile singularity -c $config
"""


sarek_variant_calling_launcher="""#!/bin/bash
#SBATCH --partition={PARTITION}
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
nextflow run nf-core/sarek -r 3.5.1 --genome GATK.GRCh38 --input $input \
    --step variant_calling --tools cnvkit,freebayes,strelka,ascat,mutect2 --joint_mutect2 true \
    --outdir $output_dir_combination -profile singularity -c $config
"""

sequenza_launcher="""#!/bin/bash
#SBATCH --partition={PARTITION}
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

nextflow run $base/main.nf -profile singularity --input $input --outdir $output_dir_combination -c $config
"""



tumourevo_launcher="""#!/bin/bash
#SBATCH --partition={PARTITION}
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

nextflow run {TUMOUREVO_BASEDIR}/main.nf --input $input \
    --tools mobster,viber,pyclone-vi,sparsesignatures,sigprofiler \
    --genome GRCh38 \
    --fasta {FASTA_PATH} \
    --download_cache_vep false \
    --vep_cache {VEP_CACHE} \
    --vep_genome GRCh38 \
    --vep_cache_version {VEP_CACHE_VERSION} \
    --vep_species Homo_sapiens \
    --filter false \
    --outdir $output_dir_combination -profile singularity -c $config
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

def write_sarek_sample_variant_calling_lines(sarek_file, SPN, seq_type, sample_name, sarek_dir,coverage,purity):
    if (seq_type == 'tumour'):
        status = 1
    elif (seq_type == 'normal'):
        status = 0
    else:
        raise RuntimeError(f'Unsupported sequence type "{seq_type}"')

    recal_cram_suffix = '.recal.cram'
    recal_crai_suffix = '.recal.cram.crai'

    line = 1
    cram_dir_base_name = f'{coverage}x_{purity}p/preprocessing/recalibrated/{sample_name}/'
    cram_filename = os.path.abspath(os.path.join(sarek_dir,cram_dir_base_name,
                                                sample_name+ recal_cram_suffix))
    crai_filename = os.path.abspath(os.path.join(sarek_dir,cram_dir_base_name,
                                                sample_name+recal_crai_suffix))
    sarek_file.write(f'\n{SPN},{subject_gender},{status},{sample_name},'
                    + f'{cram_filename},{crai_filename}')
    
def write_tumourevo_lines(tumourevo_file, SPN, sample_name, combination, coverage, purity, sarek_output_dir, cancer_type = 'PANCANCER'):
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
    parser.add_argument('-b', '--singularity_binding', type=str, default="",
                       help="The singularity binding in the form /path:/path",required=True)
    parser.add_argument('-TB', '--tumourevo_base_dir', type=str, default="",
                       help="The full path to tumourevo main.nf script, e.g. /path/cloned/tumourevo",required=True)
    parser.add_argument('-Vc', '--vep_cache', type=str, default="",
                       help="The full path to VEP cache",required=True)
    parser.add_argument('-Vv', '--vep_cache_version', type=str, default="110",
                       help="The version of VEP cache",required=True)
    parser.add_argument('-FP', '--fasta_path', type=str, default="",
                       help="The full path to fasta reference genome",required=True)
   

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
    if not os.path.exists(gender_filename):
        with open('ProCESS_subject_gender.R', 'w') as outstream:
            outstream.write(gender_R_script)

        with open('ProCESS_subject_gender.sh', 'w') as outstream:
            outstream.write(gender_shell_script)

        cmd = ['sbatch', '--account={}'.format(account),
            '--partition={}'.format(args.partition),
            ('--export=PHYLO_FOREST={},IMAGE={},DIR={}').format(args.phylogenetic_forest,
                                                args.image_path, curr_dir),
            './ProCESS_subject_gender.sh']

        subprocess.run(cmd)

    with open('ProCESS_seq.R', 'w') as outstream:
        outstream.write(R_script)

    partition = args.partition
    space_per_lot = 3 * cohorts['tumour']['max_coverage'] * 5 / num_of_lots_T
    memory_per_lot = math.ceil(args.mem_per_node*space_per_lot/args.scratch_per_node)
    memory_per_lot = max(memory_per_lot, math.ceil(args.mem_per_node/5))
    shell_script = shell_script.replace('{MEMORY}', str(memory_per_lot))

    with open('ProCESS_seq.sh', 'w') as outstream:
        outstream.write(shell_script)

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    
    sarek_dir = os.path.join(args.output_dir, 'sarek')
    tumourevo_dir = os.path.join(args.output_dir, 'tumourevo')

    tumour_type_file=os.path.join(os.path.dirname(args.phylogenetic_forest),"tumour_type.txt")

    with open('ProCESS_tumour_type.sh', 'w') as outstream:
        outstream.write(write_tumour_type_file)

    cmd = ['sbatch', '--account={}'.format(account),
          '--partition={}'.format(args.partition),
          ('--export=DIR={},SPN={},BASEDIR={},PARTITION={}').format(curr_dir,args.SPN,
                                                os.path.dirname(args.phylogenetic_forest),partition),
          './ProCESS_tumour_type.sh']

    subprocess.run(cmd)

    
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
                        '--job-name={}_{}_{}'.format(args.SPN, purity,lot_name),
                        ('--export=PHYLO_FOREST={},SPN={},LOT={},DEST={},'
                            + 'COVERAGE={},TYPE={},NODE_SCRATCH={},'
                            + 'SEED={},PURITY={},IMAGE={},DIR={},SINGULARITY_BIND={}').format(args.phylogenetic_forest,
                                                args.SPN, lot_name,
                                                output_dir, lot_coverage,
                                                seq_type, args.node_scratch_directory,
                                                i, purity, args.image_path, curr_dir,args.singularity_binding),
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

            if seq_type == 'normal':
                with open(gender_filename, "r") as gender_file:
                    subject_gender = gender_file.read().strip('\n')
                
                normal_fastq_dir = os.path.join(f'{args.output_dir}', 'normal/purity_1/FASTQ')
                lines = math.ceil(math.log10(num_of_lots+1))
                with open(f'{sarek_dir}/sarek_normal.csv', 'w') as sarek_file:
                    sarek_file.write('patient,sex,status,sample,lane,fastq_1,fastq_2')
                    write_sarek_sample_lines(sarek_file, args.SPN, 'normal', 'normal_sample', num_of_lots, normal_fastq_dir, zeros, lines)
                
                job_id=seq_type
                sarek_file_launcher_orig = sarek_file_normal_launcher
                sarek_file_normal_launcher = sarek_file_normal_launcher.replace('{ACCOUNT}', str(account))
                sarek_file_normal_launcher = sarek_file_normal_launcher.replace('{PARTITION}', str(partition))
                sarek_file_normal_launcher = sarek_file_normal_launcher.replace('{JOB_NAME}', str(job_id))
                sarek_file_normal_launcher = sarek_file_normal_launcher.replace('{INPUT_DIR}', str(sarek_dir))
                sarek_file_normal_launcher = sarek_file_normal_launcher.replace('{CONFIG}', str(config_file))
                sarek_file_normal_launcher = sarek_file_normal_launcher.replace('{SAREK_OUT}', str(args.sarek_output_dir))

                with open(f'{sarek_dir}/sarek_mapping_vc_normal.sh', 'w') as outstream:
                    outstream.write(sarek_file_normal_launcher)
                sarek_file_normal_launcher = sarek_file_launcher_orig
                with open('ProCESS_merge_rds.R', 'w') as outstream:
                    outstream.write(merging_R_script)

                with open('ProCESS_merge_rds.sh', 'w') as outstream:
                    outstream.write(merging_shell_script)
                num_of_normal_lots_list = [cohorts['normal']['max_coverage']]

                lots=''
                for l in num_of_normal_lots_list:
                  lots = lots+' '+str(l)
                cmd = ['sbatch', '--account={}'.format(account),
                    '--partition={}'.format(args.partition),
                    '--output={}/merge_normal_{}_{}_{}.log'.format(log_dir,args.SPN, 1,cohorts['normal']['max_coverage']),
                    '--job-name=merge_normal_{}_{}_{}'.format(args.SPN, 1,cohorts['normal']['max_coverage']),
                    ('--export=LOTS_LIST={},SPN={},INPUT_DIR={},PURITY={},TYPE={},IMAGE={},DIR={},MAX_COVERAGE={},TOT_LOTS={},SINGULARITY_BIND={}').format(lots,args.SPN,
                                                        args.output_dir,purity,seq_type,
                                                        args.image_path, curr_dir,cohorts['normal']['max_coverage'],num_of_lots_N,args.singularity_binding),
                    './ProCESS_merge_rds.sh']

                subprocess.run(cmd)
                    
            else:
                with open(gender_filename, "r") as gender_file:
                    subject_gender = gender_file.read().strip('\n')
                    
                tumour_fastq_dir = os.path.join(f'{args.output_dir}', f'tumour/purity_{purity}/FASTQ')
                sample_names = get_sample_names_from_FASTQ(tumour_fastq_dir)
                lines = math.ceil(math.log10(num_of_lots*(len(sample_names)+1)))
                
                for cohort_cov in cohort_coverages:
                    num_of_tumour_lots = math.ceil((cohort_cov*num_of_lots)/cohorts['tumour']['max_coverage'])
                    
                    with open(f'{sarek_dir}/sarek_{cohort_cov}x_{purity}p.csv', 'w') as sarek_file, open(f'{sarek_dir}/sarek_variant_calling_{cohort_cov}x_{purity}p.csv', 'w') as sarek_file_vc:
                        sarek_file.write('patient,sex,status,sample,lane,fastq_1,fastq_2')
                        sarek_file_vc.write('patient,sex,status,sample,cram,crai')
                        cram_normal = os.path.join(f'{args.sarek_output_dir}', f'normal/preprocessing/recalibrated/normal_sample/normal_sample.recal.cram')
                        crai_normal = os.path.join(f'{args.sarek_output_dir}', f'normal/preprocessing/recalibrated/normal_sample/normal_sample.recal.cram.crai')
                        sarek_file_vc.write(f'\n{args.SPN},{subject_gender},0,normal_sample,{cram_normal},{crai_normal}')
                        for sample_name in sample_names:
                            write_sarek_sample_lines(sarek_file, args.SPN, 'tumour', sample_name, num_of_tumour_lots, tumour_fastq_dir, zeros, lines)
                            write_sarek_sample_variant_calling_lines(sarek_file_vc, args.SPN, 'tumour', sample_name, args.sarek_output_dir,cohort_cov,purity)

                    #sarek mapping sh file
                    sarek_file_launcher_orig = sarek_file_launcher    
                    job_id=f'{cohort_cov}x_{purity}p'

                    sarek_file_launcher = sarek_file_launcher.replace('{ACCOUNT}', str(account))
                    sarek_file_launcher = sarek_file_launcher.replace('{PARTITION}', str(partition))
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
                    sarek_variant_calling_launcher = sarek_variant_calling_launcher.replace('{PARTITION}', str(partition))
                    sarek_variant_calling_launcher = sarek_variant_calling_launcher.replace('{JOB_NAME}', str(job_id))
                    sarek_variant_calling_launcher = sarek_variant_calling_launcher.replace('{INPUT_DIR}', str(sarek_dir))
                    sarek_variant_calling_launcher = sarek_variant_calling_launcher.replace('{CONFIG}', str(config_file))
                    sarek_variant_calling_launcher = sarek_variant_calling_launcher.replace('{SAREK_OUT}', str(args.sarek_output_dir))
                    
                    with open(f'{sarek_dir}/sarek_variant_calling_{cohort_cov}x_{purity}p.sh', 'w') as outstream:
                        outstream.write(sarek_variant_calling_launcher)
                    sarek_variant_calling_launcher = sarek_variant_calling_launcher_orig
                    
                    #sequenza sh file
                    sequenza_launcher_orig = sequenza_launcher
                    job_id=f'{cohort_cov}x_{purity}p'
                    process_path = '/'.join(str(config_file).split('/')[:-2]) + '/sequenza'
                    
                    sequenza_launcher = sequenza_launcher.replace('{ACCOUNT}', str(account))
                    sequenza_launcher = sequenza_launcher.replace('{PARTITION}', str(partition))
                    sequenza_launcher = sequenza_launcher.replace('{JOB_NAME}', str(job_id))
                    sequenza_launcher = sequenza_launcher.replace('{INPUT_DIR}', str(sarek_dir))
                    sequenza_launcher = sequenza_launcher.replace('{CONFIG}', str(config_file))
                    sequenza_launcher = sequenza_launcher.replace('{SAREK_OUT}', str(args.sarek_output_dir))
                    sequenza_launcher = sequenza_launcher.replace('{PROCESS_DIR}', str(process_path))

                    
                    with open(f'{sarek_dir}/sequenza_{cohort_cov}x_{purity}p.sh', 'w') as outstream:
                        outstream.write(sequenza_launcher)
                    sequenza_launcher = sequenza_launcher_orig
                    
                    #tumourevo sh file and csv file
                     with open(tumour_type_file) as cancer_type_file:
                         cancer_type = cancer_type_file.read().strip()

                    variant_callers = ['freebayes', 'strelka', 'mutect2']
                    cn_caller = ['ascat', 'sequenza']
                    combinations = []
                    for vc in variant_callers:
                        for cnc in cn_caller:
                            combinations.append([vc, cnc])
                    
                    for comb in combinations:
                        vc = comb[0]
                        cc = comb[1]
                        with open(f'{tumourevo_dir}/tumourevo_{cohort_cov}x_{purity}p_{vc}_{cc}.csv', 'w') as tumourevo_file:
                            if comb[0] == 'mutect2':
                                tumourevo_file.write('dataset,patient,tumour_sample,normal_sample,vcf,tbi,cna_segments,cna_extra,cna_caller,cancer_type')
                            else:
                                tumourevo_file.write('dataset,patient,tumour_sample,normal_sample,vcf,tbi,cna_segments,cna_extra,cna_caller,cancer_type,tumour_alignment,tumour_alignment_index')
                                
                            for sample_name in sample_names:
                                write_tumourevo_lines(tumourevo_file, args.SPN, sample_name, comb, cohort_cov, purity, args.sarek_output_dir,cancer_type)
                    
                        tumourevo_launcher_orig = tumourevo_launcher
                        job_id=f'{cohort_cov}x_{purity}p_{vc}_{cc}'
                        tumourevo_launcher = tumourevo_launcher.replace('{ACCOUNT}', str(account))
                        tumourevo_launcher = tumourevo_launcher.replace('{PARTITION}', str(partition))
                        tumourevo_launcher = tumourevo_launcher.replace('{JOB_NAME}', str(job_id))
                        tumourevo_launcher = tumourevo_launcher.replace('{INPUT_DIR}', str(tumourevo_dir))
                        tumourevo_launcher = tumourevo_launcher.replace('{CONFIG}', str(config_file))
                        tumourevo_launcher = tumourevo_launcher.replace('{TUMOUREVO_OUT}', str(args.tumourevo_output_dir))
                        tumourevo_launcher = tumourevo_launcher.replace('{TUMOUREVO_BASEDIR}', str(args.tumourevo_base_dir))
                        tumourevo_launcher = tumourevo_launcher.replace('{FASTA_PATH}', str(args.fasta_path))
                        tumourevo_launcher = tumourevo_launcher.replace('{VEP_CACHE}', str(args.vep_cache))
                        tumourevo_launcher = tumourevo_launcher.replace('{VEP_CACHE_VERSION}', str(args.vep_cache_version))



                        with open(f'{tumourevo_dir}/tumourevo_{cohort_cov}x_{purity}p_{vc}_{cc}.sh', 'w') as outstream:
                            outstream.write(tumourevo_launcher)
                        tumourevo_launcher = tumourevo_launcher_orig
                    
                with open('ProCESS_merge_rds.R', 'w') as outstream:
                  outstream.write(merging_R_script)

                with open('ProCESS_merge_rds.sh', 'w') as outstream:
                  outstream.write(merging_shell_script)
                
                num_of_tumour_lots_list = [math.ceil((cohort_cov*num_of_lots)/cohorts['tumour']['max_coverage']) for cohort_cov in cohort_coverages]
               
                lots=''
                for l in num_of_tumour_lots_list:
                  lots = lots+' '+str(l)
                
                cmd = ['sbatch', '--account={}'.format(account),
                    '--partition={}'.format(args.partition),
                    '--output={}/merge_tumour_{}_{}_{}.log'.format(log_dir,args.SPN, purity,cohorts['tumour']['max_coverage']),
                    '--job-name=merge_tumour_{}_{}_{}'.format(args.SPN, purity,cohorts['tumour']['max_coverage']),
                    ('--export=LOTS_LIST={},SPN={},INPUT_DIR={},PURITY={},TYPE={},IMAGE={},DIR={},MAX_COVERAGE={},TOT_LOTS={},SINGULARITY_BIND={}').format(lots,args.SPN,
                                                        args.output_dir,purity,seq_type,
                                                        args.image_path, curr_dir,cohorts['tumour']['max_coverage'],num_of_lots_T,args.singularity_binding),
                    './ProCESS_merge_rds.sh']
                subprocess.run(cmd)
