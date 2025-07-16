rm(list = ls())

setwd("/orfeo/cephfs/scratch/cdslab/gsantacatterina/ProCESS-examples/validation/Somatic")

require(tidyverse)
library(optparse)
source("utils/mutect_utils.R")
source("utils/process_utils.R")
source("utils/freeBayes_utils.R")
source("utils/strelka_utils.R")
source("../../getters/process_getters.R")
source("../../getters/sarek_getters.R")

get_missing_files <- function(spn,purity,coverage,caller, chr,base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"){
  missing_files <- list()
  samples <- get_sample_names(spn)
  for (sample in samples){
    for (type in c("SNV","INDEL")){
      expected_file_name <- file.path(base_path,spn,"validation/somatic/",spn,paste0(coverage,"x_",purity,"p"),caller,sample,type,paste0("chr",chr,".rds"))
      # print(expected_file_name)
      if (!file.exists(expected_file_name)){
        missing_files <- c(missing_files,expected_file_name)
        # print(expected_file_name)
      }
    }    
  }
  return(missing_files)
}

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN01'),
                    make_option(c("--purity"), type = "character", default = '0.6'),
                    make_option(c("--coverage"), type = "character", default = '100'),
                    make_option(c("--chr"), type = "character", default = '22')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
data_dir = '/orfeo/scratch/cdslab/shared/SCOUT'
spn_id = opt$spn_id
coverage = opt$coverage
purity = opt$purity
chromosome <- opt$chr
outdir <-  file.path(data_dir,spn_id,"/validation/somatic/")

gt_path = get_mutations(spn = spn_id, 
                        base_path = data_dir, 
                        coverage = coverage, 
                        purity = purity, 
                        type = "tumour")


if (length(get_missing_files(spn_id,coverage=coverage,purity=purity,caller="process",chromosome))>0){
  process_seq_results(gt_path, spn_id, purity, coverage, chromosome, base_path = data_dir, outdir = outdir)
  message("PARSING ProCESS")
} else {
  message("ProCESS preprocessing already run!")
}

if (length(get_missing_files(spn_id,coverage=coverage,purity=purity,caller="mutect2",chromosome))>0){
  process_mutect2_results(gt_path, spn_id, purity, coverage, chromosome, base_path = data_dir, outdir = outdir)
  message("PARSING MUTECT2")
} else {
  message("MUTECT2 preprocessing already run!")
}

if (length(get_missing_files(spn_id,coverage=coverage,purity=purity,caller="freebayes",chromosome))>0){
  process_freebayes_results(gt_path, spn_id, purity, coverage, chromosome, base_path = data_dir, outdir = outdir, pass_quality = 20, min_vaf = 0.01, max_normal_vaf = 0.02)
  message("PARSING FREEBAYES")
} else {
  message("FREEBAYES preprocessing already run!")
}

if (length(get_missing_files(spn_id,coverage=coverage,purity=purity,caller="strelka",chromosome))>0){
  process_strelka_results(gt_path, spn_id, purity, coverage, chromosome, base_path = data_dir, outdir = outdir)
  message("PARSING STRELKA")
} else {
  message("STRELKA preprocessing already run!")
}
