
rm(list = ls())
options(bitmapType='cairo')
require(tidyverse)
library(optparse)

source("../../getters/process_getters.R")
source("utils/utils.R")
source("utils/plot_utils.R")

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN04'),
                    make_option(c("--purity"), type = "character", default = '0.3'),
                    make_option(c("--coverage"), type = "character", default = '50'))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
data_dir = '/orfeo/scratch/cdslab/shared/SCOUT/'
spn_id = opt$spn_id
coverage = opt$coverage
purity = opt$purity

# gender_file <- read.table(file = paste0(data_dir,spn_id,"/process/subject_gender.txt"),header = FALSE,col.names = "gender")
# gender <- gender_file$gender
gender = get_process_gender(spn = spn_id)
if (gender=="XX"){
  chromosomes = c(paste0('chr',1:22), 'chrX')
} else {
  chromosomes = c(paste0('chr',1:22), 'chrX', 'chrY')
}

# INPUT PARAMATERS ####
callers = c("mutect2", "strelka", "freebayes")
min_vaf = .02
mut_types = c("INDEL", "SNV")
comb = list(PI = purity, COV = coverage)
samples = get_sample_names(spn = spn_id)

input_dir <-  paste0(data_dir,spn_id,"/validation/somatic/")
outdir <- paste0(data_dir,spn_id,"/validation/somatic/report")
dir.create(outdir, recursive = T, showWarnings = F)

# Preparing report
message("Parsing combination: purity=", purity, ", cov=", coverage)
caller = "mutect2"
sample_id = samples[1]
mut_type = "SNV"
caller =
for (caller in callers) {
  message("  Using caller: ", caller)
  for (mut_type in mut_types) {
    message("      Considering only ", mut_type)
    
    # spn <- gsub(".*SCOUT/(SPN[0-9]+).*", "\\1", path_to_seq)
    # purity <- gsub(".*purity_([0-9.]+).*", "\\1", path_to_seq)
    # coverage <- gsub(".*coverage_([0-9]+).*", "\\1", path_to_seq)
    combination = paste0(coverage, "x_", purity, "p")
    
    gt_res = lapply(samples, function(sample_id) {
      process_folder_path <- file.path(input_dir, spn_id, combination, "process", sample_id, mut_type)  
      # Get ground truth
      gt_res = lapply(chromosomes, function(chromosome){
        gt_path = file.path(process_folder_path, paste0(chromosome,".rds"))
        readRDS(gt_path)
      }) %>% do.call("bind_rows", .) %>% 
        dplyr::mutate(sample = sample_id) %>% 
        dplyr::mutate(mutationID = paste0(mutationID,":",sample_id))
      gt_res
    }) %>% do.call("bind_rows", .)
    
    caller_res = lapply(samples, function(sample_id) {
      # Get caller res
      caller_folder_path = file.path(input_dir, spn_id, combination, caller, sample_id, mut_type)
      caller_res = lapply(chromosomes, function(chromosome) {
        caller_path = file.path(caller_folder_path, paste0(chromosome,".rds"))
        readRDS(caller_path) %>% 
          dplyr::filter(FILTER == "PASS")
      }) %>% do.call("bind_rows", .) %>% 
        dplyr::filter(!is.na(VAF)) %>% 
        dplyr::mutate(sample = sample_id) %>% 
        dplyr::mutate(mutationID = paste0(mutationID,":",sample_id))
      caller_res
    }) %>% do.call("bind_rows", .)
    
    sample_info = list(caller_name=caller, mut_type=mut_type, spn=spn_id, purity=purity, coverage=coverage)
    report = get_report(seq_res_long = gt_res, 
                        caller_res = caller_res, 
                        sample_info = sample_info, 
                        min_vaf = min_vaf)
    
    results_folder_path = file.path(input_dir, spn_id, combination, caller, mut_type)
    dir.create(results_folder_path, recursive = T)
    metrics_path = file.path(results_folder_path, "metrics.rds")
    saveRDS(list(report_metrics=report$report_metrics, vaf_comparison=report$vaf_comparison), metrics_path)
    
    filename = paste(spn_id, combination, caller, mut_type, sep = '_')
    file_path = file.path(outdir, filename)
    ggsave(paste0(file_path, '.png'), plot = report$report_plot, width = 18, height = 18, units = "in", dpi = 400)
    message("        Report done!")
  }
}
