
rm(list = ls())
options(bitmapType='cairo')
require(tidyverse)
library(optparse)

source("../../getters/process_getters.R")
source("utils/utils.R")
source("utils/plot_utils.R")

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN04'),
                    make_option(c("--purity"), type = "character", default = '0.6'),
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
# chromosomes = paste0("chr", c(1:22, "X", "Y"))
min_vaf = .02
mut_types = c("INDEL", "SNV")
comb = list(PI = purity, COV = coverage)
#path_to_cna <- paste0(data_dir,spn_id,"/process/cna_data/")
#samples = gsub(pattern="_cna.rds",replacement="",x=list.files(path_to_cna,pattern="_cna.rds"))
samples = get_sample_names(spn = spn_id)

# path_to_seq = get_mutations(spn = spn_id, 
#                         base_path = data_dir, 
#                         coverage = coverage, 
#                         purity = purity, 
#                         type = "tumour")

#path_to_seq <- paste0(data_dir,spn_id,"/sequencing/tumour/purity_",purity,"/data/mutations/seq_results_muts_merged_coverage_",coverage,"x.rds")

input_dir <-  paste0(data_dir,spn_id,"/validation/somatic/")
outdir <- paste0(data_dir,spn_id,"/validation/somatic/report")
dir.create(outdir, recursive = T, showWarnings = F)

# Preparing report
message("Parsing combination: purity=", purity, ", cov=", coverage)
#sample_id = samples[1]
mut_type = "INDEL"

for (mut_type in mut_types) {
  message("      Considering only ", mut_type)
  combination = paste0(coverage, "x_", purity, "p")
  
  # Single lapply to get both ground truth and caller results
  results = lapply(samples, function(sample_id) {
    # Get ground truth
    process_folder_path <- file.path(input_dir, spn_id, combination, "process", sample_id, mut_type)  
    gt_res = lapply(chromosomes, function(chromosome){
      gt_path = file.path(process_folder_path, paste0(chromosome,".rds"))
      readRDS(gt_path)
    }) %>% do.call("bind_rows", .) %>% 
      dplyr::mutate(sample = sample_id) %>% 
      dplyr::mutate(mutationID = paste0(mutationID,":",sample_id))
    
    # Get caller results for this sample
    caller_res_list = lapply(callers, function(caller) {
      print(paste("Processing caller:", caller, "for sample:", sample_id))
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
    })
    names(caller_res_list) = callers
    
    # Return both ground truth and caller results for this sample
    list(gt = gt_res, callers = caller_res_list)
  })
  names(results) = samples
  
  # Extract ground truth and caller results
  gt_res = lapply(results, function(x) x$gt) %>% do.call("bind_rows", .)
  caller_res_list = lapply(callers, function(caller) {
    lapply(results, function(x) x$callers[[caller]]) %>% do.call("bind_rows", .)
  })
  names(caller_res_list) = callers
  
  sample_info = list(mut_type=mut_type, spn=spn_id, purity=purity, coverage=coverage)
  report_plot = get_multi_caller_report(seq_res_long = gt_res, 
                                   caller_res_list = caller_res_list, pi = as.numeric(purity),
                                   sample_info = sample_info, 
                                   min_vaf = min_vaf, 
                                   only_pass = TRUE)
  
  filename = paste(spn_id, combination, "allCaller", mut_type, sep = '_')
  file_path = file.path(outdir, filename)
  ggsave(paste0(file_path, '.png'), plot = report_plot, width = 15, height = 20, units = "in", dpi = 400)
  
  # Get Metrics multi-caller by sample
  
  metrics = lapply(names(caller_res_list), function(nc) {
    caller_res = caller_res_list[[nc]]
    
    samples_results = lapply(samples, function(sample_id) {
      sample_ground_truth_res = gt_res %>% dplyr::filter(sample == sample_id)
      sample_caller_res = caller_res %>% dplyr::filter(sample == sample_id)
      
      metrics_results = analyze_vaf_performance(sample_ground_truth_res, 
                                                sample_caller_res, 
                                                only_pass = TRUE, 
                                                min_vaf_threshold = min_vaf, 
                                                vaf_tolerance_pct = 5)  # 5% VAF tolerance
      
      metrics_results = analyze_ccf_performance(sample_ground_truth_res, 
                                                sample_caller_res, 
                                                pi = as.numeric(purity),
                                                only_pass = TRUE, 
                                                min_ccf_threshold = min_vaf, 
                                                ccf_tolerance_pct = 5)  # 5% VAF tolerance
      
      metrics_results$raw_data = NULL
      metrics_results
    })
    names(samples_results) = samples
    samples_results
  })
  names(metrics) = names(caller_res_list)
  
  results_folder_path = file.path(input_dir, spn_id, combination, "allCaller", mut_type)
  dir.create(results_folder_path, recursive = T)
  metrics_path = file.path(results_folder_path, "metrics.rds")
  saveRDS(metrics, metrics_path)
  
  message("        Report done!")
}
