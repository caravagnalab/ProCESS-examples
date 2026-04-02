setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')

out = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/assing_signature/"

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN01')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

base = '/orfeo/cephfs/scratch/cdslab/shared/SCOUT/validation_subclonal_new/tables/'
spn = opt$spn_id

get_exposure <- function(table){
  
  nmuts <- table %>% 
    group_by(patient_id, coverage, purity, cluster_id_process) %>% 
    summarise(
      tot_nmuts = n_distinct(mutation_id),
      .groups = "drop"
    )
  
  exposure_raw <- table %>% 
    group_by(patient_id, coverage, purity, causes, cluster_id_process) %>% 
    summarise(
      nmuts_cause = n_distinct(mutation_id),
      .groups = "drop"
    )
  
  exposure <- exposure_raw %>% 
    left_join(
      nmuts,
      by = c("patient_id", "coverage", "purity", "cluster_id_process")
    ) %>% 
    mutate(
      exposure = nmuts_cause / tot_nmuts
    )
  return(exposure)
}

cov = 100
pur = .9
for (cov in c(50,100,150)){
  print(cov)
  
  for (pur in c(0.3, 0.6, 0.9)){
    print(pur)
    
    out_data = paste0(out,  spn, '/process_univariate/', cov, 'x_', pur, 'p')
    dir.create(out_data, recursive = T, showWarnings = F)
    
    table <- readRDS(paste0(base, 'table_process_univariate_w_private_', spn, '_', cov, 'x_', pur, 'p_mutect2_ascat.rds')) %>% 
      filter(is_driver_process != TRUE) %>% 
      filter(!str_detect(causes, "errors")) 
  
    # correct name of clusters
    table <- table %>% 
      group_by(cluster_id_process, sample_id) %>%
      mutate(is_clonal_process=replace(FALSE, ccf_process > 0.95, TRUE)) %>% 
      ungroup() %>%
      mutate(cluster_id_process_full = cluster_id_process) %>%
      mutate(cluster_id_process = replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal')) %>% 
      separate(cluster_id_process, sep = '_', into = c('cluster_id_process', 'tmp')) %>% 
      select(-tmp)
    
    for (s in unique(table$sample_id)){
      sample = sub("^(.+)_\\1_", "\\1_", s)
        
      n_muts <- table %>% 
        filter(sample_id == s) %>% 
        group_by(cluster_id_process) %>% 
        summarise(n = n()) %>% 
        filter(n>=100)
      
      table_sbs <- table %>% 
        filter(sample_id == s) %>% 
        filter(str_detect(causes, "SBS"))
      
      table_id <- table %>% 
        filter(sample_id == s) %>% 
        filter(str_detect(causes, "ID"))
      
      
      exposure_sbs <- get_exposure(table_sbs) %>% mutate(sample_id = sample)
      exposure_id <- get_exposure(table_id) %>% mutate(sample_id = sample)
      
      saveRDS(object = exposure_sbs, paste0(out_data, '/', sample, '_exposure_SBS.rds'))
      saveRDS(object = exposure_id, paste0(out_data, '/', sample, '_exposure_ID.rds'))
    }
  }
}
