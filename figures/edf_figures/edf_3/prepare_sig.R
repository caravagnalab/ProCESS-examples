library(ggplot2)
library(tidyverse)
library(patchwork)

reticulate::use_python("/orfeo/scratch/area/lvaleriani/myconda/bin")
reticulate::py_config()
library(SigProfilerMatrixGeneratorR)
library(SigProfilerAssignmentR)

setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/')
source('getters/tumourevo_getters.R')
source('getters/process_getters.R')
source('figures/figure3/utils_plot.R')

base = "/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/"


SPN_colors <-c("SPN01"='steelblue', "SPN02"='seagreen', "SPN03"='goldenrod', 
               "SPN04"='coral', "SPN05"="magenta4","SPN06"='palevioletred', "SPN07"='indianred3')

cov = 50
pur = 0.9


for (spn in names(SPN_colors)){
  data <- readRDS(get_mutations(spn = spn, coverage = cov, purity = pur, type = 'tumour')) 
  
  somatic <- data %>%
    dplyr::filter(classes!="germinal") %>%
    dplyr::filter(!stringr::str_detect(causes, 'errors')) %>%
    ProCESS::seq_to_long()
  
  mut_class <- somatic %>% 
    ungroup() %>%
    filter(DP != 0) %>%
    filter(NV != 0) %>%
    mutate(mut_id = paste(chr,  from, ref, alt, sep = ':')) %>% 
    select(mut_id, causes)
  
  sig_data <- somatic %>%
      ungroup() %>%
      filter(DP != 0) %>%
      filter(NV != 0) %>%
      select(chr,  from, ref, alt, sample_name) %>%
      dplyr::mutate(Project = 'spn', Genome = 'GRCh38', mut_type = 'SNP', Type = 'SOMATIC', ID = sample_name, Sample = sample_name) %>%
      dplyr::rename(chrom = chr, pos_start = from) %>%
      filter(ref != alt) %>%
      rowwise() %>%
      dplyr::mutate(pos_end = pos_start + abs(str_count(ref) - str_count(alt))) %>%
      dplyr::select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type) %>%
      distinct()
  
  dir.create(paste0('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/sign/', spn), showWarnings = F, recursive = T)
  write.table(sig_data, file = paste0('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/sign/', spn, '/mutations.txt'), quote = F, sep = '\t', row.names = F)
  SigProfilerMatrixGeneratorR(project = spn, seqInfo = T, genome = "GRCh38", matrix_path = paste0('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/sign/', spn, '/'), plot=T)
  
  saveRDS(mut_class, paste0('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/sign/', spn, '/mut_class.rds'))

}

