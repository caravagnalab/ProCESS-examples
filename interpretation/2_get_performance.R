setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")

library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
library(ggalluvial)
library(purrr)
library(stringr)
library(ggtext)  

df_process <- readRDS('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/process.rds')
coverage_list = c(50, 100, 150)
purity_list = c(0.9, 0.6, 0.3)
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04','SPN05', 'SPN06', 'SPN07')
tool_list = c( 'viber', 'pyclonevi')
sig_tool_list = c('BASCULE', 'SigProfiler')
combs = expand.grid(spn = spn_list,
                    coverage=coverage_list,
                    purity=purity_list,
                    tool=tool_list,
                    sig_tool = sig_tool_list)

option_list <- list(make_option(c("--cna_caller"), type = "character", default = 'sequenza'),
                    make_option(c("--vcf_caller"), type = "character", default = 'mutect2')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

vcf_caller = opt$vcf_caller #"mutect2"
cna_caller = opt$cna_caller #"sequenza"
base = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/interpretation_', vcf_caller, "_", cna_caller, '/')

i=1
results <- lapply(1:nrow(combs), FUN = function(i){
  sp = combs[i, "spn"]
  cov = combs[i, "coverage"]
  pur = combs[i, "purity"]
  s_tool = combs[i, "sig_tool"]
  tool = combs[i, "tool"]

  data_tool = paste0(base, '/res_int/', cov, 'x_', pur, 'p_', s_tool,'_',tool)
  df_tool = readRDS(file = paste0(data_tool, '_df.rds')) %>%
    filter(spn == sp) %>%
    select(cluster, contains_driver_process, score_all, score_driver, score_tail, score_sign, score_no_driver, score_no_tail, score_no_sign, driver_label_process, is_clonal_tool)
  
  if (nrow(df_tool) > 0){
  
    tmp_process = df_process %>%
      filter(spn == sp, coverage == cov, purity == pur) %>%
      select(cluster, score_all, score_driver, score_tail, score_sign, score_no_driver, score_no_tail, score_no_sign, driver_label_process) %>%
      filter(cluster != 'Subclonal')
  
    join = left_join(tmp_process, df_tool,
                     by = join_by(driver_label_process),
                     suffix = c('_process', '_tool')) %>%
      mutate(spn = sp, coverage = cov, purity = pur, sig_tool = s_tool, sub_tool = tool)
  }

}) %>%
  bind_rows()

saveRDS(file = paste0(base,'/performance_cluster.rds'),
        object = results)



coverage_list = c(50, 100, 150)
purity_list = c(0.9, 0.6, 0.3)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04','SPN05', 'SPN06', 'SPN07')
tool_list = c( 'viber', 'pyclonevi')
sig_tool_list = c('BASCULE', 'SigProfiler')
combs = expand.grid(spn = spn_list, 
                    coverage=coverage_list,
                    purity=purity_list,
                    tool=tool_list,
                    sig_tool = sig_tool_list)

i=1

mutations <- lapply(1:nrow(combs), FUN = function(i){
  sp = combs[i, "spn"]
  cov = combs[i, "coverage"]
  pur = combs[i, "purity"]
  s_tool = combs[i, "sig_tool"]
  tool = combs[i, "tool"]
  print(paste(sp, cov, pur, s_tool, tool))
  
  df_process <- readRDS('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/process.rds') %>% 
    filter(sp == spn, 
           coverage == cov, 
           purity == pur) %>% 
    select(cluster, score_all, score_driver, score_tail, score_sign, score_no_driver, score_no_tail, score_no_sign) %>% 
    dplyr::rename(cluster_id_process = cluster)
  
  data_tool = paste0(base, '/res_int/', cov, 'x_', pur, 'p_', s_tool,'_',tool)
  df_tool = readRDS(file = paste0(data_tool, '_df.rds')) %>% 
    filter(spn == sp) %>% 
    select(cluster, score_all, score_driver, score_tail, score_sign, is_clonal_tool, score_no_driver, score_no_tail, score_no_sign) %>% 
    distinct()  %>% 
    dplyr::rename(cluster_id_tool = cluster)
  
  if (nrow(df_tool) > 0){
    if (file.exists(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal_new/tables_interpreted/', tool, '_', sp, '_', cov, 'x_', pur, 'p_', vcf_caller, '_',cna_caller, '.rds'))){
      mutations <-  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal_new/tables_interpreted/', tool, '_', sp, '_', cov, 'x_', pur, 'p_', vcf_caller, '_',cna_caller, '.rds'))
      
      if (tool == 'pyclonevi'){
        mutations <- mutations %>% mutate(cluster_id_tool = paste0('C', cluster_id_tool))
      }
      
      join <- mutations %>% 
        select(mutation_id, cluster_id_tool, cluster_id_process) %>% 
        left_join(df_tool) %>% 
        left_join(df_process, by = join_by(cluster_id_process), suffix = c('_tool', '_process')) %>% 
        distinct() %>% 
        mutate(spn = sp, coverage = cov, purity = pur, sig_tool = s_tool, sub_tool = tool) %>% 
        filter(!is.na(score_all_tool))
      
      return(join)
    }
  }
}) %>% 
  bind_rows()
saveRDS(file = paste0(base, '/performance_mutations.rds'), 
        object = mutations)
