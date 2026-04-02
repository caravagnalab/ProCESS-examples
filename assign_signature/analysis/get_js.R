setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
library(rstatix)
source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')
source("../validation/SCOUT/colors.R")

base = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/validation_subclonal_new/tables_interpreted/"


cna_caller = c('ascat', 'sequenza')
mut_caller = c('mutect2', 'strelka')

spn = 'SPN02'
cov = 100
pur = 0.9
tool = 'pyclonevi'

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

jaccard_df <- tibble()
spn_list = paste0('SPN0', c(1,2,3,4,5,6,7))
for (vcf in mut_caller){
  for (cna in cna_caller){
    for (sp in spn_list){
      for (cov in c(50, 100, 150)){
        for (pur in c(0.3, 0.6, 0.9)){
          for (tl in c('viber', 'pyclonevi', 'mobster_univariate')){
            
            file = paste0(base, tl, '_', sp, '_', cov, 'x_', pur, 'p_',vcf, '_', cna, '.rds')
            if (file.exists(file)){
              join = readRDS(file)
              
              p = join %>% filter(is_clonal_process == T) %>% pull(mutation_id)
              t = join %>% filter(is_clonal_tool == T) %>% pull(mutation_id)
              
              jaccard_df <- bind_rows(jaccard_df, tibble(spn = sp,
                                                         coverage = cov,
                                                         purity = pur,
                                                         tool = tl,
                                                         mut_caller = vcf,
                                                         cna_caller = cna, 
                                                         jaccard = jaccard(p, t)))
            }
          }
        }
      }
    }
  }
}


saveRDS(object = jaccard_df, file = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature/analysis/js_clonal.rds')

