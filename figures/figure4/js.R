library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
library(rstatix)

source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/SCOUT/colors.R')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils.R")

df = tibble()
for (spn in paste0('SPN0', 1:7)){
  for (cov in c(50, 100, 150)){
    for (pur in c(0.3, 0.6, 0.9)){
      for (signature_tool in c('SigProfiler', 'BASCULE')){
        for (tool in c('viber', 'pyclonevi')){
          data = paste0('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature/assign_res/',spn,'/',spn,'_',cov, 'x_', pur, 'p_',tool,'_',signature_tool,'.rds')  
          if (file.exists(data)){
            res = readRDS(data)
            df = bind_rows(df, res)
          }
        }
      }
    }
  }
}

df_uni = tibble()
for (spn in paste0('SPN0', 1:7)){
  for (cov in c(50, 100, 150)){
    for (pur in c(0.3, 0.6, 0.9)){
      for (signature_tool in c('SigProfiler', 'BASCULE')){
        tool = 'mobster_univariate'
        data = paste0('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature/assign_res_univariate/',spn,'/',spn,'_',cov, 'x_', pur, 'p_',tool,'_',signature_tool,'.rds')  
        if (file.exists(data)){
          res = readRDS(data)
          df_uni = bind_rows(df_uni, res)
        }
      }
    }
  }
}


classes <- readRDS('/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure5/signatures_cohort/sample_classification_clonal_heterogeneity.rds') %>% 
  select(patient_id, sample, coverage, purity, sample_class) %>% 
  distinct()


df_uni = df_uni %>% 
  filter(!is.na(CosineSimilarity)) %>% 
  mutate(cluster_process = ifelse(cluster_process == 'Clonal', 'Clonal', 'Subclonal')) %>% 
  filter(!is.na(cluster_process)) 

df_uni = df_uni  %>% 
  left_join(classes %>% dplyr::rename(pur = purity, cov = coverage))


df_uni_filt <- df_uni %>% filter(cluster != 'Tail-Subclonal')
df = df %>% 
  bind_rows(df_uni_filt) %>% 
  filter(!is.na(CosineSimilarity)) %>% 
  mutate(cluster_process = ifelse(cluster_process == 'Clonal', 'Clonal', 'Subclonal')) %>% 
  filter(!is.na(cluster_process)) %>% 
  mutate(patient_class = ifelse(spn %in% c('SPN03','SPN04', 'SPN06', 'SPN07'), 'High\nComplexity', 'Low\nComplexity')) %>% 
  filter(!(pur == 0.6 & spn == "SPN07"))

df_jaccard <- df %>% 
  select(cluster, cluster_process, cov, pur, tool, spn, jaccard, sample) %>% 
  distinct() 


plt_jaccard <- df_jaccard %>% 
  mutate(tool = case_when(
    tool =='mobster_univariate' ~ 'MOBSTER',
    tool =='pyclonevi' ~ 'PyClone-VI' ,
    tool =='viber' ~ 'VIBER'
  )) %>% 
  ggplot() +
  geom_boxplot(aes(x = cluster_process, y = jaccard, col = cluster_process, fill= cluster_process), alpha = .5, outliers = F) + 
  scale_color_manual(
    "Cluster Type",
    values = c('Clonal' = "palegreen4", 'Subclonal' =  "#9F4576")
  )  +
  scale_fill_manual(
    "Cluster Type",
    values = c('Clonal' = "palegreen4", 'Subclonal' =  "#9F4576")
  )  +
  my_ggplot_theme() +
  facet_grid(.~tool) +
  ylab('Mutation assignment (Jaccard index)')
