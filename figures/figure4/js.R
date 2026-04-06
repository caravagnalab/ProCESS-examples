setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")

library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
library(rstatix)

source('../../validation/SCOUT/colors.R')
source("../figure3/utils_plot.R")
source("../figure3/utils.R")

df_jaccard <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature/analysis/js_clonal.rds')


plt_jaccard <- df_jaccard %>% 
  mutate(tool = case_when(
    tool =='mobster_univariate' ~ 'MOBSTER',
    tool =='pyclonevi' ~ 'PyClone-VI' ,
    tool =='viber' ~ 'VIBER'
  )) %>% 
  mutate(mut_caller = case_when(
    mut_caller =='mutect2' ~ 'Mutect2',
    mut_caller =='strelka' ~ 'Strelka'
  )) %>% 
  mutate(cna_caller = case_when(
    cna_caller =='ascat' ~ 'ASCAT',
    cna_caller =='sequenza' ~ 'Sequenza'
  )) %>% 
  ggplot() +
  geom_boxplot(aes(x = mut_caller, y = jaccard, col = tool, fill= tool), alpha = .5, outliers = F, show.legend = T) + 
  scale_color_manual(
    "Tool",
    values = c('MOBSTER' = "#219ebc", 'PyClone-VI' =  "palegreen4", 'VIBER' = "tomato2")
  )  +
  scale_fill_manual(
    "Tool",
    values = c('MOBSTER' = "#219ebc", 'PyClone-VI' =  "palegreen4", 'VIBER' = "tomato2")
  )  +
  my_ggplot_theme() +
  ylab('Clonal cluster mutation assignment\nJaccard index') +
  xlab('Tool') #+
  #facet_grid(.~cna_caller)

