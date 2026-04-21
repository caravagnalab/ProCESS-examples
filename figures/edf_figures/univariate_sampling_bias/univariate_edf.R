setwd("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/univariate_sampling_bias/")
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
source('../../../getters/tumourevo_getters.R')
source('../../../getters/process_getters.R')
source("../../../validation/SCOUT/colors.R")

indir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/assign_signature/"
indir_process= "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/assign_signature/"


colors_cluster = c('indianred', 
                   'steelblue', 
                   'palevioletred',                   
                   'goldenrod', 
                   'darkorange3', 
                   'palevioletred', 
                   'mediumpurple', 
                   'cornsilk4', 
                   'olivedrab3', 
                   'steelblue4', 
                   'indianred4',
                   'aquamarine3',
                   'saddlebrown',
                   'deeppink2',
                   'cornflowerblue',
                   'black')
names(colors_cluster) = paste0('C',0:15)
tail = c('Tail' = 'gray70')
colors_cluster = c(colors_cluster, tail)

color_palette_process_uni = RColorBrewer::brewer.pal(n = 8, name = "Set1") 
names(color_palette_process_uni) <- c('Clonal', paste0('Clone ', 1:7))
color_palette_process_uni['Subclonal'] = 'gray70'

examples = tibble(spn = c('SPN06', 'SPN04', 'SPN07', 'SPN03'),
                  sample = c('SPN06_2.1', 'SPN04_2.1', 'SPN07_1.2', 'SPN03_2.1'),
                  mut_caller = c('mutect2','mutect2', 'mutect2', 'mutect2' ),
                  cna_caller =  c('ascat','ascat','sequenza', 'ascat' ),
                  title = c('Monoclonal', 'Monoclonal-Sampling Bias', 'Polyclonal',
                            'Polyclonal-Weak Evidence'),
                  coverage = c(100, 100, 100, 150))

cluster_plots_samples <- list()
tool = 'mobster_univariate'
for (i in 1:nrow(examples)){
  tmp = examples[i,]
  spn = tmp$spn
  s = tmp$sample
  sample = paste0(spn, '_', tmp$sample)
  mut_caller = tmp$mut_caller
  cna_caller = tmp$cna_caller
  cov = tmp$coverage
  pur = 0.9
  name = tmp$title
  
  table_univariate = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal_new/tables_interpreted/', tool, '_', spn, '_', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '.rds')
  mutations <-  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal_new//tables_interpreted/', tool, '_', spn, '_', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '.rds')) %>% 
    separate(cluster_id_process, sep = '_', into = c('cluster_id_process', 'tmp')) %>% 
    select(-tmp)

  mutations <- mutations %>% mutate(driver_label = sub("_.*", "", driver_label_tool))
  driver_tool <- mutations %>% select(patient_id, mutation_id, is_driver_tool, driver_label, is_driver_process)
  
  id_muts <- mutations %>% 
    dplyr::summarise(n = dplyr::n(), .by = c(patient_id, mutation_id, cluster_id_tool,
                                             sample_id)) |>
    dplyr::filter(n > 1L) 
  
  muts_tool_raw <- mutations %>% 
    filter(!mutation_id %in% id_muts$mutation_id) %>% 
    select(patient_id, mutation_id, cluster_id_tool, sample_id, vaf_tool) %>% 
    pivot_wider(values_from = vaf_tool, names_from = sample_id) %>% 
    replace(is.na(.), 0) %>% 
    mutate(cluster_id_tool = as.character(cluster_id_tool))
  muts_tool_raw <- muts_tool_raw %>% left_join(driver_tool) %>% distinct()
  
  muts_process <- mutations %>% 
    filter(!mutation_id %in% id_muts$mutation_id) %>% 
    select(patient_id, mutation_id, sample_id, vaf_process, cluster_id_process) %>%
    pivot_wider(values_from = vaf_process, names_from = sample_id) %>% 
    replace(is.na(.), 0) 
  
  driver_process <- mutations %>%   filter(!mutation_id %in% id_muts$mutation_id) %>% select(patient_id, mutation_id, is_driver_process, driver_label)
  muts_process <- muts_process %>% left_join(driver_process) %>% distinct()
  

  color_palette_process = RColorBrewer::brewer.pal(n = max(3,length(unique(muts_process$cluster_id_process))), name = "Dark2") %>%
    setNames(str_sort(unique(muts_process$cluster_id_process), numeric=T))
  color_palette_process['Subclonal'] = 'gray70'
  
  p_process = muts_process %>%
    filter(.data[[sample]] != 0) %>% 
    ggplot(aes(x =.data[[sample]], fill = cluster_id_process))+
    geom_histogram( alpha=0.8, binwidth = 0.01, position = "identity")+
    xlim(-0.01,1.01)+
    xlab(paste0('VAF ', s)) +
    theme_bw() +
    scale_fill_manual('ProCESS clusters', values = color_palette_process_uni)+
    ggtitle(subtitle = 'ProCESS', label = name) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))  +
    ylab('')+
    theme(text         = element_text(size = 9), legend.position = 'none')
  
  p_tool_raw = muts_tool_raw %>%
    filter(.data[[sample]] != 0) %>% 
    ggplot(aes(x =.data[[sample]], fill=cluster_id_tool))+
    geom_histogram( alpha=0.8, binwidth = 0.01, position = "identity")+
    xlim(-0.01,1.01)+
    xlab(paste0('VAF ', s)) +
    theme_bw() +
    ylab('') + 
    scale_fill_manual('MOBSTER clusters', values = colors_cluster)+
    theme_bw() +
    ggtitle(label = '', subtitle = 'MOBSTER') +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme(text         = element_text(size = 9), legend.position = 'none') 
    
  cluster_plots_samples[[sample]] <- p_process + p_tool_raw + plot_layout(guides = 'collect')
}

# ggsave(plot = all, 
#        filename = 'samplign_bias.pdf',
#        width = 4.7, 
#        height = 6.7)

