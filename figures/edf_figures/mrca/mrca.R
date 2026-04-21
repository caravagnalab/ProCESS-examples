setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils.R")
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/SCOUT/colors.R')

color_sbs <- c('SBS25' = 'sienna2', 'SBS31' = 'maroon' )

color_palette_process_uni = RColorBrewer::brewer.pal(n = 8, name = "Dark2") 
names(color_palette_process_uni) <- c('Clonal', paste0('Clone ', 1:7))
color_palette_process_uni['Other'] = 'gray70'

relative_to_absolute <- function(cna){
  reference_genome = CNAqc:::get_reference('GRCh38', data = x$genomic_coordinates)
  
  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr
  
  cna = cna %>%
    mutate(begin = begin + vfrom[chr],
           end = end + vfrom[chr])
  return(cna)
}


base = '/orfeo/cephfs/scratch/cdslab/shared/SCOUT/validation_subclonal_new/tables/'
sp = 'SPN06'
cov = 100
pur = 0.9
spns = paste0('SPN0', c(4,6))

all_samples <- c()
for (sp in spns){
  samples <- get_sample_names(sp)
  all_samples <- c(all_samples, paste0(sp, '_', samples))
}

plots_sample <- list()
plots_sign <-  list()
plots_bar <- list()

for (sp in spns){
  print(sp)
  gt =  readRDS(paste0(base, 'table_process_univariate_w_private_', sp, '_', cov, 'x_', pur, 'p_mutect2_ascat.rds')) %>% 
    filter(!str_detect(causes, "errors")) 
  
  samples = get_sample_names(sp)
  df_spn <- tibble()
  s=samples[[1]]
  for (s in samples){
    print(s)
    tmp_name = paste0(sp, '_', s)
    
    if (s == 'SPN06_1.1'){
      gt_f = gt %>% 
        filter(sample_id == paste0(sp, '_', s)) %>% 
        group_by(cluster_id_process, sample_id) %>%
        mutate(is_clonal_process=replace(FALSE, ccf_process > 0.97, TRUE)) %>% 
        ungroup() %>%
        mutate(cluster_id_process_full = cluster_id_process) %>%
        mutate(cluster_id_process = replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal')) %>% 
        separate(cluster_id_process, sep = '_', into = c('cluster_id_process', 'tmp')) %>% 
        select(-tmp) 
      
    } else {
      gt_f = gt %>% 
        filter(sample_id == paste0(sp, '_', s)) %>% 
        group_by(cluster_id_process, sample_id) %>%
        mutate(is_clonal_process=replace(FALSE, ccf_process > 0.95, TRUE)) %>% 
        ungroup() %>%
        mutate(cluster_id_process_full = cluster_id_process) %>%
        mutate(cluster_id_process = replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal')) %>% 
        separate(cluster_id_process, sep = '_', into = c('cluster_id_process', 'tmp')) %>% 
        select(-tmp) 
      
    }
    
    gt_f <- gt_f %>% 
      mutate(cluster_id_process = ifelse(cluster_id_process == 'Clonal', 'Clonal', 'Other'))
    
    filt_gt <- gt_f %>% 
      filter(!grepl(":(X|Y):", mutation_id))
    
    if (nrow(filt_gt)>0){
      
      if (sp == 'SPN04'){
        ymax = 150
      } else {
        ymax = 1200
      }
      
      plots_sample[[s]]  <- filt_gt %>% 
        ggplot()+
        geom_histogram(aes(x = vaf_process, fill = cluster_id_process), alpha=0.6, binwidth = 0.01, position = "identity") +
        xlim(0.05,1.01) +
        ylim(0, ymax) + 
        xlab(paste0('VAF ', s)) + 
        my_ggplot_theme() + 
        scale_fill_manual('ProCESS clusters', values = color_palette_process_uni) 
      
      
      if (sp == 'SPN04'){
        filt_gt = filt_gt %>% 
          filter(causes %in% names(sbs_colors)) %>% 
          mutate(causes = ifelse(causes =='SBS25', causes, 'Other'))
      
        plots_sign[[s]] <- filt_gt %>% 
          ggplot()+
          geom_histogram(aes(x = vaf_process, fill = causes), alpha=0.6, binwidth = 0.01, position = "identity") +
          xlim(0.05,1.01) +
          ylim(0, ymax) + 
          xlab(paste0('VAF ', s)) + 
          my_ggplot_theme() + 
          scale_fill_manual('SBS', values = c(color_sbs, 'Other'='gray70')) 
        
      } else {
        filt_gt = filt_gt %>% 
          filter(causes %in% names(sbs_colors)) %>% 
          mutate(causes = ifelse(causes =='SBS31', causes, 'Other'))
        
        plots_sign[[s]] <- filt_gt %>% 
          ggplot()+
          geom_histogram(aes(x = vaf_process, fill = causes), alpha=0.6, binwidth = 0.01, position = "identity") +
          xlim(0.05,1.01) +
          ylim(0, ymax) + 
          xlab(paste0('VAF ', s)) + 
          my_ggplot_theme() + 
          scale_fill_manual('SBS', values = c(color_sbs, 'Other'='gray70'))
      }
    }
    df_spn <- bind_rows(df_spn, filt_gt)
    
    if (sp == 'SPN06'){
      plots_bar[[sp]] <-  df_spn %>% 
        filter(sample_id %in% c("SPN06_SPN06_1.2", "SPN06_SPN06_2.1")) %>% 
        group_by(sample_id, cluster_id_process) %>% 
        summarise(n =n()) %>% 
        filter(cluster_id_process == 'Clonal') %>% 
        #mutate(sample_id = sub("^([^_]+)_\\1", "\\1", sample_id)) %>% 
        mutate(sample_id = case_when(
          grepl("SPN06_SPN06_1.2", sample_id) ~ "Pre",
          grepl("SPN06_SPN06_2.1", sample_id) ~ "Post"
        )) %>% 
        mutate(sample_id = factor(sample_id, levels = c('Pre', 'Post'))) %>% 
        ggplot()+
        geom_col(aes(x = sample_id, y = n, fill = cluster_id_process), position = "dodge")+ 
        my_ggplot_theme() + 
        scale_fill_manual('ProCESS clusters', values = color_palette_process_uni)  +
        xlab('Sample') +
        ylab('Number of mutations')
    } else {
      plots_bar[[sp]] <-  df_spn %>% 
        group_by(sample_id, cluster_id_process) %>% 
        summarise(n =n()) %>% 
        filter(cluster_id_process == 'Clonal') %>% 
        mutate(sample_id = sub("^([^_]+)_\\1", "\\1", sample_id)) %>% 
        mutate(sample_id = case_when(
          grepl("SPN04_1.1", sample_id) ~ "Pre",
          grepl("SPN04_2.1", sample_id) ~ "Post"
        )) %>% 
        mutate(sample_id = factor(sample_id, levels = c('Pre', 'Post'))) %>% 
        ggplot()+
        geom_col(aes(x = sample_id, y = n, fill = cluster_id_process), position = "dodge") + 
        my_ggplot_theme() + 
        scale_fill_manual('ProCESS clusters', values = color_palette_process_uni) +
        xlab('Sample') +
        ylab('Number of mutations')
    }
  }
}
color_sbs <- c('SBS25' = 'sienna2', 'SBS31' = 'maroon' )
sbs_colors[['SBS25']] = 'sienna2'
sbs_colors[['SBS31']] = 'maroon'

library(ProCESS)
forest = load_sample_forest(get_sample_forest(spn = 'SPN04'))
f_spn04 <- plot_forest(forest, color_map = c("Clone 1"='gainsboro',
                                             "Clone 2"='gainsboro',
                                             "Clone 3"='gainsboro',
                                             "Clone 4"='gainsboro',
                                             "Clone 5"='gainsboro'))  %>%
  annotate_forest(forest)

forest = load_sample_forest(get_sample_forest(spn = 'SPN06'))
sub_06 <- forest$get_subforest_for(c("SPN06_1.2", "SPN06_2.1"))
f_spn06 <- plot_forest(sub_06, color_map = c("Clone 1"='gainsboro',
                                             "Clone 2"='gainsboro',
                                             "Clone 3"='gainsboro',
                                             "Clone 4"='gainsboro',
                                             "Clone 5"='gainsboro'))  %>%
  annotate_forest(sub_06)


source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_5/signatures_edf.R')
sig_04 <- df_simultated_signatures %>% 
  filter(sample %in% c('SPN04_1.1', 'SPN04_2.1')) %>% 
  filter(signature %in% names(sbs_colors)) %>% 
  filter(coverage == 100, purity == 0.9) %>% 
  mutate(sample = case_when(
    grepl("SPN04_1.1", sample) ~ "Pre",
    grepl("SPN04_2.1", sample) ~ "Post"
  )) %>% 
  ggplot() + 
  geom_col(aes(y = sample, x = exposure, fill =signature )) +
  my_ggplot_theme() + 
  ylab('Exposure') +
  xlab('Sample') +
  scale_fill_manual('SBS', values = c(sbs_colors, 'Other'='gray70'))
  

sig_06 <- df_simultated_signatures %>% 
  filter(sample %in% c('SPN06_1.2', 'SPN06_2.1')) %>% 
  filter(signature %in% names(sbs_colors)) %>% 
  filter(coverage == 100, purity == 0.9) %>% 
  mutate(sample = case_when(
    grepl("SPN06_1.2", sample) ~ "Pre",
    grepl("SPN06_2.1", sample) ~ "Post"
  )) %>% 
  ggplot() + 
  geom_col(aes(y = sample, x = exposure, fill =signature )) +
  my_ggplot_theme() + 
  ylab('Exposure') +
  xlab('Sample') +
  scale_fill_manual('SBS', values = c(sbs_colors, 'Other'='gray70'))

spn04 <- free(f_spn04 + theme(legend.position = 'none')) + 
  plots_sample$SPN04_1.1 + 
  #ggtitle('Pre-Treatment') + 
  plots_sample$SPN04_2.1 + 
  #ggtitle('Post-Treatment') +
  plots_bar$SPN04 + 
  plots_sign$SPN04_2.1 + 
  #ggtitle('Post-Treatment') +
  sig_04 + 
  plot_layout(guides = 'collect', nrow = 1, design = 'AAA#BBDEE\nAAA#CCDFF') 

spn06 <- free(f_spn06 +  theme(legend.position = 'none')) + 
  plots_sample$SPN06_1.2 + 
  #ggtitle('Pre-Treatment') + 
  plots_sample$SPN06_2.1 + 
  #ggtitle('Post-Treatment') +
  plots_bar$SPN06 + 
  plots_sign$SPN06_2.1 + 
  #ggtitle('Post-Treatment') +
  sig_06 + 
  plot_layout(guides = 'collect', nrow = 1, design = 'AAA#BBDEE\nAAA#CCDFF') 


all <- wrap_plots(spn04, spn06, nrow = 2, guides = 'collect')
ggsave(filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/mrca.pdf',
       plot = all, 
       width = 8.5, 
       height = 4, 
       units = 'in' )
