setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
source('../../getters/tumourevo_getters.R')
source('../../getters/process_getters.R')
source('/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/utils_plot.R')

color_palette_process_uni = RColorBrewer::brewer.pal(n = 8, name = "Set1") 
names(color_palette_process_uni) <- c('Clonal', paste0('Clone ', 1:7))
color_palette_process_uni['Subclonal'] = 'gray70'

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
spns = paste0('SPN0', 1:7)

all_samples <- c()
for (sp in spns){
  samples <- get_sample_names(sp)
  all_samples <- c(all_samples, paste0(sp, '_', samples))
}

class <- c(NA, 'Monoclonal', NA, 
           'Monoclonal', 'Monoclonal', 
           'Monoclonal', 'Polyclonal','Monoclonal', 'Polyclonal',
           'Monoclonal', 'Monoclonal', 
           'Polyclonal','Monoclonal','Polyclonal', 
           'Monoclonal', 'Polyclonal', 'Monoclonal', 'Monoclonal', NA,
           'Monoclonal', 'Polyclonal', 'Monoclonal', 'Monoclonal', 'Monoclonal')
sub_class <- c(NA, 'Monoclonal', NA, 
               'Monoclonal', 'Sampling Bias', 
               'Monoclonal', 'Weak evidence','Sampling Bias', 'Polyclonal',
               'Monoclonal', 'Sampling Bias', 
               'Weak evidence','Monoclonal','Weak evidence', 
               'Monoclonal', 'Polyclonal', 'Monoclonal', 'Monoclonal', NA,
               'Sampling Bias', 'Polyclonal', 'Monoclonal', 'Monoclonal', 'Sampling Bias')
table = tibble(sample_id = all_samples, class =class, subclass =  sub_class)
saveRDS(object = table, 'class_table.rds')

plots_sample <- list()

for (sp in spns){
  print(sp)
  gt =  readRDS(paste0(base, 'table_process_univariate_w_private_', sp, '_', cov, 'x_', pur, 'p_mutect2_ascat.rds')) %>% 
    filter(!str_detect(causes, "errors")) 
  
  samples = get_sample_names(sp)
  
  s=samples[[1]]
  for (s in samples){
    print(s)
    tmp_name = paste0(sp, '_', s)
    # cn = readRDS(get_process_cna(spn = as.character(sp), sample = s)) %>%
    #   as.data.frame() %>%                    # <-- add this
    #   group_by(chr, begin, end) %>%
    #   group_modify(~ {
    #     if (any(.x$major != 1 | .x$minor != 1)) {
    #       .x %>%
    #         filter(major != 1 | minor != 1) %>%
    #         slice_max(ratio, n = 1, with_ties = FALSE)
    #     } else {
    #       .x %>%
    #         slice(1)
    #     }
    #   }) %>%
    #   ungroup() %>%
    #   select(-ratio)
    # 
    cn = readRDS(get_process_cna(spn = as.character(sp), sample = s)) %>%
      as.data.frame() %>%
      group_by(chr, begin, end) %>%
      filter(if (any(major != 1 | minor != 1)) (major != 1 | minor != 1) else TRUE) %>%
      slice_max(ratio, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(-ratio)
    
    gt_cn = gt %>%
      filter(sample_id == paste0(sp, '_', s)) %>% 
      select(mutation_id) %>%
      separate(mutation_id, into = c('spn', 'chr', 'pos', 'alt'), remove = F, convert = T) %>%
      left_join(cn, by = join_by(chr)) %>%
      filter(pos >= begin & pos <= end) %>%
      mutate(CN = paste(major, minor, sep = ':')) %>%
      select(mutation_id, CN) %>% 
      distinct()
    
    if (s == 'SPN06_1.1'){
      gt_f = gt %>% 
        filter(sample_id == paste0(sp, '_', s)) %>% 
        left_join(gt_cn) %>% 
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
        left_join(gt_cn) %>% 
        group_by(cluster_id_process, sample_id) %>%
        mutate(is_clonal_process=replace(FALSE, ccf_process > 0.95, TRUE)) %>% 
        ungroup() %>%
        mutate(cluster_id_process_full = cluster_id_process) %>%
        mutate(cluster_id_process = replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal')) %>% 
        separate(cluster_id_process, sep = '_', into = c('cluster_id_process', 'tmp')) %>% 
        select(-tmp) 

    }
    
    
    filt_gt <- gt_f %>% 
      filter(CN == '1:1' | is.na(CN)) 

    if (nrow(filt_gt)>0){
      tmp_class = table %>% filter(sample_id == tmp_name) %>% pull(subclass)
      
      plots_sample[[s]]  <- filt_gt %>% 
        ggplot()+
        geom_histogram(aes(x = vaf_process, fill = cluster_id_process), alpha=0.6, binwidth = 0.01, position = "identity") +
        xlim(0.05,1.01) +
        xlab(paste0('VAF ', s)) + 
        my_ggplot_theme() + 
        ggtitle(tmp_class) + 
        scale_fill_manual('ProCESS clusters', values = color_palette_process_uni) 
    }
  }
}

final = wrap_plots(plots_sample, ncol = 3)
ggsave(filename = '/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/plot_mobster_diploid.png', 
       plot = final, width = 12, height = 15, dpi = 200, units = 'in')
