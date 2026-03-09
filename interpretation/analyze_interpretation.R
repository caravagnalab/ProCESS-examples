library(ggplot2)
library(tidyverse)
library(ProCESS)

source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/SCOUT/colors.R')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils.R")


colors_cluster = c('indianred', 
                   'steelblue', 
                   'forestgreen', 
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

vcf_caller = "mutect2"
cna_caller = "ascat"

coverage_list = c(50,100,150)
purity_list = c(0.3, 0.6, 0.9)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04','SPN05', 'SPN06', 'SPN07')
tool_list = c( 'viber', 'pyclonevi')
sig_tool_list = c('BASCULE', 'SigProfiler')

base <- '/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation'
subclonal <- readRDS(paste0(base, '/subclonal_deconvolution_new.rds')) %>% 
  dplyr::rename(spn = patient_id, 
                dec_tool = tool, 
                cluster = cluster_id_tool)
signature <- readRDS(paste0(base, '/signature_deconvolution_summary.rds'))
signature$cluster <- sub("^X", "", signature$cluster)

all <- full_join(subclonal, signature)

combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    tool=tool_list,
                    sig_tool = sig_tool_list)

i=26

for (i in 1:nrow(combs)){
  cov = combs[i, "coverage"]
  pur = combs[i, "purity"]
  s_tool = combs[i, "sig_tool"]
  tool = combs[i, "tool"]
  
  f_all <- all %>% 
    filter(coverage == cov,
           purity == pur, 
           dec_tool == tool,
           sig_tool == s_tool | is.na(sig_tool), 
           spn %in% spn_list)
  
  f_all <- f_all %>% 
   group_by(across(c(-sig,-cs_nsig, -cs_exp, -match))) %>% 
   summarize(match = sum(match),
          cs_exp = mean(cs_exp)) %>% 
   mutate(match = ifelse(match >1, T, F)) %>% 
   ungroup()
  
  
  w_driver = 0.2
  w_match_sig = 0.2
  w_tail_sample = 0.2
  w_perc_tail = 0.2
  w_cs_sign = 0.2
  
  df_score <- f_all %>% 
    mutate(class_driver = ifelse(contains_driver_tool == F, 0, 1),
           class_match_sig = ifelse(match == F, 1, 0),
           class_match_sig = ifelse(is.na(match) & is_clonal_tool == F, 0, class_match_sig), 
           class_match_sig = ifelse(is.na(match) & is_clonal_tool == T, 1, class_match_sig), 
           class_tail_sample =  ifelse(one_sample_tail == F, 1, 0), #almeno un sample
           class_perc_tail = (100-percent_tail)/100,
           class_cs_sign = 1-cs_exp,
           class_cs_sign = ifelse(is.na(cs_exp) & is_clonal_tool == F, 0, class_cs_sign),
           class_cs_sign = ifelse(is.na(cs_exp) & is_clonal_tool == T, 1, class_cs_sign)) %>% 
    mutate(score_driver = 1*class_driver) %>% 
    mutate(score_all = (class_driver + 
             n_never_tail + 
             class_cs_sign)/3) %>% 
    mutate(score_sign = 0.5*class_driver + 
             0.5*class_cs_sign ) %>% 
    mutate(score_tail= 0.5*class_driver + 
             0.5*n_never_tail ) %>% 
    mutate(score_all = round(score_all, 2),
           score_tail = round(score_tail, 2),
           score_sign = round(score_sign, 2)) %>% 
    filter(total_mutations >= 100) 
  
  # df_score <- df_score %>% 
  #   mutate(class = case_when(
  #     score_all >= .95 ~ '5',
  #     score_all < .95 & score_all >= 0.8 ~ '4',
  #     score_all < .8 & score_all >= 0.6 ~ '3',
  #     score_all < .6 & score_all > .1 ~ '2',
  #     score_all < .1 ~ '1'
  #   )) %>% 
  #   mutate(status = case_when(
  #     class %in% c('5', '4', '3') & contains_driver_process == T ~ 'TP',
  #     class %in% c('5', '4', '3') & contains_driver_process == F ~ 'FP',
  #     class %in% c('2', '1') & contains_driver_process == F ~ 'TP', #TN
  #     class %in% c('2', '1') & contains_driver_process == T ~ 'FN'
  #   ))
  
  # res <- df_score %>% 
  #   select(spn, cluster, score_all, contains_driver_process, class, status) 
  
  # summary <- res %>% 
  #   group_by(spn) %>% 
  #   summarize(TP = sum(status == 'TP'),
  #             TN = sum(status == 'TN'),
  #             FN = sum(status == 'FN'),
  #             FP = sum(status == 'FP'),
  #             n = n()) %>% 
  #   mutate(accuracy =  (TP + TN) / (TP + TN + FP + FN),
  #          precision = TP / (TP + FP) ,
  #          recall = TP / (TP + FN), #Sensitivity / TPR
  #          specificity = TN / (TN + FP), #TNR
  #          f1 = 2 * (precision * recall) / (precision + recall)) %>% 
  #   mutate(coverage = cov, purity = pur, tool_sub = tool, tool_sig = s_tool)
  
  if (tool == 'pyclonevi'){
    df_score <- df_score %>% 
      mutate(cluster = paste0('C', cluster))
  }
  
  plt <- df_score %>% 
    pivot_longer(cols = c(score_driver, score_all, score_tail, score_sign)) %>% 
    mutate(name = factor(name, levels = c('score_driver', 'score_tail', 'score_sign', 'score_all'))) %>% 
    ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.95, ymax = 1, 
             fill = "palegreen3", alpha = 0.2) + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.8, ymax = 0.95, 
             fill = "goldenrod", alpha = 0.2) + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.6, ymax = 0.8, 
             fill = "chocolate3", alpha = 0.2) + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = 0.6, 
             fill = "indianred4", alpha = 0.2) + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = 0, 
             fill = "gainsboro", alpha = 0.2) + 
    geom_point(aes(x = name, y = value, col = cluster, shape = contains_driver_process), size = 4) +
    geom_line(aes(x = name, y = value, col = cluster, group = cluster), linewidth = .6)  +
    theme_minimal() +
    scale_shape_manual('Contains True Driver', values = c(4, 20)) + 
    facet_wrap(.~spn) +
    ylab('Score') +
    xlab('')+ 
    scale_color_manual('Cluster', 
                       values = colors_cluster) +
    scale_x_discrete(labels = c('score_driver' = 'Only\nDriver', 
                                'score_tail'   = 'Driver\nTail', 
                                'score_sign'   = 'Driver\nSignature', 
                                'score_all'    = 'All')) +
    ggtitle(paste0(cov, 'x_',pur, 'p_', tool, '_', s_tool)) + my_ggplot_theme()
  
  name_file = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/plot_int/', cov, 'x_', pur, 'p_', s_tool,'_',tool, '.png')
  ggsave(plot = plt, 
         filename = name_file, 
         dpi = 400, 
         width = 10, height = 10)
  
  path = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/res_int/', cov, 'x_', pur, 'p_', s_tool,'_',tool)
  saveRDS(object = summary, file = paste0(path, '_summary.rds'))
  saveRDS(object = df_score, file = paste0(path, '_df.rds'))
}


all_res <- lapply(1:nrow(combs), FUN = function(i){
  cov = combs[i, "coverage"]
  pur = combs[i, "purity"]
  s_tool = combs[i, "sig_tool"]
  tool = combs[i, "tool"]

  path = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/res_int/', cov, 'x_', pur, 'p_', s_tool,'_',tool)
  return(readRDS(file =paste0(path, '_summary.rds')))

  }) %>% bind_rows()




all_res %>% 
  pivot_longer(cols = c(accuracy, precision, recall, f1)) %>% 
  ggplot() +
  geom_boxplot(aes(x = as.factor(purity), y = value, fill = tool_sub)) + 
  my_ggplot_theme() +
  facet_grid(.~name)
