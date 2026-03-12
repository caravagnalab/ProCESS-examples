library(ggplot2)
library(tidyverse)
library(ProCESS)

source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/SCOUT/colors.R')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils.R")


df_process <- readRDS('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/process.rds')
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

results <- lapply(1:nrow(combs), FUN = function(i){
  sp = combs[i, "spn"]
  cov = combs[i, "coverage"]
  pur = combs[i, "purity"]
  s_tool = combs[i, "sig_tool"]
  tool = combs[i, "tool"]
  
  data_tool = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/res_int/', cov, 'x_', pur, 'p_', s_tool,'_',tool)
  df_tool = readRDS(file = paste0(data_tool, '_df.rds')) %>% 
    filter(spn == sp) %>% 
    select(cluster, contains_driver_process, score_all, score_driver, score_tail, score_sign, driver_label_process)
  
  tmp_process = df_process %>% 
    filter(spn == sp, coverage == cov, purity == pur) %>% 
    select(cluster, score_all, score_driver, score_tail, score_sign, driver_label_process) %>% 
    filter(cluster != 'Subclonal')
  
  join = left_join(tmp_process, df_tool, 
                   by = join_by(driver_label_process), 
                   suffix = c('_process', '_tool')) %>% 
    mutate(spn = sp, coverage = cov, purity = pur, sig_tool = s_tool, sub_tool = tool)
  
}) %>% 
  bind_rows() %>%  
  mutate(class_process = case_when(
    score_all_process >=.70 &  cluster_process!= 'Clonal' ~ 'Strong\nExpansion',
    score_all_process <.70 & score_all_process >=.2  &  cluster_process!= 'Clonal' ~ 'Medium\nExpansion',
    score_all_process <.2   &  cluster_process!= 'Clonal' ~ 'Neutral\nExpansion',
    cluster_process == 'Clonal' ~ 'Clonal'
  )) %>% 
  mutate(class_process = factor(class_process, levels = c('Clonal', 'Strong\nExpansion', 'Medium\nExpansion', 'Neutral\nExpansion'))) %>% 
  mutate(class_tool = case_when(
    score_all_tool >=.70 ~ 'Strong\nExpansion',
    score_all_tool <.70 & score_all_tool >=.2 ~ 'Medium\nExpansion',
    score_all_tool <.2   ~ 'Neutral\nExpansion'
  )) %>% 
  mutate(class_tool = factor(class_tool, levels = c('Clonal', 'Strong\nExpansion', 'Medium\nExpansion', 'Neutral\nExpansion'))) 
  

scatter <- results %>% 
  ggplot() +
 
  annotate("rect", xmin = 0.7, xmax = 1, ymin = 0.7, ymax = 1, 
           fill = "palegreen4", alpha = 0.2) + 
  annotate("rect", xmin = 0.2, xmax = 0.7, ymin = 0.2, ymax = 0.7, 
           fill = "goldenrod", alpha = 0.2) + 
  annotate("rect", xmin = 0.2, xmax = 0, ymin = 0.2, ymax = 0, 
           fill = "gainsboro", alpha = 0.2) +
  geom_point(aes(x = score_all_process, y = score_all_tool,  col = as.factor(sub_tool)), width = .05) +
  geom_smooth(aes(x = score_all_process, y = score_all_tool,  col = as.factor(sub_tool)), method = 'lm') +
  my_ggplot_theme() +
  xlab('Score ProCESS') +
  ylab('Score Tool')



plot_tool <- results %>% 
  mutate(score_all_diff = score_all_process - score_all_tool) %>% 
  ggplot() + 
  geom_boxplot(aes(x = class, y = score_all_diff, col = sub_tool, fill = sub_tool), alpha = .3) +
  my_ggplot_theme()  +
  ylab("Score All\nProCESS - Tool") +
  scale_color_manual('Tool', values = c('mediumorchid4', 'darkgoldenrod3')) + 
  scale_fill_manual('Tool', values = c('mediumorchid4', 'darkgoldenrod3')) +
  facet_grid(.~purity)


plot_score <- results %>% 
  mutate(score_all_diff = score_all_process - score_all_tool,
         score_driver_diff = score_driver_process - score_driver_tool,
         score_tail_diff = score_tail_process - score_tail_tool,
         score_sign_diff = score_sign_process - score_sign_tool) %>% 
  pivot_longer(cols = c(score_all_diff, score_driver_diff, score_tail_diff, score_sign_diff)) %>% 
  mutate(name = case_when(
    name == 'score_all_diff' ~ 'All',
    name == 'score_tail_diff' ~ 'Tail',
    name == 'score_sign_diff' ~ 'Signature',
    name == 'score_driver_diff' ~ 'Driver'
  )) %>% 
  ggplot() + 
  geom_boxplot(aes(x = class, y = value, col = name, fill = name), alpha = .3) +
  my_ggplot_theme()  +
  ylab("Score Error\nProCESS - Tool") +
  scale_color_manual('Score', values = c('plum4', 'darkseagreen4', 'cadetblue4', 'salmon2')) + 
  scale_fill_manual('Score', values = c('plum4', 'darkseagreen4', 'cadetblue4', 'salmon2'))



coverage_list = c(50, 100, 150)
purity_list = c(0.9, 0.6, 0.3)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04','SPN05', 'SPN06', 'SPN07')
tool_list = c( 'viber')
sig_tool_list = c('BASCULE', 'SigProfiler')
combs = expand.grid(spn = spn_list, 
                    coverage=coverage_list,
                    purity=purity_list,
                    tool=tool_list,
                    sig_tool = sig_tool_list)

i=1
join <- lapply(1:nrow(combs), FUN = function(i){
  sp = combs[i, "spn"]
  cov = combs[i, "coverage"]
  pur = combs[i, "purity"]
  s_tool = combs[i, "sig_tool"]
  tool = combs[i, "tool"]
  
  df_process <- readRDS('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/process.rds') %>% 
    filter(sp == spn, 
           coverage == cov, 
           purity == pur) %>% 
    select(cluster, score_all) %>% 
    mutate(class_process = case_when(
      score_all >=.70 ~ 'Strong\nExpansion',
      score_all <.70 & score_all >=.2  &  cluster!= 'Clonal' ~ 'Medium\nExpansion',
      score_all <.2   &  cluster!= 'Clonal' ~ 'Neutral\nExpansion',
    ))  %>% 
    dplyr::rename(cluster_id_process = cluster, score_all_process = score_all)
  
  data_tool = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/res_int/', cov, 'x_', pur, 'p_', s_tool,'_',tool)
  df_tool = readRDS(file = paste0(data_tool, '_df.rds')) %>% 
    filter(spn == sp) %>% 
    select(cluster, score_all) %>% 
    mutate(class_tool = case_when(
      score_all >=.70 ~ 'Strong\nExpansion',
      score_all <.70 & score_all >=.2   ~ 'Medium\nExpansion',
      score_all <.2  ~ 'Neutral\nExpansion',
    )) %>% 
    dplyr::rename(cluster_id_tool = cluster, score_all_tool = score_all)
    
  if (file.exists(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal/tables_interpreted/', tool, '_', sp, '_', cov, 'x_', pur, 'p_mutect2_ascat.rds'))){
    mutations <-  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal/tables_interpreted/', tool, '_', sp, '_', cov, 'x_', pur, 'p_mutect2_ascat.rds'))
    join <- mutations %>% 
      select(mutation_id, cluster_id_tool, cluster_id_process) %>% 
      left_join(df_tool) %>% 
      left_join(df_process) %>% 
      distinct()
    return(join)
  }
  })
saveRDS(file = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/interpretation/performance.rds', object = join)

join %>% 
  bind_rows() %>% 
  mutate(class_process = factor(class_process, levels = c('Clonal', 'Strong\nExpansion', 'Medium\nExpansion', 'Neutral\nExpansion')))  %>% 
  ggplot(aes(x = class_process, y = score_all_process - score_all_tool, col = class_process, fill = class_process)) +
  geom_boxplot(alpha = .3) + 
  theme_minimal() +
  xlab('ProCESS Class') +
  ylab("Score All\nProCESS - Tool") +
  scale_color_manual('ProCESS class', values = c('palegreen4', 'goldenrod', 'gainsboro')) + 
  scale_fill_manual('ProCESS class', values = c('palegreen4', 'goldenrod', 'gainsboro')) 


