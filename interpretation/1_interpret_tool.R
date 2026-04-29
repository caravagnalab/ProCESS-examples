setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')
source("../validation/SCOUT/colors.R")
source("../figures/figure3/utils_plot.R")
source("../figures/figure3/utils.R")


df_score <- tibble(spn = paste0('SPN0', 1:7),
                   signature = c(F,F,F,F,F,T,T),
                   driver = c(T,T,T,T,T,T,T),
                   tail = c(T,T,T,T,T,T,T))


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

option_list <- list(make_option(c("--cna_caller"), type = "character", default = 'sequenza'),
                    make_option(c("--vcf_caller"), type = "character", default = 'strelka')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

vcf_caller = opt$vcf_caller #"mutect2"
cna_caller = opt$cna_caller #"sequenza"

print(vcf_caller)
print(cna_caller)

base = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/interpretation_', vcf_caller, "_", cna_caller, '/')
dir.create(paste0(base, '/plot_int'), recursive = T, showWarnings = F)
dir.create(paste0(base, '/res_int'), recursive = T, showWarnings = F)

coverage_list = c(50, 100, 150)
purity_list = c(0.9, 0.6, 0.3)
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04','SPN05', 'SPN06', 'SPN07')
tool_list = c( 'viber', 'pyclonevi')
sig_tool_list = c('BASCULE', 'SigProfiler')

subclonal <- readRDS(paste0(base, '/subclonal_deconvolution.rds')) %>% 
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

i=2

for (i in 1:nrow(combs)){
  print(i)
  cov = combs[i, "coverage"]
  pur = combs[i, "purity"]
  s_tool = combs[i, "sig_tool"]
  tool = combs[i, "tool"]
  
  tmp_df <-  all %>% 
    filter(coverage == cov,
           purity == pur, 
           dec_tool == tool,
           spn %in% spn_list) %>% 
    group_by(spn, sig_tool) %>% 
    summarise() %>% 
    filter(!is.na(sig_tool)) %>% 
    filter(sig_tool == s_tool)
  
  f_all <- all %>%
    filter(
      coverage == cov,
      purity == pur,
      dec_tool == tool,
      spn %in% spn_list,
      sig_tool == s_tool | (is.na(sig_tool) & spn %in% tmp_df$spn)
    )
  
  f_all <- f_all %>% 
   group_by(across(c(-sig, -cs_exp, -match, -n_rel, -bg_match, -bg_cs_exp, -bg_n_rel))) %>% 
   summarize(cs_exp = min(cs_exp),
          n_rel = max(n_rel),
          bg_cs_exp = min(bg_cs_exp)) %>% 
   ungroup()
  
  f_all <- f_all %>% left_join(df_score)
  
  score <- df_score %>% 
    mutate(w_signature = as.numeric(signature),
           w_driver = as.numeric(driver),
           w_tail = as.numeric(tail)) %>% 
    rowwise() %>% 
    mutate(n = sum(w_signature, w_driver, w_tail)) %>% 
    select(spn, w_signature, w_driver, w_tail, n)
  
  f_all <- f_all %>% left_join(score)
    
  interpret_table <- f_all %>% 
    mutate(driver = ifelse(contains_driver_tool == F, 0, 1),
           bg_cs = 1 - bg_cs_exp,
           cs_sign = 1 - cs_exp,
           cs_sign = ifelse(is.na(cs_exp) & is_clonal_tool == T, 1, cs_sign),
           bg_sign = ifelse(is.na(bg_cs) & is_clonal_tool == T, 1, bg_cs),
           n_rel = ifelse(is.na(cs_exp) & is_clonal_tool == T, 1, n_rel)) %>% 
    rowwise() %>% 
    mutate(score_spn = (w_driver * driver + w_tail * n_never_tail + w_signature* ((cs_sign+n_rel+bg_sign)/3))/n,
           score_all = (driver + n_never_tail + ((cs_sign+n_rel+bg_sign)/3))/3,
           score_driver = driver, 
           score_sign = (cs_sign+n_rel+bg_sign)/3,
           score_tail = n_never_tail,
           score_no_tail = (driver + ((cs_sign+n_rel+bg_sign)/3))/2,
           score_no_driver = (((cs_sign+n_rel+bg_sign)/3) + n_never_tail)/2,
           score_no_sign = (driver + n_never_tail)/2)

  
  if (tool == 'pyclonevi'){
    interpret_table <- interpret_table %>% 
      mutate(cluster = paste0('C', cluster))
  }
  
  plt <- interpret_table %>%
    filter(total_mutations > 200) %>% 
    pivot_longer(cols = c(score_driver, score_all, score_tail, score_no_driver, score_no_tail, score_no_sign, score_sign)) %>%
    mutate(name = factor(name, levels = c('score_driver', 'score_tail', 'score_sign','score_no_driver', 'score_no_tail', 'score_no_sign', 'score_all'))) %>%
    ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.9, ymax = 1, 
             fill = "palegreen4", alpha = 0.2) + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.55, ymax = .9, 
             fill = "goldenrod", alpha = 0.2) + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.2, ymax = 0.55, 
             fill = "salmon1", alpha = 0.2) + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.2, ymax = 0, 
             fill = "gainsboro", alpha = 0.2) +
    geom_point(aes(x = name, y = value, col = cluster, shape = contains_driver_process), size = 4) +
    geom_line(data = ~ filter(.x, !is.na(value)), aes(x = name, y = value, col = cluster, group = cluster), linewidth = .6)  +
    geom_text(data = ~ filter(.x, is_clonal_tool & name == 'score_all'), aes(x = name, y = value+0.07, col = cluster, group = cluster, label = 'Clonal')) + 
    geom_text(data = ~ filter(.x, contains_driver_process & name == 'score_all'), aes(x = name, y = value+0.03, col = cluster, group = cluster, label = driver_label_process)) + 
    theme_minimal() +
    scale_shape_manual('Contains True Driver', values = c(4, 20)) +
    facet_wrap(.~spn) +
    ylab('Score') +
    xlab('')+
    scale_color_manual('Cluster',
                       values = colors_cluster) +
    scale_x_discrete(labels = c('score_driver' = 'Driver',
                                'score_tail'   = 'Tail',
                                'score_sign'   = 'Signature',
                                'score_no_driver' = 'Signature\nTail',
                                'score_no_tail' = 'Driver\nSignature',
                                'score_no_sign' = 'Driver\nTail',
                                'score_all'    = 'All')) +
    ggtitle(paste0(cov, 'x_',pur, 'p_', tool, '_', s_tool)) + 
    my_ggplot_theme() 
  
  
  name_file = paste0(base,'/plot_int/', cov, 'x_', pur, 'p_', s_tool,'_',tool, '.png')
  ggsave(plot = plt,
         filename = name_file,
         dpi = 400,
         width = 10, height = 10)

  path = paste0(base, '/res_int/', cov, 'x_', pur, 'p_', s_tool,'_',tool)
  saveRDS(object = interpret_table %>% filter(total_mutations > 200), file = paste0(path, '_df.rds'))
}
