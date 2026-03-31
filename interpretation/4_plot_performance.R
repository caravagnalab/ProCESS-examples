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

source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils.R")

out = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/"
base = '/orfeo/cephfs/scratch/cdslab/shared/SCOUT/validation_subclonal/tables/'

color_palette_process = RColorBrewer::brewer.pal(n = 8, name = "Dark2") 
names(color_palette_process) <- c('Clonal', paste0('Clone ', 1:7))
color_palette_process['Subclonal'] = 'gray70'

color_class <- c("palegreen4","darkseagreen", "goldenrod", "indianred", "gainsboro")
names(color_class) <- c('Clonal', 'Strong\nExpansion', 'Medium\nExpansion', 'Low\nExpansion', 'Neutral\nExpansion')

color_score <- c('plum4', 'darkseagreen4', 'cadetblue4', 'salmon2')
score_types <- c("all", 'no_tail', 'no_driver', 'no_sign')
names(color_score) <- score_types

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

sp = 'SPN01'
cov = 50
pur = 0.9
tool = 'pyclonevi'
s_tool = 'SigProfiler'

coverage_list = c(50, 100, 150)
purity_list = c(0.9, 0.6, 0.3)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04','SPN05', 'SPN06', 'SPN07')
tool_list = c( 'viber', 'pyclonevi')
sig_tool_list = c('BASCULE', 'SigProfiler')

# combs = expand.grid(spn = spn_list, 
#                     coverage=coverage_list,
#                     purity=purity_list,
#                     tool=tool_list,
#                     sig_tool = sig_tool_list)
# results <- lapply(1:nrow(combs), FUN = function(i){
#   sp = combs[i, "spn"]
#   cov = combs[i, "coverage"]
#   pur = combs[i, "purity"]
#   s_tool = combs[i, "sig_tool"]
#   tool = combs[i, "tool"]

full_class_mutations <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/interpretation/performance_mutations_class.rds')
for (sp in spn_list){

  df_process <- readRDS(paste0(out, 'process.rds')) %>% 
    filter(spn == sp, coverage == cov, purity == pur)
  plt_process <- df_process %>% 
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
    geom_point(aes(x = name, y = value, col = cluster), size = 4, shape = 20) +
    geom_line(data = ~ filter(.x, !is.na(value)), aes(x = name, y = value, col = cluster, group = cluster), linewidth = .6)  +
    geom_text(data = ~ filter(.x, contains_driver & name == 'score_all'), aes(x = name, y = value+0.03, col = cluster, group = cluster, label = driver_label_process)) + 
    ylab('Score') +
    xlab('')+ 
    scale_color_manual('Cluster', 
                       values = color_palette_process) +
    scale_x_discrete(labels = c('score_driver' = 'Driver',
                                'score_tail'   = 'Tail',
                                'score_sign'   = 'Signature',
                                'score_no_driver' = 'Signature\nTail',
                                'score_no_tail' = 'Driver\nSignature',
                                'score_no_sign' = 'Driver\nTail',
                                'score_all'    = 'All'))  + 
    my_ggplot_theme()  +
    ggtitle(paste0(cov, 'x_',pur, 'p_ProCESS'))
  
  df_tool = readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/res_int/', cov, 'x_', pur, 'p_', s_tool,'_',tool, '_df.rds')) %>% 
    filter(spn == sp)
  
  plt_tool <- df_tool %>%
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
  
  
  class_mutations <- full_class_mutations %>% 
    filter(spn == sp, coverage == cov, purity == pur, sig_tool == s_tool, sub_tool == tool)
  
  plt_score <- list()
  for (score in score_types){
    score_process = paste0('class_all_process')
    score_tool = paste0('class_', score, '_tool')
    
    df <- class_mutations %>%
      count(cluster_id_process, cluster_id_tool, .data[[score_process]], .data[[score_tool]], name = "freq")
    
    df_long <- df %>%
      mutate(flow = row_number()) %>%   # unique id for each flow
      pivot_longer(
        cols = c(cluster_id_process, cluster_id_tool),
        names_to = "axis",
        values_to = "cluster"
      ) %>%
      mutate(
        axis = recode(axis,
                      cluster_id_process = "Process",
                      cluster_id_tool = "Tool"),
        class = if_else(axis == "Process", .data[[score_process]], .data[[score_tool]])
      )
    
    plt_alluvial <- ggplot(df_long,
           aes(x = axis,
               stratum = cluster,
               alluvium = flow,
               y = freq)) +
      
      geom_alluvium(
        aes(fill = .data[[score_process]]),
        width = 0.35,
        alpha = 0.75,
        knot.pos = 0.4,
        curve_type = "sigmoid"
      ) +
      
      geom_stratum(
        aes(fill = class),
        width = 0.35,
        linewidth = 0.1
      ) +
      
      geom_text(
        stat = "stratum",
        aes(label = after_stat(stratum)),
        size = 3,
        lineheight = 0.9
      ) +
      scale_fill_manual(values = color_class) +
      scale_y_continuous(
        labels = scales::comma,
        expand = expansion(mult = c(0.02, 0.12))
      ) +
      labs(
        title    = score,
        y        = "Mutation count",
        x        = NULL,
        fill     = "Process class"
      ) +
      my_ggplot_theme()
    plt_score[[score]] <- plt_alluvial
  }
  
  plt <- plt_process + plt_tool + wrap_plots(plt_score)+ plot_layout(design = 'AB\nCC\nCC')
  ggsave(plot =  plt, 
         filename = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/plot_performance/', sp, '_', cov, 'x_', pur, 'p_', s_tool,'_',tool, '.png'),
         height = 9, width = 11)
}
