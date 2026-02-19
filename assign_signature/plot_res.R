library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)

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



df = df %>% 
  filter(!is.na(CosineSimilarity)) %>% 
  mutate(cluster_process = ifelse(cluster_process == 'Clonal', 'Clonal', 'Subclonal')) %>% 
  filter(!is.na(cluster_process)) 

df_uni = df_uni %>% 
  filter(!is.na(CosineSimilarity)) %>% 
  mutate(cluster_process = ifelse(cluster_process == 'Clonal', 'Clonal', 'Subclonal')) %>% 
  filter(!is.na(cluster_process)) 

df_all = df %>% 
  bind_rows(df_uni) %>% 
  mutate(class = ifelse(spn %in% c('SPN03', 'SPN06', 'SPN07'), 'High\nComplexity', 'Low\nComplexity'))

# 1. Run tests per cluster
stat.test <- df_all %>%
  rstatix::wilcox_test(CosineSimilarity ~ cluster_process) %>%
  rstatix::adjust_pvalue(method = "BH") %>%       
  rstatix::add_significance("p.adj") %>%          
  rstatix::add_xy_position(x = "cluster_process")

# 2. Plot
pd <- position_dodge(0.8)

sign_cluster <- df_all %>% 
  ggplot(aes(y = CosineSimilarity,
             x = cluster_process,
             fill = cluster_process,
             colour = cluster_process)) +
  geom_boxplot(position = pd,
               alpha = 0.2,
               outlier.size = 1, show.legend = F) +
  stat_pvalue_manual(
    stat.test,
    label = "p.adj.signif",   # "***", "**", etc.
    tip.length = 0.01
  ) +
  scale_color_manual('Cluster Type', values = c('lightcoral', 'cadetblue4')) + 
  scale_fill_manual('Cluster Type', values = c('lightcoral', 'cadetblue4')) + 
  my_ggplot_theme()+
  xlab('Cluster Type')  #+
  #coord_flip()

stat.within <- df_all %>%
  group_by(class) %>%
  wilcox_test(CosineSimilarity ~ cluster_process) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "class", dodge = 0.8)

pairwise_by_process <- function(df, process, y_offset) {
  df %>%
    filter(cluster_process == process) %>%
    pairwise_wilcox_test(
      CosineSimilarity ~ class,
      p.adjust.method = "BH"
    ) %>%
    add_significance("p.adj") %>%
    mutate(
      y.position = max(df_all$CosineSimilarity, na.rm = TRUE) * y_offset
    )
}

stat.subclonal <- pairwise_by_process(
  df_all, process = "Subclonal", y_offset = 1.05
)

stat.clonal <- pairwise_by_process(
  df_all, process = "Clonal", y_offset = 1.10
)

df_all %>%
  ggplot(aes(
    y = CosineSimilarity,
    x = class,
    fill = cluster_process,
    colour = cluster_process
  )) +
  geom_boxplot(
    position = pd,
    alpha = 0.2,
    outlier.size = 1,
    show.legend = TRUE
  ) +
  
  ## within-class (Clonal vs Subclonal)
  stat_pvalue_manual(
    stat.within,
    label = "p.adj.signif",
    tip.length = 0.01,
    hide.ns = TRUE
  ) +
  
  ## between-class (Subclonal only)
  stat_pvalue_manual(
    stat.subclonal,
    label = "p.adj.signif",
    tip.length = 0.01,
    inherit.aes = FALSE,
    color = "cadetblue4"
  ) +
  
  ## between-class (Clonal only)
  stat_pvalue_manual(
    stat.clonal,
    label = "p.adj.signif",
    tip.length = 0.01,
    inherit.aes = FALSE,
    color = "lightcoral"
  ) +
  
  scale_color_manual(
    "Cluster Type",
    values = c("lightcoral", "cadetblue4")
  ) +
  scale_fill_manual(
    "Cluster Type",
    values = c("lightcoral", "cadetblue4")
  ) +
  my_ggplot_theme() +
  xlab("Patient Class")

ggsave(filename = 'plot.pdf', width = 3, height = 4, units = 'in')
