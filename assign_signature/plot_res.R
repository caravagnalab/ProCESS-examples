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
# 
# df_uni %>% 
#   filter(cluster != 'Tail-Subclonal') %>% 
#   ggplot() +
#   geom_boxplot(aes(x = sample_class, y = CosineSimilarity, fill = cluster_process, col = cluster_process), alpha = .4) + 
#   #scale_fill_manual('Cluster Type', values = c('plum4', 'cadetblue4')) + 
#   #scale_color_manual('Cluster Type', values = c('plum4', 'cadetblue4')) + 
#   my_ggplot_theme() +
#   xlab('Purity') +
#   ggtitle('Univariate MOBSTER')


# pd <- position_dodge(width = 0.3)
# 
# plt_uni <- df_uni %>%
#   group_by(pur, cluster_process, sample_class) %>%
#   summarise(
#     q25 = quantile(CosineSimilarity, 0.25, na.rm = TRUE),
#     q50 = quantile(CosineSimilarity, 0.50, na.rm = TRUE),
#     q75 = quantile(CosineSimilarity, 0.75, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   ggplot(aes(
#     x = as.factor(pur),
#     y = q50,
#     group = interaction(cluster_process, sample_class),
#     color = cluster_process,
#     linetype = sample_class,
#     shape = sample_class
#   )) +
#   geom_line(position = pd, linewidth = .3) +
#   geom_point(position = pd, size = 3) +
#   geom_errorbar(
#     aes(ymin = q25, ymax = q75),
#     position = pd,
#     width = 0.15,
#     alpha = 0.6
#   ) +
#   scale_color_manual(
#     "Cluster Type",
#     values = c('Clonal' = "plum4", 'Subclonal' =  "cadetblue4")
#   ) +
#   scale_shape_manual('Sample Class', values = c(16,17)) + 
#   scale_linetype_manual('Sample Class', values = c(1,2)) + 
#   my_ggplot_theme() +
#   xlab("Purity") +
#   ylab("Median Cosine Similarity") +
#   ggtitle("Univariate")
# 
# plt_uni

df_uni_filt <- df_uni %>% filter(cluster != 'Tail-Subclonal')
df = df %>% 
  bind_rows(df_uni_filt) %>% 
  filter(!is.na(CosineSimilarity)) %>% 
  mutate(cluster_process = ifelse(cluster_process == 'Clonal', 'Clonal', 'Subclonal')) %>% 
  filter(!is.na(cluster_process)) %>% 
  mutate(patient_class = ifelse(spn %in% c('SPN03', 'SPN04', 'SPN06', 'SPN07'), 'High\nComplexity', 'Low\nComplexity')) 

# 1. Run tests per cluster
stat.test <- df %>%
  group_by(cluster_process) %>%
  wilcox_test(CosineSimilarity ~ patient_class) %>%
  rstatix::adjust_pvalue(method = "BH") %>%       
  rstatix::add_significance("p.adj") %>%          
  rstatix::add_xy_position(x = "cluster_process") 

# 2. Plot
pd <- position_dodge(0.8)


stat.within <- df %>%
  group_by(patient_class) %>%
  wilcox_test(CosineSimilarity ~ cluster_process) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "patient_class", dodge = 0.8)%>% mutate(y.position  = 1.02)

pairwise_by_process <- function(df, process, y_offset) {
  df %>%
    filter(cluster_process == process) %>%
    pairwise_wilcox_test(
      CosineSimilarity ~ patient_class,
      p.adjust.method = "BH"
    ) %>%
    add_significance("p.adj") %>%
    mutate(
      y.position = max(df$CosineSimilarity, na.rm = TRUE) * y_offset
    )
}

stat.subclonal <- pairwise_by_process(
  df, process = "Subclonal", y_offset = 1.06
)

stat.clonal <- pairwise_by_process(
  df, process = "Clonal", y_offset = 1.10
)


df$patient_class <- factor(
  df$patient_class,
  levels = c("High\nComplexity", "Low\nComplexity")
)

df$cluster_process <- factor(
  df$cluster_process,
  levels = c( "Clonal", "Subclonal")
)

plt_multi <- df %>%
  ggplot(aes(
    y = CosineSimilarity,
    x = patient_class,
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
  #hide.ns = TRUE
) +

## between-class (Subclonal only)
stat_pvalue_manual(
  stat.subclonal,
  label = "p.adj.signif",
  tip.length = 0.01,
  inherit.aes = FALSE,
  color = "#9F4576"
) +

## between-class (Clonal only)
stat_pvalue_manual(
  stat.clonal,
  label = "p.adj.signif",
  tip.length = 0.01,
  inherit.aes = FALSE,
  color = "palegreen4"
) +  scale_y_continuous(
  breaks = c(0.3, 0.6, 0.9, 1),
  labels = c("0.3", "0.6", "0.9", "1")
) +

  scale_color_manual(
    "Cluster Type",
    values = c('Clonal' = "palegreen4", 'Subclonal' =  "#9F4576")
  ) +
  scale_fill_manual(
    "Cluster Type",
    values = c('Clonal' = "palegreen4", 'Subclonal' =  "#9F4576")
  ) +
  my_ggplot_theme() +
  xlab("Patient Class") +
  ylab('Exposure Accuracy (Cosine Similarity)')


#plt_multi
#ggsave(filename = 'plot.pdf', width = 4, height = 4, units = 'in')
