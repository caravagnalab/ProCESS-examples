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



df = df %>% 
  filter(!is.na(CosineSimilarity)) %>% 
  mutate(cluster_process = ifelse(cluster_process == 'Clonal', 'Clonal', 'Subclonal')) %>% 
  filter(!is.na(cluster_process)) 


# 1. Run tests per cluster
stat.test <- df %>%
  rstatix::wilcox_test(CosineSimilarity ~ cluster_process) %>%
  rstatix::adjust_pvalue(method = "BH") %>%       
  rstatix::add_significance("p.adj") %>%          
  rstatix::add_xy_position(x = "cluster_process")

# 2. Plot
pd <- position_dodge(0.8)

sign_cluster <- df %>% 
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

#ggsave('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature/cluster_sig.pdf',
#       width = 4, height = 3, units = 'in')


# nmi <- readRDS('/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/metrics_tables/new_table_clusters_metrics_final.rds')
# nmi <- nmi %>% select(spn, purity, coverage, nmi_raw, nmi_interpreted, nmi_interpreted_driver, tool)
# 
# merge <- df %>% left_join( nmi %>% dplyr::rename(pur = purity, cov = coverage))
# 
# merge %>% 
#   ggplot() +
#   geom_point(aes(x = nmi_raw, y = CosineSimilarity, col = tool)) +
#   #geom_smooth(aes(x = nmi_raw, y = CosineSimilarity, col = tool, fill = tool), method = 'glm') + 
#   scale_color_manual(values = c('plum4', 'aquamarine4')) + 
#   scale_fill_manual(values = c('plum4', 'aquamarine4')) + 
#   theme_bw() 
