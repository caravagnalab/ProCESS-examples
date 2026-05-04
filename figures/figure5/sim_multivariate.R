setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")


color_palette = RColorBrewer::brewer.pal(n = 8, name = "Dark2") 
names(color_palette) <- c('Clonal', paste0('C', 1:7))
color_palette['Subclonal'] = 'gray70'

data_orig <- readRDS("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure5/table_viber_interpreted.rds")

data <- readRDS("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure5/table_viber_interpreted.rds") %>% 
  select(sample_id, mutation_id, VAF, cluster_id_tool, is_driver_true, driver_fake) %>% 
  pivot_wider(
    names_from  = sample_id,
    values_from = VAF
  )

plt_multi <- data %>% 
  ggplot() + 
  geom_point(aes(x = `Sample 1`, y = `Sample 2`, col = cluster_id_tool), size = .5) +
  my_ggplot_theme() +
  scale_color_manual('Cluster', values = color_palette)

plt_multi <- plt_multi + ggrepel::geom_label_repel(
  data = data %>% filter(is_driver_true == TRUE),
  aes(
    x = `Sample 1`,
    y = `Sample 2`,
    label = driver_fake,
    colour = cluster_id_tool,
  ),
  show.legend = F,
  inherit.aes = FALSE,
  size = 3,
  min.segment.length = 0,
  box.padding = 1)

plt_sample1 <- data_orig %>% 
  filter(sample_id == 'Sample 1') %>% 
  ggplot() +
  geom_histogram(aes(x = VAF, fill = cluster_id_tool), alpha=0.6, binwidth = 0.01, position = "identity") +
  xlim(0.05,1.01)+
  xlab(paste0('VAF Sample1')) + 
  my_ggplot_theme() +
  scale_fill_manual('Cluster', values = color_palette) 

plt_sample2 <- data_orig %>% 
  filter(sample_id == 'Sample 2') %>% 
  ggplot() +
  geom_histogram(aes(x = VAF, fill = cluster_id_tool), alpha=0.6, binwidth = 0.01, position = "identity") +
  xlim(0.05,1.01)+
  xlab(paste0('VAF Sample2')) + 
  my_ggplot_theme() +
  scale_fill_manual('Cluster', values = color_palette) 

p = plt_sample1 + theme_minimal() + 
  plt_sample2 + theme_minimal() + 
  plt_multi + theme_minimal() + 
  plot_layout(ncol = 3)

ggsave(filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure5/plot_multi.pdf', 
       plot = p, width = 12, height = 3, dpi = 300, units = 'in')
