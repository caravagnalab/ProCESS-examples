library(dplyr)
library(ggplot2)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/multivariate.R')
plt_subclone <- results %>% 
  pivot_longer(cols = c(viber, pyclone)) %>% 
  ggplot() +
  geom_boxplot(aes(x = as.factor(purity), y = value, col = type, fill = type), alpha =.2) +
  my_ggplot_theme() + 
  xlab('Purity') +
  scale_fill_manual('SPN', values = c('darkseagreen4', 'orangered3')) + 
  scale_color_manual('SPN', values = c('darkseagreen4', 'orangered3')) + 
  ylab('Relative number of inferred cluster\nClusters / Samples') + 
  ggh4x::facet_nested( vcf_caller + cna_caller ~ name + coverage)

ggsave("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/supp_figures/Final_SD_Multivariate_SCOUT_Validation.pdf",
       width = 9, height = 5)
