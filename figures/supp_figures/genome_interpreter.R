library(dplyr)
library(ggplot2)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

color_class <- c("palegreen4", "goldenrod", "indianred", "gray")
names(color_class) <- c('Tier 1', 'Tier 2', 'Tier 3', 'Tier 4')


pairs <- c('mutect2_ascat', 'mutect2_sequenza', 'strelka_ascat', 'strelka_sequenza')
confusion_per_group_class <- lapply(pairs, FUN = function(p){
  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/interpretation_',p, '/confusion_cluster.rds')) %>% 
    mutate(pair = p) %>% 
    tidyr::separate(pair, into = c('mut_caller', 'cna_caller'), sep = '_')
}) %>% bind_rows()

confusion_per_group_class_mutations <- lapply(pairs,FUN = function(p){
  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/interpretation_',p, '/confusion_mutations.rds')) %>% 
    mutate(pair = p) %>% 
    tidyr::separate(pair, into = c('mut_caller', 'cna_caller'), sep = '_')
}) %>% bind_rows()

cluster <- confusion_per_group_class %>%
  mutate(class = case_when(
    class == 'Clonal' ~ 'Tier 1',
    class == 'Medium\nExpansion' ~ 'Tier 2',
    class == 'Low\nExpansion' ~ 'Tier 3',
    class == 'Neutral\nExpansion' ~ 'Tier 4'
  )) %>% 
  pivot_longer(cols = c(Precision, Recall, Specificity, F1, Accuracy)) %>%
  filter(name == 'Accuracy') %>% 
  filter(score_type == 'all') %>% 
  mutate(type = 'Cluster')

mutations <- confusion_per_group_class_mutations %>% 
  mutate(class = case_when(
    class == 'Clonal' ~ 'Tier 1',
    class == 'Medium\nExpansion' ~ 'Tier 2',
    class == 'Low\nExpansion' ~ 'Tier 3',
    class == 'Neutral\nExpansion' ~ 'Tier 4'
  )) %>% 
  pivot_longer(cols = c(Precision, Recall, Specificity, F1, Accuracy)) %>%
  filter(name == 'Accuracy') %>% 
  filter(score_type == 'all') %>% 
  mutate(type = 'Mutations')

plt_cluster <- cluster %>% 
  ggplot() +
  geom_boxplot(aes(x = class, y = value, col = class, fill = class), alpha = .6) +
  scale_color_manual('ProCESS Tier', values = color_class) +
  scale_fill_manual('ProCESS Tier', values = color_class) +
  my_ggplot_theme() +
  ggh4x::facet_nested(coverage+purity~sub_tool+mut_caller+cna_caller) +
  ylab('Tier Assignment\nAccuracy')  +
  xlab('') + 
  ylim(0,1) +
  ggtitle('Cluster')


plt_mutatations <- mutations %>% 
  ggplot() +
  geom_boxplot(aes(x = class, y = value, col = class, fill = class), alpha = .6) +
  scale_color_manual('ProCESS Tier', values = color_class) +
  scale_fill_manual('ProCESS Tier', values = color_class) +
  my_ggplot_theme() +
  ggh4x::facet_nested(coverage+purity~sub_tool+mut_caller+cna_caller) +
  ylab('Tier Assignment\nAccuracy')  +
  xlab('') + 
  ylim(0,1) +
  ggtitle('Mutations')

p <- plt_cluster + plt_mutatations + plot_layout(nrow = 2, guides = 'collect')
ggsave("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/supp_figures/Final_GenomeInterpreter_SCOUT_Validation.pdf",
       width = 11, height = 12)
