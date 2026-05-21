library(tidyr)
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/SCOUT/colors.R')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils.R")

color_class <- c("palegreen4", "goldenrod", "indianred", "gray")
names(color_class) <- c('Tier1', 'Tier2', 'Tier3', 'Tier4')

color_score <- c('plum4', 'darkseagreen4', 'cadetblue4', 'salmon2')

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
    class == 'Clonal' ~ 'Tier1',
    class == 'Medium\nExpansion' ~ 'Tier2',
    class == 'Low\nExpansion' ~ 'Tier3',
    class == 'Neutral\nExpansion' ~ 'Tier4'
  )) %>% 
  pivot_longer(cols = c(Precision, Recall, Specificity, F1, Accuracy)) %>%
  filter(name == 'Accuracy') %>% 
  filter(score_type == 'all') %>% 
  mutate(type = 'Cluster')

mutations <- confusion_per_group_class_mutations %>% 
  mutate(class = case_when(
    class == 'Clonal' ~ 'Tier1',
    class == 'Medium\nExpansion' ~ 'Tier2',
    class == 'Low\nExpansion' ~ 'Tier3',
    class == 'Neutral\nExpansion' ~ 'Tier4'
  )) %>% 
  pivot_longer(cols = c(Precision, Recall, Specificity, F1, Accuracy)) %>%
  filter(name == 'Accuracy') %>% 
  filter(score_type == 'all') %>% 
  mutate(type = 'Mutations')


kw_test <- mutations %>%
  kruskal_test(value ~ class)

# Step 2 — pairwise Wilcoxon post-hoc
stat_test <- mutations %>%
  wilcox_test(value ~ class, comparisons = NULL) %>%  # NULL = all pairs
  adjust_pvalue(method = "BH") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "class")
  
plt <- #cluster %>% 
  mutations %>% 
  ggplot() +
  geom_boxplot(aes(x = class, y = value, col = class, fill = class), alpha = .6) +
  scale_color_manual('Ground truth', values = color_class) +
  scale_fill_manual('Ground truth', values = color_class) +
  stat_pvalue_manual(
    stat_test,
    label = "p.adj.signif",
    #hide.ns = TRUE,
    tip.length = 0.01,
    size = 3
  ) +
  my_ggplot_theme() +
  #ggh4x::facet_grid2(type~.) +
  #ggh4x::facet_grid2(type~sub_tool+mut_caller) +
  ylab('Class Assignment\nAccuracy')  +
  xlab('') + 
  ylim(0,1)



ggsave(filename = paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure5/fig5_boxplot.pdf"), 
       plot = plt, device="pdf", width=3, height=2, units="in")

# source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/spn02/edf_spn02.R')
# 
# free(p_tool + theme(legend.position = 'none')) + plt + plot_layout(nrow = 2)

p_all <-  ggplot() + labs(tag = "A") + 
  plt + labs(tag = "B") +
  ggplot() + labs(tag = "C") +
  ggplot()+  labs(tag = "D") + 
  ggplot()+  labs(tag = "E") + 
  plot_layout( design = 'AAABBB\nAAABBB\nCCDDEE\nCCDDEE\nCCDDEE') & theme(legend.position = 'bottom')

ggsave(filename = paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure5/fig5_main.pdf"), 
       plot = p_all, device="pdf", width=6.5, height=9, units="in")
