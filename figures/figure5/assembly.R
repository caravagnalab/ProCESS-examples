source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/SCOUT/colors.R')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils.R")

color_class <- c("palegreen4","darkseagreen", "goldenrod", "indianred", "gray")
names(color_class) <- c('Clonal', 'Strong\nExpansion', 'Medium\nExpansion', 'Low\nExpansion', 'Neutral\nExpansion')


color_score <- c('plum4', 'darkseagreen4', 'cadetblue4', 'salmon2')

pairs <- c('mutect2_ascat','strelka_ascat', 'mutect2_sequenza')
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
  pivot_longer(cols = c(Precision, Recall, Specificity, F1, Accuracy)) %>%
  filter(name == 'Accuracy') %>% 
  filter(score_type == 'all') %>% 
  mutate(type = 'Cluster')

mutations <- confusion_per_group_class_mutations %>% 
  pivot_longer(cols = c(Precision, Recall, Specificity, F1, Accuracy)) %>%
  filter(name == 'Accuracy') %>% 
  filter(score_type == 'all') %>% 
  mutate(type = 'Mutations')


  
plt <- cluster %>% 
  bind_rows(mutations) %>% 
  ggplot() +
  geom_boxplot(aes(x = class, y = value, col = class, fill = class), alpha = .6) +
  scale_color_manual('ProCESS Class', values = color_class) +
  scale_fill_manual('ProCESS Class', values = color_class) +
  my_ggplot_theme() +
  facet_grid(type~.) +
  ylab('Class Assignment\nAccuracy')  +
  xlab('') + 
  ylim(0,1)

p_all <-  ggplot() + labs(tag = "A") + 
  plt + labs(tag = "B") +
  ggplot() + labs(tag = "C") +
  ggplot()+  labs(tag = "D") + 
  ggplot()+  labs(tag = "E") + 
  plot_layout( design = 'AAABBB\nAAABBB\nCCDDEE\nCCDDEE\nCCDDEE') & theme(legend.position = 'bottom')

ggsave(filename = paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure5/fig5_main.pdf"), 
       plot = p_all, device="pdf", width=6.5, height=9, units="in")
