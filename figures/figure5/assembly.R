source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/SCOUT/colors.R')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils.R")

color_class <- c("palegreen4","darkseagreen", "goldenrod", "indianred", "gray")
names(color_class) <- c('Clonal', 'Strong\nExpansion', 'Medium\nExpansion', 'Low\nExpansion', 'Neutral\nExpansion')


color_score <- c('plum4', 'darkseagreen4', 'cadetblue4', 'salmon2')

confusion_per_group_class <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/interpretation/confusion_cluster.rds')
confusion_per_group_class_mutations <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/interpretation/confusion_mutations.rds')

all_metrics <- confusion_per_group_class %>% 
  pivot_longer(cols = c(Precision, Recall, Specificity, F1, Accuracy)) %>% 
  ggplot()+
  geom_boxplot(aes(x = class, y = value, col = class, fill = class), alpha = .6) +
  scale_color_manual('ProCESS Class', values = color_class) +
  scale_fill_manual('ProCESS Class', values = color_class) +
  my_ggplot_theme() +
  facet_grid(name~.)  +
  ggtitle('Clusters') + 
  
confusion_per_group_class_mutations %>% 
  pivot_longer(cols = c(Precision, Recall, Specificity, F1, Accuracy)) %>% 
  ggplot()+
  geom_boxplot(aes(x = class, y = value, col = class, fill = class), alpha = .6) +
  scale_color_manual('ProCESS Class', values = color_class) +
  scale_fill_manual('ProCESS Class', values = color_class) +
  my_ggplot_theme() +
  facet_grid(name~.) +
  ggtitle('Mutations') 
ggsave(filename = paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/interpretation/plot_metrics.pdf"), 
       plot = all_metrics, device="pdf", width=8.5, height=6, units="in")


all_metrics_score <- confusion_per_group_class %>% 
  pivot_longer(cols = c(Precision, Recall, Specificity, F1, Accuracy)) %>% 
  ggplot()+
  geom_boxplot(aes(x = class, y = value, col = score_type, fill = score_type), alpha = .6) +
  scale_color_manual('Score', values = color_score) +
  scale_fill_manual('Score', values = color_score) +
  my_ggplot_theme() +
  facet_grid(name~.)  +
  ggtitle('Clusters') + 
  
  confusion_per_group_class_mutations %>% 
  pivot_longer(cols = c(Precision, Recall, Specificity, F1, Accuracy)) %>% 
  ggplot()+
  geom_boxplot(aes(x = class, y = value, col = score_type, fill = score_type), alpha = .6) +
  scale_color_manual('Score', values = color_score) +
  scale_fill_manual('Score', values = color_score) +
  my_ggplot_theme() +
  facet_grid(name~.) +
  ggtitle('Mutations') 

ggsave(filename = paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/interpretation/plot_score.pdf"), 
       plot = all_metrics_score, device="pdf", width=8.5, height=6, units="in")


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
