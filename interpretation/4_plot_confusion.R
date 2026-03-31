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

source('../validation/SCOUT/colors.R')
source("../figures/figure3/utils_plot.R")
source("../figures/figure3/utils.R")

color_class <- c("palegreen4","darkseagreen", "goldenrod", "indianred", "gray")
names(color_class) <- c('Clonal', 'Strong\nExpansion', 'Medium\nExpansion', 'Low\nExpansion', 'Neutral\nExpansion')
color_score <- c('plum4', 'darkseagreen4', 'cadetblue4', 'salmon2')

vcf_caller = "mutect2"
cna_caller = "sequenza"
base = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/interpretation_', vcf_caller, "_", cna_caller, '/')


confusion_per_group_class <- readRDS(paste0(base, '/confusion_cluster.rds'))
confusion_per_group_class_mutations <- readRDS(paste0( base, 'confusion_mutations.rds'))

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
ggsave(filename = paste0(base, "/plot_metrics.pdf"), 
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

ggsave(filename = paste0(base, "/plot_score.pdf"), 
       plot = all_metrics_score, device="pdf", width=8.5, height=6, units="in")
