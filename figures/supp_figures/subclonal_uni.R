library(dplyr)
library(ggplot2)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/subclone_mobster.R')
plt_subclone <- results %>% 
  filter(!is.na(gt)) %>% 
  mutate(TP = gt & pred,
         FP = !gt & pred,
         FN = gt & !pred,
         TN = !gt & !pred) %>% 
  group_by(subclass, purity, coverage, vcf_caller, cna_caller) %>% 
  summarise(Polyclonal = sum(TP) / (sum(TP) + sum(FN)),  #recall
            Monoclonal = sum(TN) / (sum(TN) + sum(FP))) %>%  #Specificity
  pivot_longer(cols = c(Polyclonal, Monoclonal),
               names_to = "metric",
               values_to = "value") %>%
  ggplot(aes(x = as.factor(purity), 
             y = value, 
             color = subclass, 
             group = subclass)) +
  geom_line() +
  geom_point(size = 2) +
  scale_color_manual('Subclass', values = c('Polyclonal' = 'goldenrod', 
                                            'Monoclonal'='#645394', 
                                            'Weak evidence' = 'khaki3',
                                            'Sampling Bias' = 'plum')) +
  ylab("Sample composition detection\nTPR") +
  xlab("Purity") +
  my_ggplot_theme() +
  ylim(0,1) + 
  ggh4x::facet_nested(vcf_caller + cna_caller~metric + coverage )

ggsave("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/supp_figures/Final_SD_Univariate_SCOUT_Validation.pdf",
    width = 9, height = 5)
