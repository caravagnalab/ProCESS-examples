library(dplyr)
library(ggplot2)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

info <- tibble(spn = paste0('SPN0',1:7),
               class = c('Other', 'Hypermutant', 'Other', 'Other', 'Hypermutant', 'Other', 'Hypermutant'))

data <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/drivers/performance_driver.rds') %>% 
  left_join(info)

plt_driver <- data %>% 
  mutate(class = factor(class, levels = c('Hypermutant', 'Other'))) %>% 
  ggplot(aes(x = class, y = F1)) +
  #geom_boxplot(aes(x = class, y = F1), alpha = .4) + 
  geom_boxplot(outliers = F, col = 'gray40', fill = 'gainsboro', alpha = .4) +
  stat_summary(
    aes(x = class, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3),
    size = .4, show.legend = T
  ) +
  stat_summary(
    aes(x = class, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3),
    geom = 'line',
    linewidth=.35, show.legend = F
  ) +
  my_ggplot_theme() +
  xlab('Mutation Rate Class') +
  ylab('Driver Detection\nF1 Score') +
  scale_color_manual(values=SPN_colors, name='SPN')

plt_driver
