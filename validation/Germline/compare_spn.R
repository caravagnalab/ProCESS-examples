rm(list = ls())
options(bitmapType='cairo')
library(tidyverse)
library(vcfR)
library(optparse)
library(caret)
library(dplyr)
library(patchwork)

source('utils.R')
source('../../getters/sarek_getters.R')
source('../../getters/process_getters.R')

dir <- '/orfeo/scratch/cdslab/shared/SCOUT/'

all_metric <- tibble()
all_baf <- tibble()

spn <- paste0('SPN0', 1:7)
for (s in spn){

  input <- paste0(dir, s, "/validation/germline/report")
  if (file.exists(file.path(input, 'normal_metrics.rds'))){
    
    rds <- readRDS(file.path(input, 'normal_metrics.rds'))
    metric <- rds$report_metrics
    baf <- rds$baf_metric

    all_metric <- bind_rows(all_metric, metric) 
    all_baf <- bind_rows(all_baf, baf) 
    
  }
}

colors = c("ProCESS"="deepskyblue3", "haplotypecaller"="coral3", "strelka"="palegreen4", "freebayes"="maroon")
plt_metric <- all_metric %>% 
  pivot_longer(cols = c(Accuracy, Sensitivity, Precision, Recall, F1_Score)) %>% 
  ggplot() +
  geom_point(aes(x = spn, y = value, col = tool), size = 3) +
  scale_color_manual(values = colors) +
  xlab('metric') + 
  facet_grid(name~.) + 
  theme_bw()

plt_baf <- all_baf %>% 
  ggplot() +
  geom_boxplot(aes(x = spn, y = RMSE, fill = tool)) +
  scale_fill_manual(values = colors) +
  theme_bw() + 
  all_baf %>% 
  ggplot() +
  geom_boxplot(aes(x = spn, y = cor_coeff, fill = tool)) +
  scale_fill_manual(values = colors) +
  theme_bw() +
  plot_layout(guides = 'collect')
plt_baf
