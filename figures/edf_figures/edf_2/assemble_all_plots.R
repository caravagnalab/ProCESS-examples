library(dplyr)
library(lubridate)
library(ggplot2)
source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_2/stats_plot.R")
source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_2/timelime_scout.R")


plt <- wrap_plots(list(plt_muts_count,
                       stats_space_plot1,
                       plt_gannt + guides(colour = guide_legend(nrow = 3, bycol = TRUE))
                       ), 
                  design = "AABB\nCCCC") + 
  plot_annotation(tag_levels = 'A') & 
  theme(legend.position = "bottom")

ggsave(plot = plt, 
       filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_2/edf_space.pdf', 
       width = 8, height = 7, units = 'in')






