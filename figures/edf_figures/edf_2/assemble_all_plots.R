library(dplyr)
library(lubridate)
library(ggplot2)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/edf_figures/edf_2/stats_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/edf_figures/edf_2/timelime_scout.R")


wrap_plots(list(plt_muts_count,stats_space_plot1,plt_gannt),design = "AABB\nAABB\nCCCC") + 
  plot_annotation(tag_levels = 'A') & theme(legend.position = "bottom")
