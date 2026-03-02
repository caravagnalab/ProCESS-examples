library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")
source('/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/figures/figure3/utils_plot.R')

panelB = readRDS(file = paste0(save_path, "plots/metrics/cluster_analysis/rds/precision_blind_vs_interpreted_w_tool.rds"))

panelC = readRDS(file = paste0(save_path, "plots/metrics/driver_analysis/rds/precision_recall_clonal_vs_subclonal_only_TPc.rds"))

final_plot = patchwork::wrap_plots(
  ggplot(),
  panelB,
  panelC,
  design="abb\nacc")&
  my_ggplot_theme()&
  plot_annotation(tag_levels = 'A')

final_plot

ggsave(filename = paste0(save_path, "plots/figure4.pdf"), plot = final_plot,
       device="pdf", width=6, height=5, units="in")

