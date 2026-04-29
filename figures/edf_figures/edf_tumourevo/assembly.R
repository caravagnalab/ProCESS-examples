
source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/univariate_sampling_bias/univariate_edf.R")
source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_5/spn03_edf.R")
source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_5/spn07_edf.R")

plt_spn03 <- spn03 #+ labs(tag = "A")
plt_spn07 <- spn07 #+ labs(tag = "B")
plt_sb <- wrap_plots(cluster_plots_samples, nrow = 2) + labs(tag = 'C')

p <- wrap_plots(
  plt_spn03,
  plt_spn07,
  plt_sb,
  design = 'AABB\nCCCC')

setwd('/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_tumourevo/')
ggsave(plot = p, 
       filename = 'edf_tumourevo.pdf',
       width = 8.5, height = 9, units = 'in', dpi = 300) #device = cairo_pdf()
