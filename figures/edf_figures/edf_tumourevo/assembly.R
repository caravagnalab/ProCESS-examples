
source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/univariate_sampling_bias/univariate_edf.R")
source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_5/spn03_edf.R")

plt_cartoon <- ggplot() + labs(tag = "A")
plt_spn03 <- spn03
plt_sb <- wrap_plots(cluster_plots_samples, nrow = 4) + labs(tag = 'C')

p <- wrap_plots(
  free(plt_cartoon),
  free(plt_spn03),
  free(plt_sb),
  design = 'AC\nBC\nBC')

setwd('/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_tumourevo/')
ggsave(plot = p, 
       filename = 'edf_tumourevo.pdf',
       width = 9.5, height = 7, units = 'in', dpi = 300) #device = cairo_pdf()
