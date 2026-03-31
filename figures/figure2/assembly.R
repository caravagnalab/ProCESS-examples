setwd('/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure2')
source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/tmb/tmb.R")
source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/cn_heatmap/cn_profile.R")

plt_tmb <- plt_tmb + labs(tag = "A")
plt_onco <- ggplot() + labs(tag = "B")
cn_plot <- cn_plot + labs(tag = 'C')

p <- wrap_plots(
  free(plt_tmb),
  free(plt_onco),
  free(cn_plot),
  design = 'AABBB\nAABBB\nAABBB\nAABBB\nCCCCC\nCCCCC')


ggsave(plot = p, filename = paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure2/", 'fig2_main.pdf'), 
       width = 10, height = 10, units = 'in', dpi = 300) #device = cairo_pdf()
