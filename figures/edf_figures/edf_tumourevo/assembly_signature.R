
source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_5/spn03_edf.R")
spn03 <- wrap_plots(list(pp2+theme(legend.position = 'none'),
                         pp1+theme(legend.position = 'none'),
                         pp3+theme(legend.position = 'none'),
                         p_sankey+theme(strip.text.x = element_blank(),strip.text.y = element_blank(), legend.position = 'none')),
                    design = "AAAA\nBBBC\nBBBC\nDDDD")

source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_tumourevo/phylo_03.R")


plt_spn03 <- wrap_plots(ggplot(), phylo_spn03 + theme(legend.position = 'none'), nrow = 1) #+ plot_annotation(tag_levels = 'A')
setwd('/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_tumourevo/')
#plt_spn03 <- ggrastr::rasterize(plt_spn03, layers=c('Point',  'Edge'), dpi=300)
ggsave(plot = plt_spn03, 
       filename = 'edf_sig_spn03.png',
       width = 8.5, height = 4.5, units = 'in', dpi = 600) #device = cairo_pdf()


source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_5/spn07_edf.R")
source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_tumourevo/phylo_07.R")
plt_spn07 <- wrap_plots(ggplot(), phylo_spn07 + theme(legend.position = 'none'), nrow = 1) #+ plot_annotation(tag_levels = 'A')
setwd('/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_tumourevo/')
ggsave(plot = plt_spn07, 
       filename = 'edf_sig_spn07.png',
       width = 8.5, height = 4.5, units = 'in', dpi = 600) #device = cairo_pdf()


