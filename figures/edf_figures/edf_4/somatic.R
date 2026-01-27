source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/somatic/somatic.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/germline/germline.R')


coverage_line_boxes  = coverage_line_boxes + labs(tag = "A") 
purity_line_boxes = purity_line_boxes + labs(tag = "B") 
caller_line_germline = caller_line_germline + guides(
    col  = guide_legend(nrow = 3),
    fill = guide_legend(nrow = 3)
  ) + theme(legend.position = 'bottom')+ labs(tag = "C") 
empty = ggplot() + labs(tag = "D") 
cn_heatmap = ggplot() + labs(tag = "E") 


wrap_plots(
  free(coverage_line_boxes),
  free(purity_line_boxes), 
  free(caller_line_germline),
  free(empty),
  free(cn_heatmap),
  design="aac\nbbd\neee"
)

ggsave(filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_4/somatic.pdf', 
       width = 210, height = 250, units = 'mm', dpi = 300)
