setwd('/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure2')
source('stats_plot.R')

setwd('/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure2')
ht <- readRDS('signature.rds')
ht_grob <- grid.grabExpr(draw(ht))


free(plt_muts_count + labs(tag = "A")) +  
  free(stats_space_plot1 + labs(tag = "B")) + 
  free(wrap_elements(ht_grob) + labs(tag = "D")) + 
  free(ggplot() + labs(tag = "C")) + plot_layout(design = 'aabbdd\nccccdd') & theme(legend.position = "bottom")

base = "/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/"
ggsave(filename = paste0(base, 'fig2_main.pdf'), 
       width = 260, height = 150, units = 'mm', dpi = 300) #device = cairo_pdf()
