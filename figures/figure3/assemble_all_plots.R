### ASSEMBLE FINAL PANEL
base = "/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/"
source(paste0(base, 'figure3/cna/cna.R'))
source(paste0(base, 'figure3/somatic/somatic.R'))
base = "/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/"
source(paste0(base, 'figure3/segm_cna/plot_driver_cna.R'))

cartoon <- ggplot() + labs(tag = "A")
caller_line_boxes <- caller_line_boxes 
scatter_purity_ploidy <- scatter_purity_ploidy + labs(tag = "B") 
alluvial_plot <- alluvial_plot + labs(tag = "C") 
plt_breakpoint <- plt_breakpoint + labs(tag = "D")
driver_gene_plot <- driver_gene_plot + labs(tag = "E")


wrap_plots(
  free(cartoon),
  free(caller_line_boxes), 
  free(scatter_purity_ploidy),
  free(alluvial_plot),
  free(plt_breakpoint),
  free(driver_gene_plot),
  design="aaaccc\nbbbccc\nbbbccc\nbbbddd\nbbbddd\neeefff\neeefff\neeefff"
)


ggsave(filename = paste0(base, 'fig3_main.pdf'), 
        width = 260, height = 210, units = 'mm', dpi = 300) #device = cairo_pdf()
