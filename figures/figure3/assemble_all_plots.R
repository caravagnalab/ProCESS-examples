### ASSEMBLE FINAL PANEL
base = "/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/"
source(paste0(base, 'figure3/cna/cna.R'))
source(paste0(base, 'figure3/somatic/somatic.R'))
base = "/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/"
source(paste0(base, 'figure3/segm_cna/panel_E.R'))
source(paste0(base, 'figure3/driver_cna/driver_cna_final.R'))

caller_line_boxes <- caller_line_boxes + labs(tag = "A")
plt_purity_ploidy <- plt_purity_ploidy + labs(tag = "B")
cna_corr <- cna_corr + labs(tag = "C")
plt_breakpoint <- plt_breakpoint + labs(tag = "D")
driver_cna_vaf <- driver_cna_vaf + labs(tag = "E")

# wrap_plots(
#   wrap_plots(caller_line_boxes, purity_line_boxes, coverage_line_boxes, ncol=1, guides = 'collect'),
#   plt_purity_ploidy,
#   plt_breakpoint,
#   driver_cna_vaf,
#   design="aaaabbb
#           aaaabbb
#           aaaaccd
#           "
# )  

wrap_plots(
  free(wrap_plots(caller_line_boxes, purity_line_boxes, coverage_line_boxes, ncol=1,guides = 'collect')),
  free(plt_purity_ploidy), 
  free(cna_corr),
  free(plt_breakpoint),
  free(driver_cna_vaf),
  design="aaaaccc
          aaaaccc
          aaaaddd
          aaaaddd
          aaaaddd
          bbbbeee
          bbbbeee
          bbbbeee
          "
)

ggsave(filename = paste0(base, 'fig3_main.pdf'), 
       width = 230, height = 180, units = 'mm', dpi = 300,device = cairo_pdf())
