### ASSEMBLE FINAL PANEL

caller_line_boxes <- caller_line_boxes+labs(tag = "A")
plt_purity_ploidy <- plt_purity_ploidy+labs(tag = "B")
plt_breakpoint <- plt_breakpoint+labs(tag = "C")
driver_cna_vaf <- driver_cna_vaf+labs(tag = "D")

wrap_plots(
  wrap_plots(caller_line_boxes, purity_line_boxes, coverage_line_boxes, ncol=1),
  plt_purity_ploidy, plt_breakpoint,driver_cna_vaf,
  design="aaaabbb
          aaaabbb
          aaaaccd
          "
)  

wrap_plots(
  wrap_plots(caller_line_boxes, purity_line_boxes, coverage_line_boxes, ncol=1),
  plt_purity_ploidy, plt_breakpoint,driver_cna_vaf,
  design="aaaaccc
          aaaaccc
          aaaaddd
          bbbbddd
          bbbb###
          "
)

ggsave(filename = '/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/fig3_main.pdf', 
       width = 210, height = 150, units = 'mm', dpi = 300,device = cairo_pdf())
