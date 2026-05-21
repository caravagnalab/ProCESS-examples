### ASSEMBLE FINAL PANEL
base = "/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/"
source(paste0(base, 'figure3/cna/cna.R'))
source(paste0(base, 'figure3/somatic/somatic.R'))
base = "/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/"
source(paste0(base, 'figure3/segm_cna/plot_driver_cna.R'))

library(rstatix)
stat_test <- df_all_SPN_somatic %>%
  mutate(CCF_bin = factor(CCF_bin, levels = c("0-0.05","0.05-0.10","0.10-0.25",
                                              "0.25-0.50","0.50-0.75","0.75-0.95","0.95-1"))) %>%
  group_by(CCF_bin, muts) %>%                        # add muts if faceted
  wilcox_test(sensitivity ~ caller) %>%              # pairwise wilcoxon
  adjust_pvalue(method = "BH") %>%                   # FDR correction
  add_significance(p.col = "p") 
stat_test %>% View()

cartoon <- ggplot() + labs(tag = "A")
caller_line_boxes <- caller_line_boxes  + labs(tag = "B")
scatter_purity_ploidy <- scatter_purity_ploidy + labs(tag = "C")  + theme(legend.position = 'none')
alluvial_plot <- alluvial_plot + labs(tag = "D")  + theme(legend.position = 'none')
plt_breakpoint <- plt_breakpoint + labs(tag = "E") +   guides(
  col  = guide_legend(nrow = 2),
  fill = guide_legend(nrow = 3),
  shape = guide_legend(nrow = 2)
) + theme(
  legend.box = "horizontal",
  legend.key.justification = "top",
  legend.spacing.x = unit(0.1, "cm"),
  panel.grid.minor = element_blank(),
  legend.box.spacing = unit(0, "cm")
)
driver_gene_plot <- driver_gene_plot + 
  guides(
    col  = guide_legend(nrow = 2),
    fill = guide_legend(nrow = 3),
    shape = guide_legend(nrow = 2)
) +
theme(
    legend.box = "horizontal",
    legend.key.justification = "top",
    legend.spacing.x = unit(0.1, "cm"),
    panel.grid.minor = element_blank(),
    legend.box.spacing = unit(0, "cm")
  )+labs(tag = "F")

setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3')
ggsave(plot = alluvial_plot, filename = 'plt_alluvial.pdf', 
       width = 130, height = 50, units = 'mm')
ggsave(plot = scatter_purity_ploidy, filename = 'plt_scatter.pdf', 
       width = 130, height = 60, units = 'mm')
ggsave(plot = caller_line_boxes, filename = 'plt_somatic.pdf', 
       width = 130, height = 120, units = 'mm')
ggsave(plot = plt_breakpoint, filename = 'plt_breakpoint.pdf', 
       width = 130, height = 80, units = 'mm')
ggsave(plot = driver_gene_plot, filename = 'plt_driver.pdf', 
       width = 130, height = 80, units = 'mm')


# wrap_plots(
#   free(cartoon),
#   free(caller_line_boxes), 
#   free(scatter_purity_ploidy),
#   free(alluvial_plot),
#   free(plt_breakpoint),
#   free(driver_gene_plot),
#   # design="aaaccc\nbbbccc\nbbbccc\nbbbddd\nbbbddd\neeefff\neeefff\neeefff"
#   # design="aaaccc\naaaccc\nbbbccc\nbbbddd\nbbbddd\neeefff\neeefff\neeefff"
#   design = 'aaaaaa\nbbbccc\nbbbccc\nbbbddd\neeefff\neeefff\neeefff'
# )
# 
# 
# ggsave(filename = paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/", 'fig3_main_new.pdf'), 
#         width = 260, height = 230, units = 'mm', dpi = 300, device = cairo_pdf) #device = cairo_pdf()
