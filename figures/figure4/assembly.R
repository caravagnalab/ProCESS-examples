source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/drivers/plot_driver.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/old_figure5/signatures_cohort/signatures.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/subclone_mobster.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/old_figure4/js.R')

p_all <- 
  ggplot() + labs(tag = "A") +
  plt_signatures +  labs(tag = "B") +
  #free(plt_signatures +  labs(tag = "B")) + 
  plt_driver + labs(tag = "C") + 
  plt_subclone + labs(tag = "D") + guides(
    color = guide_legend(nrow = 2),
    fill  = guide_legend(nrow = 2)
  )+
  plt_jaccard + labs(tag = "E") +  theme(axis.text.x = element_blank(),
                                         axis.ticks.x = element_blank()) +
  guides(
    color = guide_legend(nrow = 2),
    fill  = guide_legend(nrow = 2)
  )+
  plot_layout( design = 'AAAABBB\nAAAABBB\nCCDDDDE\nCCDDDDE') & theme(legend.position = 'bottom')

ggsave(filename = paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/fig4_main.pdf"), 
       plot = p_all, device="pdf", width=8.5, height=7, units="in")
