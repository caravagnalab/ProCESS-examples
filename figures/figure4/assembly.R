source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/drivers/plot_driver.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/signatures.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/subclone_mobster.R')
#source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/js.R')

p_all <- plt_driver + labs(tag = "A") + 
  plt_signatures + ylim(0.6,1.05) + labs(tag = "B") +
  plt_subclone + labs(tag = "C") + guides(
    color = guide_legend(nrow = 2),
    fill  = guide_legend(nrow = 2)
  ) + 
  guides(
    color = guide_legend(nrow = 2),
    fill  = guide_legend(nrow = 2)
  )+
  plot_layout( design = 'ABBCC') & theme(legend.position = 'bottom')

ggsave(filename = paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/fig4_main.pdf"), 
       plot = p_all, device="pdf", width=9.5, height=4, units="in")
