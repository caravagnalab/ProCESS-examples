source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/drivers/plot_driver.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/signatures.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/subclone_mobster.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/multivariate.R')

p_all <- free(plt_driver) + 
  labs(tag = "A") + 
  ggplot() + 
  labs(tag = "B") + 
  free(plt_signatures + 
  ylim(0.6,1.05) + 
  theme(legend.position = 'bottom')) + 
  labs(tag = "C") +
  ggplot() + 
  labs(tag = "D") + 
free(plt_subclone + 
  theme(legend.position = 'bottom') + 
  labs(tag = "E") + 
  guides(
    color = guide_legend(nrow = 2),
    fill  = guide_legend(nrow = 2)
  ) + 
  guides(
    color = guide_legend(nrow = 2),
    fill  = guide_legend(nrow = 2)
  )) +
free(plt_multivariate + 
  theme(legend.position = 'bottom')) + 
  labs(tag = "F") + 
  plot_layout( design = 'AAABBBCCCC\nDDDDEEEFFF') 

ggsave(filename = paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/fig4_main.pdf"), 
       plot = p_all, device="pdf", width=9.5, height=7, units="in")
