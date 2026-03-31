source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure5/signatures_cohort/signatures.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature/plot_res.R')


p_all <- ggplot() + labs(tag = "A") + 
  free(plt_signatures +  labs(tag = "B")) + 
  free(plt_multi + labs(tag = "C")) + 
  plot_layout( design = 'AAABBBCC') & theme(legend.position = 'bottom')

ggsave(filename = paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure5/fig5_main.pdf"), 
       plot = p_all, device="pdf", width=10, height=4, units="in")

