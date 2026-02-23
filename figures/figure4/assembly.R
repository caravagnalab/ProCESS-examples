source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/plot_cluster_analysis.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/plot_driver_analysis.R')


p_all <- ggplot() + labs(tag = "A") + 
  free(p_cluster +  labs(tag = "B")) + 
  free(p_driver + labs(tag = "C")) + 
  free(ggplot() + labs(tag = "D")) +
  free(ggplot() + labs(tag = "E")) +
  plot_layout( design = 'AAABBBC\nDDDEEEE') 
  #theme(legend.position = 'bottom')

ggsave(filename = paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/all.pdf"), 
       plot = p_all, device="pdf", width=10, height=7, units="in")

