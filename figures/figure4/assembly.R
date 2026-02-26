source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/plot_cluster_analysis.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/plot_driver_analysis.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/js.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/muts_rate.R')

p_all <- ggplot() + labs(tag = "A") + 
  free(p_cluster +  labs(tag = "B")) + 
  free(p_driver + labs(tag = "C")) + 
  free(plt_jaccard + labs(tag = "D") + theme(legend.position = 'right') + xlab('Cluster Type')) +
  free(plt_mut_rate + labs(tag = "E")) +
  plot_layout( design = 'AAABBBB\nCDDDEEE') 

ggsave(filename = paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/all.pdf"), 
       plot = p_all, device="pdf", width=10, height=7, units="in")

