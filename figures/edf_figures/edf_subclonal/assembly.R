source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/plot_nmi.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/plot_cluster_analysis.R')

p_all <- free(p + labs(tag = "A")+ theme(legend.position = 'bottom')) + 
  free(p_cluster_purity +  labs(tag = "B") + theme(legend.position = 'bottom')) + 
  free(nmi_plot + labs(tag = "C")) + 
  plot_layout( design = 'AABB\nCCCC') 

ggsave(filename = 
         paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_subclonal/edf_subclonal.pdf"), 
       plot = p_all, device="pdf", width=8, height=6, units="in")
