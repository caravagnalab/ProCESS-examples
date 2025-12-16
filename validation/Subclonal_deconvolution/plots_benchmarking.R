library(ggplot2)
library(tidyverse)
library(ProCESS)
library(ggpubr)
library(patchwork)
library(randnet)
library(scales)

# coverage_list = c(50, 100, 150, 200)
# purity_list = c(0.3, 0.6, 0.9)
# vcf_caller_list = c("mutect2", "strelka", "freebayes")
# cna_caller_list = c("ascat", "sequenza", "battenberg")
# spn_list = paste("SPN", 3:7, sep="0")

coverage_list = c(50,100,150)
purity_list = c(0.3, 0.6, 0.9)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN06', 'SPN07')
tool = 'viber'

combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller=cna_caller_list,
                    spn=spn_list)

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

source(file.path(save_path, "utils_plots_final.R"))

# "sample_id" "cluster_id_process" "cluster_id_tool" "ccf_process" "ccf_tool" 
table_ccf_final = data.frame()

# if (!is.na(i)) {
for(i in 1:nrow(combs)){
  coverage = combs[i, "coverage"]
  purity = combs[i, "purity"]
  vcf_caller = combs[i, "vcf_caller"]
  cna_caller = combs[i, "cna_caller"]
  spn = combs[i, "spn"]
  
  simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
  
  table = readRDS(file.path(main_path, "subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds")))
  
  table = table %>% 
      group_by(cluster_id_process, sample_id) %>%
      mutate(ccf_process=mean(ccf_process, na.rm=TRUE)) %>%
      ungroup()
  
  # Plot CCF ####
  table_ccf = table %>% 
    filter(is_driver_process==TRUE)%>% 
    select(sample_id, cluster_id_process, cluster_id_tool, ccf_process, ccf_tool, purity, coverage) 
  
  table_ccf_final = rbind(table_ccf_final, table_ccf)
  
# View(unique(table_ccf[c('ccf_process', 'sample_id', 'cluster_id_process', 'cluster_id_tool', 'ccf_tool')]))
}

saveRDS(table_ccf_final, file.path(save_path, "metrics_tables/table_ccf.rds"))

p = ggplot(data = table_ccf_final, aes(x = ccf_process, y=ccf_tool)) + 
  geom_abline()+
  geom_point()+
  # theme_minimal()+
  xlim(0,1)+
  ylim(0,1)+
  facet_grid(~purity)
  # facet_grid(coverage~purity)

p

ggsave(file.path(save_path, "plots/metrics/ccf.png"), p)






