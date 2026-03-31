library(ggplot2)
library(tidyverse)
library(ProCESS)
library(ggpubr)
library(patchwork)
library(randnet)
library(scales)
library(ggrepel)

coverage_list = c(50,100,150)
purity_list = c(0.3, 0.6, 0.9)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04','SPN05', 'SPN06', 'SPN07')
# spn_list = c('SPN01')
tool_list = c('mobster', 'pyclonevi', 'viber')

combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller=cna_caller_list,
                    spn=spn_list)

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

source(file.path(save_path, "generate_table_main.R"))
source(file.path(save_path, "utils_tables.R"))

for(tool in tool_list){
  for(i in 1:nrow(combs)){
    
    coverage = combs[i, "coverage"]
    purity = combs[i, "purity"]
    vcf_caller = combs[i, "vcf_caller"]
    cna_caller = combs[i, "cna_caller"]
    spn = combs[i, "spn"]
    
    simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
    
    # print(paste0(spn, '_', simulation_id))
    
    table_tool = tryCatch(
      readRDS(get_table_path(paste0(main_path, "validation_subclonal_new/"), tool, spn, simulation_id)),
      error = function(e) {
        # print(paste0(spn, '_', simulation_id))
        message(paste0(tool, '_', spn, '_', simulation_id))
        return(NULL)
      }
    )
    
  }
}



