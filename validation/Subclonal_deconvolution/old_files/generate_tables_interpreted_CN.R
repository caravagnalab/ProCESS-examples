library(ggplot2)
library(tidyverse)
library(ProCESS)
library(ggpubr)
library(patchwork)
library(randnet)
library(scales)
library(ggrepel)

spn = 'SPN07'
purity=0.6
vcf_caller = "mutect2"
cna_caller = "ascat"
coverage = 150

coverage_list = c(50, 100, 150)
purity_list = c(0.3, 0.6, 0.9)
vcf_caller_list = c("mutect2")#, "strelka", "freebayes")
cna_caller_list = c("ascat")#, "sequenza", "battenberg")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN06', 'SPN07')
tool = 'pyclonevi'

combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller=cna_caller_list,
                    spn=spn_list)
# tool = 'viber'

# for(spn in spn_list){
for(i in 1:nrow(combs)){
  
  coverage = combs[i, "coverage"]
  purity = combs[i, "purity"]
  vcf_caller = combs[i, "vcf_caller"]
  cna_caller = combs[i, "cna_caller"]
  spn = combs[i, "spn"]
  
  simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
  
  # simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
  github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
  main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
  save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

  source(file.path(save_path, "generate_table_main.R"))
  source(file.path(save_path, "utils_plots.R"))
  
  # True process drivers
  true_drivers_table = readRDS(file.path(main_path,"drivers", spn, "process_drivers.rds")) %>% as_tibble()
  
  true_drivers_table = true_drivers_table  %>% rowwise() %>% 
    mutate(code=replace(code, is.na(code), paste0(type, '_', CNA_type, '_', chr, '_', start, '_', end)),
           mutation_id=paste0(spn, ":", chr, ":", start,  ":", alt)) %>% 
    ungroup()

  
  # Get process table
  mut_process = get_mutations(spn=spn, type="tumour", coverage=coverage, purity=purity)
  table_process = readRDS(get_table_path(save_path, 'process', spn, simulation_id)) # process table in folder tables/

  # Join process table with drivers
    # now in process_table we have a column "code" with the drivers gene names
  table_process = table_process %>% left_join(true_drivers_table, by = 'mutation_id') %>% 
    select(-driver_label_process) %>% 
    mutate(driver_label_process=code)


  # SPN03:9:136496197: deve diventare SPN03:9:136496196:C
  if(spn=='SPN03' & tool != 'process_viber'){
    table_process = table_process %>% mutate(mutation_id = if_else(
      mutation_id=="SPN03:9:136496197:",
      'SPN03:9:136496196:C',
      mutation_id
    ))
  }

  table_tool = readRDS(get_table_path(save_path, tool, spn, simulation_id))

  if(tool =='viber_heuristics'){
    table_tool = table_tool %>%
      mutate(cluster_id_tool = replace_na(cluster_id_tool, 'NA'))
  }
  if(tool =='process_viber'){
    table_tool = table_tool %>%
      mutate(sample_id = paste0(spn,"_",sample_id))
  }

  join_table_tool = table_tool %>% left_join(table_process) # keep all mut in tool
  join_table_final = join_table_tool %>% filter(!is.na(cluster_id_process)) # only mutations present in both

  nmi_complete = randnet::NMI(as.factor(join_table_final$cluster_id_tool),
                      as.factor(join_table_final$cluster_id_process))
  
  ari_complete = aricode::ARI(as.factor(join_table_final$cluster_id_tool), 
              as.factor(join_table_final$cluster_id_process))

  ### Find cluster/driver in tool and add column cluster_id_tool_interpreted
  # driver_clusters_tool = join_table_tool %>%
  #   distinct(cluster_id_tool, is_driver_tool) %>% 
  #   filter(is_driver_tool == TRUE) %>% 
  #   pull(cluster_id_tool)
  
  check_CN_driver = function(df, true_drivers_table){
    CNA_drivers = true_drivers_table %>% filter(type=="CNA")
    
    a = lapply(1:nrow(df), function(j) {
      mutation_id = df[j, "mutation_id"]
      minor_cn = df[j, "minor_cn"]
      major_cn = df[j, "major_cn"]
      
      splt = strsplit(mutation_id, split=":")[[1]]
      chr = splt[[2]]
      pos = splt[[3]]
      
      lapply(1:nrow(CNA_drivers), function(i) {
        if (chr != CNA_drivers[i,"chr"]) return(NA)
        if (pos > CNA_drivers[i, "start"] & pos < CNA_drivers[i, "end"]) return(CNA_drivers[i, "code"])
        if (CNA_drivers[i, "CNA_type"] == "D")
          if (minor_cn < 1 | major_cn < 1) return(CNA_drivers[i, "code"])
        return(NA)
      }) %>% unlist()
    }) %>% unlist() %>% unique() %>% purrr::discard(function(x) is.na(x))
    
    if (length(a) == 0) return(NA)
    a
    
  }
  
  join_table_tool = join_table_tool %>%
    group_by(cluster_id_tool) %>% 
    mutate(cna_driver=check_CN_driver(data.frame(mutation_id=mutation_id, minor_cn=minor_cn, major_cn=major_cn), true_drivers_table)) %>% 
    ungroup()
  
  driver_clusters_tool = join_table_tool %>% 
    filter(is_driver_tool == TRUE | !is.na(cna_driver)) %>% 
    pull(cluster_id_tool) %>% unique()
  
  final_table = join_table_tool %>%
    mutate(cluster_id_tool_interpreted = if_else(
      !(cluster_id_tool %in% driver_clusters_tool),
      'Subclonal',
      as.character(cluster_id_tool)
    ))
  
  ### Find cluster/driver in process and add column cluster_id_tool_interpreted_driver
  if(tool != 'mobster'){
  final_table_interpreted = join_table_final %>%
    mutate(cluster_id_tool_interpreted = if_else(
      !(cluster_id_tool %in% driver_clusters_tool),
      'Subclonal',
      as.character(cluster_id_tool)
    )) %>%
    mutate(cluster_id_tool_interpreted_driver = if_else(
      !(cluster_id_tool %in% driver_clusters_process),
      'Subclonal',
      as.character(cluster_id_tool)
    ))

  nmi_interpreted = randnet::NMI(as.factor(final_table_interpreted$cluster_id_tool_interpreted),
                               as.factor(final_table_interpreted$cluster_id_process))
  
  ari_interpreted = aricode::ARI(as.factor(final_table_interpreted$cluster_id_tool_interpreted),
                              as.factor(final_table_interpreted$cluster_id_process))
  }else{
    final_table_interpreted = join_table_final
  }
  # La tabella da salvare è final_table_interpreted, dove devo salvare le colonne:
  # patient_id, sample_id, coverage, purity, tool, mutation_id, driver_label_tool,
  # is_driver_tool, cluster_id_tool, vaf_tool, is_driver_process, cluster_id_process,
  # vaf_process, driver_label_process, cluster_id_tool_interpreted
  # table_to_save = final_table_interpreted %>% select(patient_id, sample_id,coverage, purity, 
  #                                                    tool, mutation_id, driver_label_tool,
  #                                                    is_driver_tool, cluster_id_tool, vaf_tool, 
  #                                                    is_driver_process, cluster_id_process,
  #                                                    vaf_process, driver_label_process, 
  #                                                    cluster_id_tool_interpreted)
  
  table_to_save = final_table_interpreted 
  saveRDS(table_to_save, file.path(main_path, "subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds")))
  saveRDS(table_to_save, file.path(save_path, "tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds")))

}



