library(ggplot2)
library(tidyverse)
library(ProCESS)
library(dplyr)
library(ggpubr)
library(patchwork)
library(randnet)
library(scales)
library(ggrepel)

spn = 'SPN07'
purity=0.3
vcf_caller = "mutect2"
cna_caller = "ascat"
coverage = 100

coverage_list = c(50,100,150)
purity_list = c(0.3, 0.6, 0.9)
vcf_caller_list = c("mutect2")#, "strelka", "freebayes")
cna_caller_list = c("ascat")#, "sequenza", "battenberg")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN06', 'SPN07')
spn_list = c('SPN07')

tool = 'viber'
univariate = F

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
    ungroup() %>% 
    select(mutation_id, code)

  
  # Get process table
  mut_process = get_mutations(spn=spn, type="tumour", coverage=coverage, purity=purity)
  if(univariate==F){
    table_process = readRDS(get_table_path(save_path, 'process_new', spn, simulation_id)) # process table in folder tables/
  }else{
    table_process = readRDS(get_table_path(save_path, 'process_univariate', spn, simulation_id)) # process table in folder tables/
  }
  
  # Join process table with drivers
    # now in process_table we have a column "code" with the drivers gene names
  table_process = table_process %>% left_join(true_drivers_table, by = 'mutation_id') %>% 
    select(-driver_label_process) %>% 
    mutate(driver_label_process=code)%>% 
    select(-code) 
  
  # SPN03:9:136496197: deve diventare SPN03:9:136496196:C
  if(spn=='SPN03' & tool != 'process_viber'){
    table_process = table_process %>% mutate(mutation_id = if_else(
      mutation_id=="SPN03:9:136496197:",
      'SPN03:9:136496196:C',
      mutation_id
    ))
  }
  
  # ----------- Create Clonal clusters process ------------- ##
  # Save the total number of sample to select clonal clusters (i.e. ccf > 0.95 in n_samples clusters)
  n_samples = length(table_process$sample_id %>% unique())
  
  c = table_process %>%
    distinct(cluster_id_process, sample_id) %>%   # keep unique pairs
    dplyr::count(cluster_id_process) %>% filter(n==n_samples) %>% 
    select(cluster_id_process)
  
  table_process = table_process %>% 
    group_by(cluster_id_process) %>% 
    mutate(
      is_clonal_process = if_else(
        dplyr::first(cluster_id_process) %in% c$cluster_id_process, # first returns the first value in the current group
        all(ccf_process > 0.95),
        FALSE
      )
    ) %>% 
    # mutate(is_clonal_process=replace(FALSE, all(ccf_process > 0.95), TRUE)) %>% ungroup() %>% 
    dplyr::mutate(cluster_id_process_full = cluster_id_process) %>% 
    dplyr::mutate(cluster_id_process=replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal'))
  
  # ----------- Create Clonal clusters process ------------- ##
  
  
  table_tool = readRDS(get_table_path(save_path, tool, spn, simulation_id))
  
  join_table_tool = table_tool %>% left_join(table_process) # keep all mut in tool
  join_table_final = join_table_tool %>% filter(!is.na(cluster_id_process)) # only mutations present in both
  join_table_tool = join_table_tool %>% filter(!is.na(cluster_id_process))
  
  
  ### Find cluster/driver in tool and add column cluster_id_tool_interpreted
  driver_clusters_tool = join_table_tool %>%
    distinct(cluster_id_tool, is_driver_tool) %>%
    filter(is_driver_tool == TRUE) %>%
    pull(cluster_id_tool)
  
  ### Find cluster/driver in process and add column cluster_id_tool_interpreted_driver
  driver_clusters_process = join_table_tool %>%
    distinct(cluster_id_tool, is_driver_process) %>%
    filter(is_driver_process == TRUE) %>%
    pull(cluster_id_tool)
  
  if(tool != 'mobster'){
    
    # ---------- Find private clusters due to drivers tool ---------- #
    # join_table_final %>% select(sample_id, ccf_tool, cluster_id_tool)
    # View(unique(join_table_final[, c("sample_id", "ccf_tool", "cluster_id_tool")]))
    no_drivers = join_table_final %>% filter(is_driver_tool == F)
    uniques = unique(no_drivers[, c("sample_id", "ccf_tool", "cluster_id_tool")])
    
    tmp = join_table_final %>% filter(is_driver_tool==T)
    # Find couples sample_id/ccf_process
    ccf_tool_drivers = unique(tmp[, c("sample_id", "ccf_tool", "cluster_id_tool")])
    
    # Private clusters without tool drivers
    private_clonal_clusters = uniques %>%
      group_by(cluster_id_tool) %>%
      summarise(n = sum(ccf_tool > 0.01, na.rm = TRUE), .groups = "drop") %>%
      filter(n == 1) %>%
      pull(cluster_id_tool)
    # non vogliamo solo questo, vogliamo anche che la ccf coincida con una di quelle nei driver clusters
    # quindi vogliamo che la coppia ccf_cluster_private/sample coincida ad una coppia ccf_cluster_driver/sample
    # quindi dopo aver filtrato n == 1 che prende i privati, guardiamo se ccf/sample esiste nei ccf_driver/sample 
    # private_clonal_clusters = 
    
    # Also save ccf of private clusters
    a = uniques %>% filter(cluster_id_tool %in% private_clonal_clusters)
    
    thr = 0.1 # threshold to choose private clonal clusters
    
    a_filt = a %>%
      # it's ok here to have a many-to-many relationship, we are matching by sample_id in order to check any possible relationship
      inner_join(
        ccf_tool_drivers  %>% 
          dplyr::mutate(driver_ccf = ccf_tool)%>%
          select(sample_id, driver_ccf),
        by = "sample_id"
      ) %>%
      filter(abs(ccf_tool - driver_ccf) <= thr) %>%
      distinct(sample_id, cluster_id_tool, ccf_tool, .keep_all = TRUE)
    
    # driver_clusters_tool = unique(c(private_clonal_clusters, driver_clusters_tool))
    driver_clusters_tool = unique(c(a_filt$cluster_id_tool, driver_clusters_tool))
    
    if(tool =='viber' & spn=='SPN07'){
      c_07 = join_table_final %>%
        distinct(cluster_id_tool, sample_id) %>%   # keep unique pairs
        dplyr::count(cluster_id_tool) %>% filter(n==n_samples) %>% 
        select(cluster_id_tool)
      
      clonal_clusters = join_table_final %>%
        group_by(cluster_id_tool) %>%
        summarise(
          is_clonal_process = if_else(
            dplyr::first(cluster_id_tool) %in% c_07$cluster_id_tool,
            all(ccf_tool > 0.9),
            FALSE
          ),
          .groups = "drop"
        ) %>%
        filter(is_clonal_process) %>%
        pull(cluster_id_tool)
      
      driver_clusters_tool = unique(c(clonal_clusters, driver_clusters_tool))
      driver_clusters_process = unique(c(clonal_clusters, driver_clusters_process))
    }
    
    # ---------- Find private clusters due to drivers tool ---------- #
    
    # ---------- Find private clusters due to drivers process (i.e. true drivers) ---------- #
    
    # join_table_final %>% select(sample_id, ccf_tool, cluster_id_tool)
    # View(unique(join_table_final[, c("sample_id", "ccf_tool", "cluster_id_tool")]))
    no_drivers_p = join_table_final %>% filter(is_driver_process == F)
    uniques_p = unique(no_drivers_p[, c("sample_id", "ccf_tool", "cluster_id_tool")])
    
    tmp_p = join_table_final %>% filter(is_driver_process==T)
    ccf_process_drivers = unique(tmp_p[, c("sample_id", "ccf_tool", "cluster_id_tool")])
    
    # private_clonal_clusters_p (i.e. private clusters without drivers) now are different because there are different drivers
    private_clonal_clusters_p = uniques_p %>%
      group_by(cluster_id_tool) %>%
      summarise(n = sum(ccf_tool > 0.01, na.rm = TRUE), .groups = "drop") %>%
      filter(n == 1) %>%
      pull(cluster_id_tool)
   
    # Also ccf of private clusters are the same as before
    a_p = uniques %>% filter(cluster_id_tool %in% private_clonal_clusters_p)
    
    thr = 0.1 # threshold to choose private clonal clusters
    
    a_filt_p = a_p %>%
      # it's ok here to have a many-to-many relationship, we are matching by sample_id in order to check any possible relationship
      inner_join(
        ccf_process_drivers  %>% 
          dplyr::mutate(driver_ccf = ccf_tool)%>%
          select(sample_id, driver_ccf),
        by = "sample_id"
      ) %>%
      filter(abs(ccf_tool - driver_ccf) <= thr) %>%
      distinct(sample_id, cluster_id_tool, ccf_tool, .keep_all = TRUE)
    
    # driver_clusters_tool = unique(c(private_clonal_clusters, driver_clusters_tool))
    driver_clusters_process = unique(c(a_filt_p$cluster_id_tool, driver_clusters_process))
    
    # ---------- Find private clusters due to drivers process (i.e. true drivers) ---------- #
    
    final_table_interpreted = join_table_final %>%
      mutate(cluster_id_tool_interpreted = if_else(
        !(cluster_id_tool %in% driver_clusters_tool),
        'Subclonal',
        as.character(cluster_id_tool)
      )) %>%
    # also add cluster_id_tool_interpreted_driver
    mutate(cluster_id_tool_interpreted_driver = if_else(
      !(cluster_id_tool %in% driver_clusters_process),
      'Subclonal',
      as.character(cluster_id_tool)
    ))
    
    nmi_interpreted = randnet::NMI(as.factor(final_table_interpreted$cluster_id_tool_interpreted),
                                 as.factor(final_table_interpreted$cluster_id_process))
    
    ari_interpreted = aricode::ARI(as.factor(final_table_interpreted$cluster_id_tool_interpreted),
                                as.factor(final_table_interpreted$cluster_id_process))
    
    nmi_interpreted_drivers = randnet::NMI(as.factor(final_table_interpreted$cluster_id_tool_interpreted_driver),
                                   as.factor(final_table_interpreted$cluster_id_process))
    
    ari_interpreted_drivers = aricode::ARI(as.factor(final_table_interpreted$cluster_id_tool_interpreted_driver),
                                   as.factor(final_table_interpreted$cluster_id_process))
  }else{
    final_table_interpreted = join_table_final
  }
  
  table_to_save = final_table_interpreted 
  if(univariate == F){
    # saveRDS(table_to_save, file.path(main_path, "subclonal/tables_interpreted_new_branches", paste0(tool, "_", spn, "_", simulation_id, ".rds")))
    saveRDS(table_to_save, file.path(save_path, "tables_interpreted_new_branches", paste0(tool, "_", spn, "_", simulation_id, ".rds")))
  }else{
    saveRDS(table_to_save, file.path(main_path, "subclonal/tables_interpreted_new_branches", paste0(tool, "_univariate_", spn, "_", simulation_id, ".rds")))
    # saveRDS(table_to_save, file.path(save_path, "tables_interpreted", paste0(tool, "_univariate_", spn, "_", simulation_id, ".rds")))
    
  }
}



