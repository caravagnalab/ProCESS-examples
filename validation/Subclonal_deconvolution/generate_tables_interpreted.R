library(ggplot2)
library(tidyverse)
library(ProCESS)
library(ggpubr)
library(patchwork)
library(randnet)
library(scales)
library(ggrepel)

# DA AGGIUNGERE MOBSTER:
  # - se un cluster finisce nella tail, mettilo come "subclonal"/"Other" ???
  # - fare interpreted + interpreted_drivers clusters mobster
  # interpreted: if a subclone has no driver --> "Other"
spn = 'SPN03'
purity=0.9
vcf_caller = "mutect2"
cna_caller = "ascat"
coverage = 100

coverage_list = c(50,100,150)
purity_list = c(0.3, 0.6, 0.9)
vcf_caller_list = c("mutect2")#, "strelka", "freebayes")
cna_caller_list = c("ascat")#, "sequenza", "battenberg")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN06', 'SPN07')
spn_list = c('SPN05')

tool = 'viber'
if(tool == 'mobster'){
  univariate = T
}else{
  univariate = F
}
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

# for(spn in spn_list){
for(i in 1:nrow(combs)){
  
  coverage = combs[i, "coverage"]
  purity = combs[i, "purity"]
  vcf_caller = combs[i, "vcf_caller"]
  cna_caller = combs[i, "cna_caller"]
  spn = combs[i, "spn"]

  simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)

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
    table_process = readRDS(get_table_path(paste0(main_path, "validation_subclonal/"), 'process', spn, simulation_id)) # process table in folder tables/
  }else{
    table_process = readRDS(get_table_path(paste0(main_path, "validation_subclonal/"), 'process_univariate_w_private', spn, simulation_id)) # process table in folder tables/
    }
  
  # Join process table with drivers
    # now in process_table we have a column "code" with the drivers gene names
  table_process = table_process %>% left_join(true_drivers_table, by = 'mutation_id') %>% 
    select(-driver_label_process) %>% 
    mutate(driver_label_process=code)


  # SPN03:9:136496197: deve diventare SPN03:9:136496196:C
  if(spn=='SPN03'){
    table_process = table_process %>% mutate(mutation_id = if_else(
      mutation_id=="SPN03:9:136496197:",
      'SPN03:9:136496196:C',
      mutation_id
    ))
  }
  
  table_tool = tryCatch(
    readRDS(get_table_path(paste0(main_path, "validation_subclonal/"), tool, spn, simulation_id)),
    error = function(e) {
      message("Skipping simulation_id: ", simulation_id,
              " (", e$message, ")")
      return(NULL)
    }
  )
  
  if (is.null(table_tool)) {
    next
  }
  # table_tool = readRDS(get_table_path(paste0(main_path, "validation_subclonal/"), tool, spn, simulation_id))

  ### Find cluster/driver in process and add column cluster_id_tool_interpreted_driver
  if(tool != 'mobster'){
    
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
      dplyr::mutate(cluster_id_process=replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal')) %>% 
      dplyr::mutate(cluster_id_process=ifelse(cluster_id_process=='Truncal', 'Clonal', cluster_id_process))
    
    # ----------- Create Clonal clusters process ------------- ##
    
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
    
    final_table_interpreted = final_table_interpreted %>%
      group_by(cluster_id_tool) %>%
      mutate(is_clonal_tool = all(ccf_tool > 0.9)) %>%
      ungroup()
    
  }else{
    join_table_tool = table_tool %>% left_join(table_process) # keep all mut in tool
    join_table_final = join_table_tool %>% filter(!is.na(cluster_id_process)) # only mutations present in both
    join_table_tool = join_table_tool %>% filter(!is.na(cluster_id_process))
    
    table = join_table_final%>% 
      group_by(cluster_id_process, sample_id) %>%
      mutate(ccf_process=mean(ccf_process, na.rm=TRUE)) %>%
      ungroup()
    
    final_table = table %>% group_by(cluster_id_process, sample_id) %>% 
      mutate(is_clonal_process=replace(FALSE, ccf_process > 0.9, TRUE)) %>% ungroup() %>% 
      mutate(cluster_id_process_full = cluster_id_process) %>% 
      mutate(cluster_id_process=replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal'))
    
    # ---- Adjust is_clonal_tool column --- #
    # final_table = table %>% group_by(cluster_id_tool, sample_id) %>% 
    #   mutate(is_clonal_tool=replace(FALSE, ccf_tool > 0.9, TRUE)) %>% ungroup() 
    # ---- Adjust is_clonal_tool column --- #
    
    sample_names = final_table$sample_id %>% unique()
    
    final_table_interpreted = data.frame()
    
    for(sample in sample_names){
      table_sample = final_table %>% 
        filter(sample_id == sample)
      
      ### Find cluster/driver in tool and add column cluster_id_tool_interpreted
      driver_clusters_tool = table_sample %>%
        distinct(cluster_id_tool, is_driver_tool) %>%
        filter(is_driver_tool == TRUE) %>%
        pull(cluster_id_tool)
      
      ### Find cluster/driver in process and add column cluster_id_tool_interpreted_driver
      driver_clusters_process = table_sample %>%
        distinct(cluster_id_tool, is_driver_process) %>%
        filter(is_driver_process == TRUE) %>%
        pull(cluster_id_tool)
      
      table_sample_final = table_sample %>%
        mutate(cluster_id_tool_interpreted = if_else(
          (!(cluster_id_tool %in% driver_clusters_tool)) & (cluster_id_tool!='Tail'),
          'Other',
          as.character(cluster_id_tool)
        )) %>%
        # also add cluster_id_tool_interpreted_driver
        mutate(cluster_id_tool_interpreted_driver = if_else(
          (!(cluster_id_tool %in% driver_clusters_process))& (cluster_id_tool!='Tail'),
          'Other',
          as.character(cluster_id_tool)
        ))
      
      final_table_interpreted = bind_rows(final_table_interpreted, table_sample_final)
    }
    # final_table_interpreted = join_table_final
  }
  
  table_to_save = final_table_interpreted 
  if(univariate == F){
    saveRDS(table_to_save, file.path(main_path, "validation_subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds")))
    # saveRDS(table_to_save, file.path(save_path, "tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds")))
  }else{
    saveRDS(table_to_save, file.path(main_path, "validation_subclonal/tables_interpreted", paste0(tool, "_univariate_", spn, "_", simulation_id, ".rds")))
    # saveRDS(table_to_save, file.path(save_path, "/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds")))
    
  }
}



