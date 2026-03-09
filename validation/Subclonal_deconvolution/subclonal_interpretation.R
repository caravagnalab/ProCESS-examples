library(ggplot2)
library(tidyverse)
library(ProCESS)


spn = 'SPN01'
purity=0.9
vcf_caller = "mutect2"
cna_caller = "ascat"
coverage = 100

coverage_list = c(50,100,150)
purity_list = c(0.3, 0.6, 0.9)
vcf_caller_list = c("mutect2")#, "strelka", "freebayes")
cna_caller_list = c("ascat")#, "sequenza", "battenberg")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04','SPN05', 'SPN06', 'SPN07')
tool_list = c( 'viber', 'pyclonevi')

tool = 'viber'

combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller=cna_caller_list,
                    spn=spn_list,
                    tool=tool_list)

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

source(file.path(save_path, "generate_table_main.R"))
source(file.path(save_path, "utils_tables.R"))

df_to_save = data.frame()

for(i in 1:nrow(combs)){
  
  coverage = combs[i, "coverage"]
  purity = combs[i, "purity"]
  vcf_caller = combs[i, "vcf_caller"]
  cna_caller = combs[i, "cna_caller"]
  spn = combs[i, "spn"]
  tool = combs[i, "tool"]
  
  simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
  print(paste0(tool, '_' ,spn,'_' ,simulation_id))
    
  # Extract tool table
  tool_table = tryCatch(
    readRDS(file.path(main_path, "validation_subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds"))),
    error = function(e) {
      message("Table not found: ", simulation_id,
              " (", e$message, ")")
      return(NULL)
    }
  )
  
  if (is.null(tool_table)) {
    next
  }
  
  # Extract mobster table and only select mutation_id and mobster_cluster = cluster_id_tool
   
  mobster_table = tryCatch(
    readRDS(file.path(main_path, "validation_subclonal/tables_interpreted", paste0("mobster_univariate_", spn, "_", simulation_id, ".rds"))),
    error = function(e) {
      message("Table not found: ", simulation_id,
              " (", e$message, ")")
      return(NULL)
    }
  )
  
  if (is.null(mobster_table)) {
    next
  }
  
  mobster_table = mobster_table %>% 
    dplyr::rename(mobster_cluster = cluster_id_tool) %>% 
    select(mutation_id, mobster_cluster)
  
  # Create a join with also mobster clusters
  
  final_table = tool_table %>% left_join(mobster_table)
  
  final_table = final_table %>% 
    mutate(mobster_cluster = ifelse(((vaf_tool > 0) & (vaf_tool <= 0.05) & (is.na(mobster_cluster))), 'Tail', mobster_cluster)) %>% 
    mutate(mobster_cluster = ifelse(((vaf_tool == 0) & (is.na(mobster_cluster))), 'Private', mobster_cluster))
  
  # I want to create a final table with:
    # spn, purity, coverage, cna_caller, vcf_caller, cluster_id_tool, contains_driver_tool, contains_driver_process, %of mobster tail mutations
    # I do not want sample_id because in the multivariate there is no distinction
  
  drivers_tool = final_table %>% 
    filter(is_driver_tool == TRUE) %>% 
    select(mutation_id, cluster_id_tool, cluster_id_process, 
           is_driver_tool, is_driver_process) %>% 
    distinct(mutation_id, cluster_id_tool, 
             is_driver_process, cluster_id_process) %>% 
    pull(cluster_id_tool)
  
  # Table of process drivers
  drivers_process = final_table %>% 
    filter(is_driver_process == TRUE) %>% 
    select(mutation_id, cluster_id_tool, cluster_id_process, 
           is_driver_tool, is_driver_process) %>% 
    distinct(mutation_id, cluster_id_tool, 
             is_driver_tool, cluster_id_process)%>% 
    pull(cluster_id_tool)
  
  
  # contains_tail = final_table %>%
  #   group_by(cluster_id_tool) %>%
  #   summarise(
  #     total_mutations = n(),
  #     tail_mutations = sum(mobster_cluster == "Tail", na.rm = TRUE),
  #     percent_tail = round(100 * tail_mutations / total_mutations, 3)
  #   ) %>%
  #   ungroup()
  
  contains_tail_patient = final_table %>%
    group_by(cluster_id_tool, sample_id) %>%
    summarise(
      total_mutations = n(),
      tail_mutations = sum(mobster_cluster == "Tail", na.rm = TRUE),
      percent_tail = round(100 * tail_mutations / total_mutations, 3)
    ) %>%
    ungroup() %>% 
    group_by(cluster_id_tool) %>%
    summarise(percent_tail_sample = list(percent_tail)) %>%
    ungroup()
  
  never_tail = final_table %>% 
    group_by(mutation_id) %>% 
    mutate(never_tail = all(mobster_cluster != "Tail")) %>% 
    group_by(cluster_id_tool) %>% 
    summarise(n_never_tail = sum(never_tail)/n())
  
  single_sample_tail = final_table %>%
    group_by(cluster_id_tool, sample_id) %>%
    summarise(
      total_mutations = n(),
      tail_mutations = sum(mobster_cluster == "Tail", na.rm = TRUE),
      percent_tail = round(100 * tail_mutations / total_mutations, 3)
    ) %>%
    ungroup() %>%
    group_by(cluster_id_tool) %>%
    summarise(one_sample_tail = sum(percent_tail > 80) >= 1)%>%
    ungroup()
  
  contains_tail = final_table %>%
    group_by(cluster_id_tool) %>%
    summarise(
      total_mutations = n(),
      tail_mutations = sum(mobster_cluster == "Tail", na.rm = TRUE),
      percent_tail = round(100 * tail_mutations / total_mutations, 3)
    ) %>%
    ungroup()
  
  contains_tail = contains_tail %>% inner_join(contains_tail_patient)
  contains_tail = contains_tail %>% inner_join(single_sample_tail)
  
  final_df = final_table %>%
    distinct(cluster_id_tool, purity, coverage, patient_id, is_clonal_tool) %>%
    arrange(purity, coverage, patient_id, cluster_id_tool, is_clonal_tool) 
  
  # if(tool=='viber'){
  #   final_df = final_df%>% 
  #     mutate(cluster_id_tool = as.character(cluster_id_tool))
  # }
  

  final_df = final_df %>% 
    mutate(contains_driver_tool=ifelse(cluster_id_tool %in% drivers_tool, TRUE, FALSE),
           contains_driver_process=ifelse(cluster_id_tool %in% drivers_process, TRUE, FALSE),
           tool = tool) %>% 
    inner_join(contains_tail) %>% 
    inner_join(never_tail) %>% 
    mutate(cluster_id_tool = as.character(cluster_id_tool))
  
  df_to_save = bind_rows(df_to_save, final_df)
  
}

saveRDS(df_to_save, '/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/subclonal_deconvolution_new.rds')
