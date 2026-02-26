library(ggplot2)
library(tidyr)
library(dplyr)
library(ProCESS)
library(stringr)

spn = 'SPN02'
purity=0.6
vcf_caller = "mutect2"
cna_caller = "ascat"
coverage = 150

coverage_list = c(50,100,150)
purity_list = c(0.3, 0.6, 0.9)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN05', 'SPN06', 'SPN07')
# spn_list = c('SPN01')

combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller=cna_caller_list,
                    spn=spn_list)


tool = 'viber'
tools = c('viber', 'pyclonevi', 'mobster')

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

metrics_drivers = data.frame()

for(tool in tools){
  print(tool)
  for(i in 1:nrow(combs)){
    
    coverage = combs[i, "coverage"]
    purity = combs[i, "purity"]
    vcf_caller = combs[i, "vcf_caller"]
    cna_caller = combs[i, "cna_caller"]
    spn = combs[i, "spn"]
    
    simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
    
    print(paste0(spn,'_', simulation_id))
    
    
    if(tool != 'mobster'){
      # Multivariate ####
      # Get interpreted table
      final_table = tryCatch(
        readRDS(file.path(main_path, "validation_subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds"))),
        error = function(e) {
          message("Skipping simulation_id: ", simulation_id,
                  " (", e$message, ")")
          return(NULL)
        }
      )
      
      if (is.null(final_table)) {
        next
      }
      
      # final_table = readRDS(file.path(main_path, "validation_subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds"))) # process table in folder tables/
      n_cluster_blind = length(unique(final_table$cluster_id_tool))
      n_cluster_interpreted = length(setdiff(unique(final_table$cluster_id_tool_interpreted), 'Subclonal'))
      n_cluster_interpreted_driver = length(setdiff(unique(final_table$cluster_id_tool_interpreted_driver), 'Subclonal'))
      n_cluster_process = length(setdiff(unique(final_table$cluster_id_process), 'Subclonal'))
      
      if(nrow(final_table%>%
              pivot_wider(
                names_from = "sample_id",
                values_from = "vaf_process"
              )) < 5000){
        next
      }
      
      
      # Table of tool drivers
      drivers_tool = final_table %>% 
        filter(is_driver_tool == TRUE) %>% 
        select(mutation_id, cluster_id_tool, cluster_id_process, 
               is_driver_tool, is_driver_process) %>% 
        distinct(mutation_id, cluster_id_tool, 
                 is_driver_process, cluster_id_process)
      
      # Table of process drivers
      drivers_process = final_table %>% 
        filter(is_driver_process == TRUE) %>% 
        select(mutation_id, cluster_id_tool, cluster_id_process, 
               is_driver_tool, is_driver_process) %>% 
        distinct(mutation_id, cluster_id_tool, 
                 is_driver_tool, cluster_id_process)
      
      ## Clusters analysis ####
      
      ### Interpreted ####
      #### TP #### 
      # clusters which are both in cluster_id_tool_interpreted and in cluster_id_tool_interpreted_drivers
      c_interpreted_drivers = setdiff(unique(final_table$cluster_id_tool_interpreted_driver), "Subclonal")
      
      which_TP_c = intersect(c_interpreted_drivers, unique(final_table$cluster_id_tool_interpreted))
      TP_c = length(which_TP_c)
      
      TP_table = final_table %>% filter(cluster_id_tool %in% which_TP_c)
      
      nmi_interpreted = randnet::NMI(as.factor(TP_table$cluster_id_tool_interpreted),
                                     as.factor(TP_table$cluster_id_process))
      
      if ((length(unique(TP_table$cluster_id_tool_interpreted)) == 1 | 
           length(unique(TP_table$cluster_id_process)) == 1) & 
          (is.na(nmi_interpreted) | (nmi_interpreted == 0))) {
        
        nmi_interpreted = 1
      }
      #### FP ####
      # clusters which are in cluster_id_tool_interpreted but not in cluster_id_tool_interpreted_drivers
      which_FP_c = setdiff(unique(final_table$cluster_id_tool_interpreted), c_interpreted_drivers)
      which_FP_c = setdiff(which_FP_c, "Subclonal")
      FP_c = length(which_FP_c)
      
      #### FN ####
      # setdiff(x, y) restituisce gli elementi presenti nel vettore x ma non nel vettore y
      which_FN_c = setdiff(c_interpreted_drivers, unique(drivers_tool$cluster_id_tool)) 
      FN_c = length(which_FN_c)
      
      ### Blind ####
      #### TP #### 
      # clusters which are both in cluster_id_tool_interpreted and in cluster_id_tool_interpreted_drivers
      c_interpreted_drivers = setdiff(unique(final_table$cluster_id_tool_interpreted_driver), "Subclonal")
      which_TP_c_blind = intersect(c_interpreted_drivers, unique(final_table$cluster_id_tool))
      TP_c_blind = length(which_TP_c_blind)
      
      TP_table = final_table %>% filter(cluster_id_tool %in% which_TP_c_blind)
      
      nmi_blind = randnet::NMI(as.factor(TP_table$cluster_id_tool),
                               as.factor(TP_table$cluster_id_process))
      
      if ((length(unique(TP_table$cluster_id_tool)) == 1 | 
           length(unique(TP_table$cluster_id_process)) == 1) & 
          (is.na(nmi_blind) | (nmi_blind == 0))) {
        
        nmi_blind = 1
      }
      
     
      #### FP ####
      # clusters which are in cluster_id_tool_interpreted but not in cluster_id_tool_interpreted_drivers
      which_FP_c_blind = setdiff(unique(final_table$cluster_id_tool), c_interpreted_drivers)
      FP_c_blind = length(which_FP_c_blind)
      
      #### FN ####
      FN_c_blind = 0
      
      ## Drivers analysis ####
        # Only consider TP clusters
     
      # TP_table = final_table %>% 
      #   dplyr::filter(cluster_id_tool_interpreted %in% which_TP_c)
      
      TP_table = final_table 
      
      # Only select clonal clusters for TP clusters
      clonal_TP_table = final_table %>% 
        dplyr::filter(is_clonal_tool==TRUE) %>% 
        dplyr::filter(cluster_id_tool_interpreted %in% which_TP_c,
                      is_clonal_tool==TRUE)
      
      # Only select subclonal clusters for TP clusters
      subclonal_TP_table = final_table %>% 
        dplyr::filter(!(is_clonal_tool==TRUE)) %>% 
        dplyr::filter(cluster_id_tool_interpreted %in% which_TP_c,
                      !(is_clonal_tool==TRUE))
      
      # Only select clonal clusters among all clusters
      clonal_table = final_table %>% 
        dplyr::filter(is_clonal_tool==TRUE)
      
      # Only select subclonal clusters among all clusters
      subclonal_table = final_table %>% 
        dplyr::filter(!(is_clonal_tool==TRUE))
      
      # Only TP 
      # Retrieve tool clonal drivers only for TP clusters
      clonal_drivers_TP_tool = clonal_TP_table %>% 
        filter(is_driver_tool == TRUE) %>% 
        select(mutation_id, cluster_id_tool, cluster_id_process, 
               is_driver_tool, is_driver_process) %>% 
        distinct(mutation_id, cluster_id_tool, 
                 is_driver_process, cluster_id_process)
      
      # Retrieve process clonal drivers only for TP clusters
      clonal_drivers_TP_process = clonal_TP_table %>% 
        filter(is_driver_process == TRUE) %>% 
        select(mutation_id, cluster_id_tool, cluster_id_process, 
               is_driver_tool, is_driver_process) %>% 
        distinct(mutation_id, cluster_id_tool, 
                 is_driver_tool, cluster_id_process)
      
      # Retrieve tool subclonal drivers only for TP clusters
      subclonal_drivers_TP_tool = subclonal_TP_table %>% 
        filter(is_driver_tool == TRUE) %>% 
        select(mutation_id, cluster_id_tool, cluster_id_process, 
               is_driver_tool, is_driver_process) %>% 
        distinct(mutation_id, cluster_id_tool, 
                 is_driver_process, cluster_id_process)
      
      # Retrieve process subclonal drivers only for TP clusters
      subclonal_drivers_TP_process = subclonal_TP_table %>% 
        filter(is_driver_process == TRUE) %>% 
        select(mutation_id, cluster_id_tool, cluster_id_process, 
               is_driver_tool, is_driver_process) %>% 
        distinct(mutation_id, cluster_id_tool, 
                 is_driver_tool, cluster_id_process)
      
      
      # All clusters
      # Retrieve tool clonal drivers for all clusters
      clonal_drivers_tool = clonal_table %>% 
        filter(is_driver_tool == TRUE) %>% 
        select(mutation_id, cluster_id_tool, cluster_id_process, 
               is_driver_tool, is_driver_process) %>% 
        distinct(mutation_id, cluster_id_tool, 
                 is_driver_process, cluster_id_process)
      
      # Retrieve process clonal drivers for all clusters
      clonal_drivers_process = clonal_table %>% 
        filter(is_driver_process == TRUE) %>% 
        select(mutation_id, cluster_id_tool, cluster_id_process, 
               is_driver_tool, is_driver_process) %>% 
        distinct(mutation_id, cluster_id_tool, 
                 is_driver_tool, cluster_id_process)
      
      # Retrieve tool subclonal drivers for all clusters
      subclonal_drivers_tool = subclonal_table %>% 
        filter(is_driver_tool == TRUE) %>% 
        select(mutation_id, cluster_id_tool, cluster_id_process, 
               is_driver_tool, is_driver_process) %>% 
        distinct(mutation_id, cluster_id_tool, 
                 is_driver_process, cluster_id_process)
      
      # Retrieve process subclonal drivers for all clusters
      subclonal_drivers_process = subclonal_table %>% 
        filter(is_driver_process == TRUE) %>% 
        select(mutation_id, cluster_id_tool, cluster_id_process, 
               is_driver_tool, is_driver_process) %>% 
        distinct(mutation_id, cluster_id_tool, 
                 is_driver_tool, cluster_id_process)
      
      ### Clonal only TP ####
      #### TP ####
      # TP: how many true drivers I have in the TP clusters
      
      clonal_TP_d = clonal_drivers_TP_tool %>%
        summarise(count = sum(is_driver_process, na.rm = TRUE))
      
      clonal_which_TP_d = clonal_drivers_TP_tool %>%
        filter(is_driver_process) %>%
        pull(cluster_id_tool)
      
      #### FP ####
      # FP: how many wrong drivers I have in the TP clusters (i.e. in the interpreted clusters)
        # tool clusters where all rows have is_driver_process == FALSE (i.e., zero TRUE in that cluster).
      clonal_FP_d = clonal_drivers_TP_tool %>%
        summarise(count = sum(!is_driver_process, na.rm = TRUE))
      
      clonal_which_FP_d = clonal_drivers_TP_tool %>%
        filter(!is_driver_process) %>%
        pull(cluster_id_tool)
      
      #### FN ####
      clonal_FN_d = clonal_drivers_TP_process %>%
        summarise(count = sum(!is_driver_tool, na.rm = TRUE))
      
      clonal_which_FN_d = clonal_drivers_TP_process %>%
        filter(!is_driver_tool) %>%
        pull(cluster_id_tool)
      
      ### Sublonal only TP ####
      #### TP ####
      # TP: how many true drivers I have in the TP clusters
      
      subclonal_TP_d = subclonal_drivers_TP_tool %>%
        summarise(count = sum(is_driver_process, na.rm = TRUE))
      
      subclonal_which_TP_d = subclonal_drivers_TP_tool %>%
        filter(is_driver_process) %>%
        pull(cluster_id_tool)
      
      #### FP ####
      # FP: how many wrong drivers I have in the TP clusters (i.e. in the interpreted clusters)
      # tool clusters where all rows have is_driver_process == FALSE (i.e., zero TRUE in that cluster).
      subclonal_FP_d = subclonal_drivers_TP_tool %>%
        summarise(count = sum(!is_driver_process, na.rm = TRUE))
      
      subclonal_which_FP_d = subclonal_drivers_TP_tool %>%
        filter(!is_driver_process) %>%
        pull(cluster_id_tool)
      
      #### FN ####
      subclonal_FN_d = subclonal_drivers_TP_process %>%
        summarise(count = sum(!is_driver_tool, na.rm = TRUE))
      
      subclonal_which_FN_d = subclonal_drivers_TP_process %>%
        filter(!is_driver_tool) %>%
        pull(cluster_id_tool)
      
      
      ### Clonal all ####
      #### TP ####
      # TP: how many true drivers I have in the TP clusters
      
      all_clonal_TP_d = clonal_drivers_tool %>%
        summarise(count = sum(is_driver_process, na.rm = TRUE))
      
      all_clonal_which_TP_d = clonal_drivers_tool %>%
        filter(is_driver_process) %>%
        pull(cluster_id_tool)
      
      #### FP ####
      # FP: how many wrong drivers I have in the TP clusters (i.e. in the interpreted clusters)
      # tool clusters where all rows have is_driver_process == FALSE (i.e., zero TRUE in that cluster).
      all_clonal_FP_d = clonal_drivers_tool %>%
        summarise(count = sum(!is_driver_process, na.rm = TRUE))
      
      all_clonal_which_FP_d = clonal_drivers_tool %>%
        filter(!is_driver_process) %>%
        pull(cluster_id_tool)
      
      #### FN ####
      all_clonal_FN_d = clonal_drivers_process %>%
        summarise(count = sum(!is_driver_tool, na.rm = TRUE))
      
      all_clonal_which_FN_d = clonal_drivers_process %>%
        filter(!is_driver_tool) %>%
        pull(cluster_id_tool)
      
      ### Sublonal all ####
      #### TP ####
      # TP: how many true drivers I have in the TP clusters
      
      all_subclonal_TP_d = subclonal_drivers_tool %>%
        summarise(count = sum(is_driver_process, na.rm = TRUE))
      
      all_subclonal_which_TP_d = subclonal_drivers_tool %>%
        filter(is_driver_process) %>%
        pull(cluster_id_tool)
      
      #### FP ####
      # FP: how many wrong drivers I have in the TP clusters (i.e. in the interpreted clusters)
      # tool clusters where all rows have is_driver_process == FALSE (i.e., zero TRUE in that cluster).
      all_subclonal_FP_d = subclonal_drivers_tool %>%
        summarise(count = sum(!is_driver_process, na.rm = TRUE))
      
      all_subclonal_which_FP_d = subclonal_drivers_tool %>%
        filter(!is_driver_process) %>%
        pull(cluster_id_tool)
      
      #### FN ####
      all_subclonal_FN_d = subclonal_drivers_process %>%
        summarise(count = sum(!is_driver_tool, na.rm = TRUE))
      
      all_subclonal_which_FN_d = subclonal_drivers_process %>%
        filter(!is_driver_tool) %>%
        pull(cluster_id_tool)
      
      if(is.na(nmi_interpreted)){
        nmi_interpreted = 1
      }
      
      if(is.na(nmi_blind)){
        nmi_blind = 1
      }
      
      df <- data.frame(
        spn = spn, purity = purity, coverage = coverage,
        vcf_caller = vcf_caller, cna_caller = cna_caller,
        tool = tool,
        
        TP_c_blind = TP_c_blind,
        TP_c_blind_list = I(list(unname(which_TP_c_blind))),
        FP_c_blind = FP_c_blind,
        FP_c_blind_list = I(list(unname(which_FP_c_blind))),
        FN_c_blind = FN_c_blind,
        FN_c_blind_list = NA,
        
        TP_c = TP_c,
        TP_c_list = I(list(unname(which_TP_c))),
        FP_c = FP_c,
        FP_c_list = I(list(unname(which_FP_c))),
        FN_c = FN_c,
        FN_c_list = I(list(unname(which_FN_c))),
        
        clonal_TP_d = clonal_TP_d$count,
        clonal_TP_d_list = I(list(unname(clonal_which_TP_d))),
        clonal_FP_d = clonal_FP_d$count,
        clonal_FP_d_list = I(list(unname(clonal_which_FP_d))),
        clonal_FN_d = clonal_FN_d$count,
        clonal_FN_d_list = I(list(unname(clonal_which_FN_d))),
        
        subclonal_TP_d = subclonal_TP_d$count,
        subclonal_TP_d_list = I(list(unname(subclonal_which_TP_d))),
        subclonal_FP_d = subclonal_FP_d$count,
        subclonal_FP_d_list = I(list(unname(subclonal_which_FP_d))),
        subclonal_FN_d = subclonal_FN_d$count,
        subclonal_FN_d_list = I(list(unname(subclonal_which_FN_d))),
        
        all_clonal_TP_d = all_clonal_TP_d$count,
        all_clonal_TP_d_list = I(list(unname(all_clonal_which_TP_d))),
        all_clonal_FP_d = all_clonal_FP_d$count,
        all_clonal_FP_d_list = I(list(unname(all_clonal_which_FP_d))),
        all_clonal_FN_d = all_clonal_FN_d$count,
        all_clonal_FN_d_list = I(list(unname(all_clonal_which_FN_d))),
        
        all_subclonal_TP_d = all_subclonal_TP_d$count,
        all_subclonal_TP_d_list = I(list(unname(all_subclonal_which_TP_d))),
        all_subclonal_FP_d = all_subclonal_FP_d$count,
        all_subclonal_FP_d_list = I(list(unname(all_subclonal_which_FP_d))),
        all_subclonal_FN_d = all_subclonal_FN_d$count,
        all_subclonal_FN_d_list = I(list(unname(all_subclonal_which_FN_d))),
        
        nmi_blind = nmi_blind,
        nmi_interpreted = nmi_interpreted,
        
        n_cluster_blind = n_cluster_blind,
        n_cluster_interpreted = n_cluster_interpreted,
        n_cluster_interpreted_driver = n_cluster_interpreted_driver,
        n_cluster_process = n_cluster_process,
        
        stringsAsFactors = FALSE
      )
      
    }else{
      # Univariate ####
      # Get interpreted table
      
      final_table = tryCatch(
        readRDS(file.path(main_path, "validation_subclonal/tables_interpreted", paste0(tool, "_univariate_", spn, "_", simulation_id, ".rds"))),
        error = function(e) {
          message("Skipping simulation_id: ", simulation_id,
                  " (", e$message, ")")
          return(NULL)
        }
      )
      
      if (is.null(final_table)) {
        next
      }
      # final_table = readRDS(file.path(main_path, "validation_subclonal/tables_interpreted", paste0(tool, "_univariate_", spn, "_", simulation_id, ".rds"))) # process table in folder tables/
      
      if(nrow(final_table%>%
              pivot_wider(
                names_from = "sample_id",
                values_from = "vaf_process"
              )) < 5000){
        next
      }
      
      df = data.frame()
      samples = final_table$sample_id %>% unique()
      # sample = samples[[1]]
      
      for(sample in samples){
        
        table_sample = final_table %>% filter(sample_id==sample)
        
        n_cluster_blind = length(unique(table_sample$cluster_id_tool))
        n_cluster_interpreted = length(setdiff(unique(table_sample$cluster_id_tool_interpreted), c("Other", "Tail")))
        n_cluster_interpreted_driver = length(setdiff(unique(table_sample$cluster_id_tool_interpreted_driver), c("Other", "Tail")))
        n_cluster_process = length(setdiff(unique(table_sample$cluster_id_process), 'Subclonal'))
        
        
      # table_sample = final_table
        # Table of tool drivers
        drivers_tool = table_sample %>% 
          filter(is_driver_tool == TRUE) %>% 
          select(mutation_id, cluster_id_tool, cluster_id_process, 
                 is_driver_tool, is_driver_process) %>% 
          distinct(mutation_id, cluster_id_tool, 
                   is_driver_process, cluster_id_process)
        
        # Table of process drivers
        drivers_process = table_sample %>% 
          filter(is_driver_process == TRUE) %>% 
          select(mutation_id, cluster_id_tool, cluster_id_process, 
                 is_driver_tool, is_driver_process) %>% 
          distinct(mutation_id, cluster_id_tool, 
                   is_driver_tool, cluster_id_process)
        
        ## Clusters analysis ####
        
        ### Interpreted ####
        #### TP #### 
        c_interpreted_drivers = setdiff(unique(table_sample$cluster_id_tool_interpreted_driver), c("Other", "Tail"))
        
        which_TP_c = intersect(c_interpreted_drivers, unique(table_sample$cluster_id_tool_interpreted))
        TP_c = length(which_TP_c)
        
        TP_table = table_sample %>% filter(cluster_id_tool %in% which_TP_c)
        
        nmi_interpreted = randnet::NMI(as.factor(TP_table$cluster_id_tool_interpreted),
                                       as.factor(TP_table$cluster_id_process))
        
        
        if ((length(unique(TP_table$cluster_id_tool_interpreted)) == 1 | 
             length(unique(TP_table$cluster_id_process)) == 1) & 
            (is.na(nmi_interpreted) | (nmi_interpreted == 0))) {
          
          nmi_interpreted = 1
        }
        
        
        #### FP ####
        # clusters which are in cluster_id_tool_interpreted but not in cluster_id_tool_interpreted_drivers
        which_FP_c = setdiff(unique(table_sample$cluster_id_tool_interpreted), c_interpreted_drivers)
        which_FP_c = setdiff(which_FP_c, c("Other", "Tail"))
        FP_c = length(which_FP_c)
        
        #### FN ####
        # setdiff(x, y) restituisce gli elementi presenti nel vettore x ma non nel vettore y
        # which_FN_c = setdiff(c_interpreted_drivers, unique(drivers_tool$cluster_id_tool)) 
        # FN_c = length(which_FN_c)
        
        # Take drivers_process and check if some if "cluster_id_tool" is 'Tail'
        FN_c = drivers_process %>%
          dplyr::summarise(n_tail = sum(cluster_id_tool == "Tail")) %>% 
          pull(n_tail)
        which_FN_c = c('Tail')
        
        ### Blind ####
        #### TP #### 
        # clusters which are both in cluster_id_tool_interpreted and in cluster_id_tool_interpreted_drivers
        c_interpreted_drivers = setdiff(unique(table_sample$cluster_id_tool_interpreted_driver), c("Other", "Tail"))
        which_TP_c_blind = intersect(c_interpreted_drivers, unique(table_sample$cluster_id_tool))
        TP_c_blind = length(which_TP_c_blind)
        
        TP_table = table_sample %>% filter(cluster_id_tool %in% TP_c_blind)
        
        nmi_blind = randnet::NMI(as.factor(TP_table$cluster_id_tool),
                                 as.factor(TP_table$cluster_id_process))
        
        if ((length(unique(TP_table$cluster_id_tool)) == 1 | 
             length(unique(TP_table$cluster_id_process)) == 1) & 
            (is.na(nmi_blind) | (nmi_blind == 0))) {
          
          nmi_blind = 1
        }
        #### FP ####
        # clusters which are in cluster_id_tool_interpreted but not in cluster_id_tool_interpreted_drivers
        no_tail = setdiff(unique(table_sample$cluster_id_tool), 'Tail')
        which_FP_c_blind = setdiff(no_tail, c_interpreted_drivers)
        FP_c_blind = length(which_FP_c_blind)
        
        #### FN ####
        FN_c_blind = 0
        
        
        
        ## Drivers analysis ####
        # Only consider TP clusters
        
        # TP_table = table_sample 
        
        # Only select clonal clusters for TP clusters
        clonal_TP_table = table_sample %>% 
          dplyr::filter(is_clonal_tool==TRUE) %>% 
          dplyr::filter(cluster_id_tool_interpreted %in% which_TP_c,
                        is_clonal_tool==TRUE)
        
        # Only select subclonal clusters for TP clusters
        subclonal_TP_table = table_sample %>% 
          dplyr::filter(!(is_clonal_tool==TRUE)) %>% 
          dplyr::filter(cluster_id_tool_interpreted %in% which_TP_c,
                        !(is_clonal_tool==TRUE))
        
        # Only select clonal clusters among all clusters
        clonal_table = table_sample %>% 
          dplyr::filter(is_clonal_tool==TRUE)
        
        # Only select subclonal clusters among all clusters
        subclonal_table = table_sample %>% 
          dplyr::filter(!(is_clonal_tool==TRUE)) %>% 
          filter(cluster_id_tool != 'Tail')
        
        # Only TP 
        # Retrieve tool clonal drivers only for TP clusters
        clonal_drivers_TP_tool = clonal_TP_table %>% 
          filter(is_driver_tool == TRUE) %>% 
          select(mutation_id, cluster_id_tool, cluster_id_process, 
                 is_driver_tool, is_driver_process) %>% 
          distinct(mutation_id, cluster_id_tool, 
                   is_driver_process, cluster_id_process)
        
        # Retrieve process clonal drivers only for TP clusters
        clonal_drivers_TP_process = clonal_TP_table %>% 
          filter(is_driver_process == TRUE) %>% 
          select(mutation_id, cluster_id_tool, cluster_id_process, 
                 is_driver_tool, is_driver_process) %>% 
          distinct(mutation_id, cluster_id_tool, 
                   is_driver_tool, cluster_id_process)
        
        # Retrieve tool subclonal drivers only for TP clusters
        subclonal_drivers_TP_tool = subclonal_TP_table %>% 
          filter(is_driver_tool == TRUE) %>% 
          select(mutation_id, cluster_id_tool, cluster_id_process, 
                 is_driver_tool, is_driver_process) %>% 
          distinct(mutation_id, cluster_id_tool, 
                   is_driver_process, cluster_id_process)
        
        # Retrieve process subclonal drivers only for TP clusters
        subclonal_drivers_TP_process = subclonal_TP_table %>% 
          filter(is_driver_process == TRUE) %>% 
          select(mutation_id, cluster_id_tool, cluster_id_process, 
                 is_driver_tool, is_driver_process) %>% 
          distinct(mutation_id, cluster_id_tool, 
                   is_driver_tool, cluster_id_process)
        
        
        # All clusters
        # Retrieve tool clonal drivers for all clusters
        clonal_drivers_tool = clonal_table %>% 
          filter(is_driver_tool == TRUE) %>% 
          select(mutation_id, cluster_id_tool, cluster_id_process, 
                 is_driver_tool, is_driver_process) %>% 
          distinct(mutation_id, cluster_id_tool, 
                   is_driver_process, cluster_id_process)
        
        # Retrieve process clonal drivers for all clusters
        clonal_drivers_process = clonal_table %>% 
          filter(is_driver_process == TRUE) %>% 
          select(mutation_id, cluster_id_tool, cluster_id_process, 
                 is_driver_tool, is_driver_process) %>% 
          distinct(mutation_id, cluster_id_tool, 
                   is_driver_tool, cluster_id_process)
        
        # Retrieve tool subclonal drivers for all clusters
        subclonal_drivers_tool = subclonal_table %>% 
          filter(is_driver_tool == TRUE) %>% 
          select(mutation_id, cluster_id_tool, cluster_id_process, 
                 is_driver_tool, is_driver_process) %>% 
          distinct(mutation_id, cluster_id_tool, 
                   is_driver_process, cluster_id_process)
        
        # Retrieve process subclonal drivers for all clusters
        subclonal_drivers_process = subclonal_table %>% 
          filter(is_driver_process == TRUE) %>% 
          select(mutation_id, cluster_id_tool, cluster_id_process, 
                 is_driver_tool, is_driver_process) %>% 
          distinct(mutation_id, cluster_id_tool, 
                   is_driver_tool, cluster_id_process)
        
        ### Clonal only TP ####
        #### TP ####
        # TP: how many true drivers I have in the TP clusters
        
        clonal_TP_d = clonal_drivers_TP_tool %>%
          summarise(count = sum(is_driver_process, na.rm = TRUE))
        
        clonal_which_TP_d = clonal_drivers_TP_tool %>%
          filter(is_driver_process) %>%
          pull(cluster_id_tool)
        
        #### FP ####
        # FP: how many wrong drivers I have in the TP clusters (i.e. in the interpreted clusters)
        # tool clusters where all rows have is_driver_process == FALSE (i.e., zero TRUE in that cluster).
        clonal_FP_d = clonal_drivers_TP_tool %>%
          summarise(count = sum(!is_driver_process, na.rm = TRUE))
        
        clonal_which_FP_d = clonal_drivers_TP_tool %>%
          filter(!is_driver_process) %>%
          pull(cluster_id_tool)
        
        #### FN ####
        clonal_FN_d = clonal_drivers_TP_process %>%
          summarise(count = sum(!is_driver_tool, na.rm = TRUE))
        
        clonal_which_FN_d = clonal_drivers_TP_process %>%
          filter(!is_driver_tool) %>%
          pull(cluster_id_tool)
        
        ### Sublonal only TP ####
        #### TP ####
        # TP: how many true drivers I have in the TP clusters
        subclonal_TP_d = subclonal_drivers_TP_tool %>%
          summarise(count = sum(is_driver_process, na.rm = TRUE))
        
        subclonal_which_TP_d = subclonal_drivers_TP_tool %>%
          filter(is_driver_process) %>%
          pull(cluster_id_tool)
        
        #### FP ####
        # FP: how many wrong drivers I have in the TP clusters (i.e. in the interpreted clusters)
        # tool clusters where all rows have is_driver_process == FALSE (i.e., zero TRUE in that cluster).
        subclonal_FP_d = subclonal_drivers_TP_tool %>%
          summarise(count = sum(!is_driver_process, na.rm = TRUE))
        
        subclonal_which_FP_d = subclonal_drivers_TP_tool %>%
          filter(!is_driver_process) %>%
          pull(cluster_id_tool)
        
        #### FN ####
        subclonal_FN_d = subclonal_drivers_TP_process %>%
          summarise(count = sum(!is_driver_tool, na.rm = TRUE))
        
        subclonal_which_FN_d = subclonal_drivers_TP_process %>%
          filter(!is_driver_tool) %>%
          pull(cluster_id_tool)
        
        ### Clonal all ####
        #### TP ####
        # TP: how many true drivers I have in the TP clusters
        
        all_clonal_TP_d = clonal_drivers_tool %>%
          summarise(count = sum(is_driver_process, na.rm = TRUE))
        
        all_clonal_which_TP_d = clonal_drivers_tool %>%
          filter(is_driver_process) %>%
          pull(cluster_id_tool)
        
        #### FP ####
        # FP: how many wrong drivers I have in the TP clusters (i.e. in the interpreted clusters)
        # tool clusters where all rows have is_driver_process == FALSE (i.e., zero TRUE in that cluster).
        all_clonal_FP_d = clonal_drivers_tool %>%
          summarise(count = sum(!is_driver_process, na.rm = TRUE))
        
        all_clonal_which_FP_d = clonal_drivers_tool %>%
          filter(!is_driver_process) %>%
          pull(cluster_id_tool)
        
        #### FN ####
        all_clonal_FN_d = clonal_drivers_process %>%
          summarise(count = sum(!is_driver_tool, na.rm = TRUE))
        
        all_clonal_which_FN_d = clonal_drivers_process %>%
          filter(!is_driver_tool) %>%
          pull(cluster_id_tool)
        
        ### Sublonal all ####
        #### TP ####
        # TP: how many true drivers I have in the TP clusters
        
        all_subclonal_TP_d = subclonal_drivers_tool %>%
          summarise(count = sum(is_driver_process, na.rm = TRUE))
        
        all_subclonal_which_TP_d = subclonal_drivers_tool %>%
          filter(is_driver_process) %>%
          pull(cluster_id_tool)
        
        #### FP ####
        # FP: how many wrong drivers I have in the TP clusters (i.e. in the interpreted clusters)
        # tool clusters where all rows have is_driver_process == FALSE (i.e., zero TRUE in that cluster).
        all_subclonal_FP_d = subclonal_drivers_tool %>%
          summarise(count = sum(!is_driver_process, na.rm = TRUE))
        
        all_subclonal_which_FP_d = subclonal_drivers_tool %>%
          filter(!is_driver_process) %>%
          pull(cluster_id_tool)
        
        #### FN ####
        all_subclonal_FN_d = subclonal_drivers_process %>%
          summarise(count = sum(!is_driver_tool, na.rm = TRUE))
        
        all_subclonal_which_FN_d = subclonal_drivers_process %>%
          filter(!is_driver_tool) %>%
          pull(cluster_id_tool)
        
        if(is.na(nmi_interpreted)){
          nmi_interpreted = 1
        }
        
        if(is.na(nmi_blind)){
          nmi_blind = 1
        }
        
        
        df_tmp <- data.frame(
          spn = spn, purity = purity, coverage = coverage,
          vcf_caller = vcf_caller, cna_caller = cna_caller,
          tool = tool,
          sample = sample,
          
          TP_c_blind = TP_c_blind,
          TP_c_blind_list = I(list(unname(which_TP_c_blind))),
          FP_c_blind = FP_c_blind,
          FP_c_blind_list = I(list(unname(which_FP_c_blind))),
          FN_c_blind = FN_c_blind,
          FN_c_blind_list = NA,
          
          TP_c = TP_c,
          TP_c_list = I(list(unname(which_TP_c))),
          FP_c = FP_c,
          FP_c_list = I(list(unname(which_FP_c))),
          FN_c = FN_c,
          FN_c_list = I(list(unname(which_FN_c))),
          
          clonal_TP_d = clonal_TP_d$count,
          clonal_TP_d_list = I(list(unname(clonal_which_TP_d))),
          clonal_FP_d = clonal_FP_d$count,
          clonal_FP_d_list = I(list(unname(clonal_which_FP_d))),
          clonal_FN_d = clonal_FN_d$count,
          clonal_FN_d_list = I(list(unname(clonal_which_FN_d))),
          
          subclonal_TP_d = subclonal_TP_d$count,
          subclonal_TP_d_list = I(list(unname(subclonal_which_TP_d))),
          subclonal_FP_d = subclonal_FP_d$count,
          subclonal_FP_d_list = I(list(unname(subclonal_which_FP_d))),
          subclonal_FN_d = subclonal_FN_d$count,
          subclonal_FN_d_list = I(list(unname(subclonal_which_FN_d))),
          
          all_clonal_TP_d = all_clonal_TP_d$count,
          all_clonal_TP_d_list = I(list(unname(all_clonal_which_TP_d))),
          all_clonal_FP_d = all_clonal_FP_d$count,
          all_clonal_FP_d_list = I(list(unname(all_clonal_which_FP_d))),
          all_clonal_FN_d = all_clonal_FN_d$count,
          all_clonal_FN_d_list = I(list(unname(all_clonal_which_FN_d))),
          
          all_subclonal_TP_d = all_subclonal_TP_d$count,
          all_subclonal_TP_d_list = I(list(unname(all_subclonal_which_TP_d))),
          all_subclonal_FP_d = all_subclonal_FP_d$count,
          all_subclonal_FP_d_list = I(list(unname(all_subclonal_which_FP_d))),
          all_subclonal_FN_d = all_subclonal_FN_d$count,
          all_subclonal_FN_d_list = I(list(unname(all_subclonal_which_FN_d))),
          
          nmi_blind = nmi_blind,
          nmi_interpreted = nmi_interpreted,
          
          n_cluster_blind = n_cluster_blind,
          n_cluster_interpreted = n_cluster_interpreted,
          n_cluster_interpreted_driver = n_cluster_interpreted_driver,
          n_cluster_process = n_cluster_process,
          
          stringsAsFactors = FALSE
        )
      
        df = bind_rows(df_tmp, df)
        }
      
      
    } 
    
    df = df %>%
      mutate(
        precision_c = ifelse(TP_c + FP_c > 0, TP_c / (TP_c + FP_c), NA_real_),
        recall_c    = ifelse(TP_c + FN_c > 0, TP_c / (TP_c + FN_c), NA_real_),
        precision_c_blind = ifelse(TP_c_blind + FP_c_blind > 0, TP_c_blind / (TP_c_blind + FP_c_blind), NA_real_),
        recall_c_blind  = 1)
      
    metrics_drivers = bind_rows(metrics_drivers, df)
  }
}

saveRDS(metrics_drivers, file.path(save_path, "metrics_tables/metrics_drivers_clonal_vs_subclonal.rds"))

metrics_drivers= readRDS(file.path(save_path, "metrics_tables/metrics_drivers_clonal_vs_subclonal.rds"))

