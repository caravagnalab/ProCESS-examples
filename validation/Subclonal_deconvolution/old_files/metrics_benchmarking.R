library(tidyverse)
library(randnet)
library(scales)

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

source(file.path(save_path, "utils_plots_final.R"))

coverage_list = c(50,100,150)
purity_list = c(0.3,0.6,0.9)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN06', 'SPN07')
tools = c('viber', 'pyclonevi', 'mobster')
tools = c('viber', 'pyclonevi')

combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller=cna_caller_list,
                    spn=spn_list)

spn = 'SPN07'
coverage=50
purity=0.6
vcf_caller = "mutect2"
cna_caller = "ascat"
simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)

metrics_table = data.frame()

for(tool in tools){
  print(tool)
  for(i in 1:nrow(combs)){
    coverage = combs[i, "coverage"]
    purity = combs[i, "purity"]
    vcf_caller = combs[i, "vcf_caller"]
    cna_caller = combs[i, "cna_caller"]
    spn = combs[i, "spn"]
    
    simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
    
    print(paste0(spn, "_", coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller))
    
    table = readRDS(file.path(main_path, "subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds")))
    
    if(nrow(table%>%
            pivot_wider(
              names_from = "sample_id",
              values_from = "vaf_process"
            )) < 5000){
      next
    }
    if(tool == 'mobster'){
      
      table = table %>% 
        group_by(cluster_id_process, sample_id) %>%
        mutate(ccf_process=mean(ccf_process, na.rm=TRUE)) %>%
        ungroup()
      
      table = table %>% group_by(cluster_id_process, sample_id) %>% 
        mutate(is_clonal_process=replace(FALSE, ccf_process > 0.95, TRUE)) %>% ungroup() %>% 
        mutate(cluster_id_process_full = cluster_id_process) %>% 
        mutate(cluster_id_process=replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal'))
      
      sample_names = table$sample_id %>% unique()
      
      # Create the df to fill sample per sample
      df = data.frame(matrix(nrow = 0, ncol = length(columns)))
      colnames(df) = columns
      
      for(sample_name in sample_names){
        
        table_sample = table %>% filter(sample_id == sample_name)
       
        # Kolmogorov-Smirnov
        ccf_tool = table_sample %>% 
          select(sample_id, cluster_id_tool, ccf_tool) %>% 
          unique() %>% arrange(sample_id) %>% pull(ccf_tool)
        
        ccf_process = table_sample %>% 
          select(sample_id, cluster_id_process, ccf_process) %>% 
          unique() %>% arrange(sample_id) %>% pull(ccf_process)
        
        ccf_ks = ks.test(ccf_tool, ccf_process)$statistic
        
        
        # now I can compute the NMI and the metrics and put them in the dataframe
        # Number of clusters tool and process
        n_raw_clones_tool = length(unique(table_sample$cluster_id_tool))
        n_true_driver_clones_process = length(unique(table_sample$cluster_id_process))
        
        # NMI and ARI
        nmi_raw = randnet::NMI(as.factor(table$cluster_id_tool),
                               as.factor(table$cluster_id_process))
        
        ari_raw = aricode::ARI(as.factor(table$cluster_id_tool),
                               as.factor(table$cluster_id_process))
        
        # Absolute error
        # Relative error
        
        df_samples = data.frame(spn = spn, purity = purity, 
                        coverage = coverage, vcf_caller = vcf_caller, 
                        cna_caller=cna_caller,
                        nmi_raw=nmi_raw,
                        nmi_interpreted=NA,
                        nmi_interpreted_driver=NA,
                        ari_interpreted=NA,
                        ari_raw=ari_raw,
                        n_raw_tool=n_raw_clones_tool,
                        n_true_driver_process = n_true_driver_clones_process,
                        n_interpreted_tool = NA,
                        n_interpreted_driver_tool = NA,
                        tool=tool,
                        n_interpreted_tool_no_tail=NA,
                        n_true_driver_process_no_tail=NA,
                        nmi_interpreted_no_tail=NA,
                        ari_interpreted_no_tail=NA,
                        Kolmogorov_distance = ccf_ks,
                        wasserstein_raw=NA,
                        wasserstein_interpreted=NA)
        
        df = rbind(df, df_samples)
      }
    }
    else{
      table = table %>% 
        group_by(cluster_id_process, sample_id) %>%
        mutate(ccf_process=mean(ccf_process, na.rm=TRUE)) %>%
        ungroup() %>% 
        add_count(cluster_id_process, name="n_mutations_process") %>%
        add_count(cluster_id_tool, name="n_mutations_tool") %>%
        add_count(cluster_id_tool_interpreted, name="n_mutations_tool_interpreted")
      
      table = table %>% group_by(cluster_id_process) %>% 
        mutate(is_clonal_process=replace(FALSE, all(ccf_process > 0.95), TRUE)) %>% ungroup() %>% 
        mutate(cluster_id_process_full = cluster_id_process) %>% 
        mutate(cluster_id_process=replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal'))
      
      # --------- Salva tabella con nuovi cluster --------- #
      # saveRDS(table, file.path(save_path, "tables_interpreted_new_clusters", paste0(tool, "_", spn, "_", simulation_id, ".rds")))
      # ------------ #
      
      # Kolmogorov - Smirnov
      ccf_tool = table %>% 
        select(sample_id, cluster_id_tool, ccf_tool) %>% 
        unique() %>% arrange(sample_id) %>% pull(ccf_tool)
      
      
      ccf_tool_interpreted = table %>% 
        filter(cluster_id_tool_interpreted != 'Subclonal') %>% 
        select(sample_id, cluster_id_tool_interpreted, ccf_tool) %>% 
        unique() %>% arrange(sample_id) %>% pull(ccf_tool)
      
      ccf_process = table %>% 
        select(sample_id, cluster_id_process, ccf_process) %>% 
        unique() %>% arrange(sample_id) %>% pull(ccf_process)
      
      ccf_process_interpreted = table %>% 
        filter(cluster_id_process != 'Subclonal') %>% 
        select(sample_id, cluster_id_process, ccf_process) %>% 
        unique() %>% arrange(sample_id) %>% pull(ccf_process)
      
      ccf_ks = ks.test(ccf_tool, ccf_process)$statistic
      # ccf_ks_int = ks.test(ccf_tool_interpreted, ccf_process_interpreted)$statistic
      
      was_raw = transport::wasserstein1d(ccf_tool, ccf_process)
      
      if(length(ccf_tool_interpreted)>0){
      was_interpreted = transport::wasserstein1d(ccf_tool_interpreted, ccf_process_interpreted)
      }else{
        was_interpreted=NA
      }
      # ccf_ks = table %>%
      #   group_by(patient_id, coverage, purity, tool) %>%
      #   dplyr::select(patient_id, sample_id, coverage, purity, tool, contains("ccf")) %>%
      #   unique() %>%
      #   summarise(ccf_ks=ks.test(ccf_tool, ccf_process)$statistic) %>%
      #   ungroup() %>% pull(ccf_ks)
      
      
      # NMI and ARI
      
      n_raw_clones_tool = length(unique(table$cluster_id_tool))
      n_interpreted_driver_clones_tool = length(unique(table$cluster_id_tool_interpreted))
      n_interpreted_driver_clones_process = length(unique(table$cluster_id_tool_interpreted_driver))
      n_true_driver_clones_process = length(unique(table$cluster_id_process))
        
      nmi_interpreted = randnet::NMI(as.factor(table$cluster_id_tool_interpreted),
                                     as.factor(table$cluster_id_process))
      
      ari_interpreted = aricode::ARI(as.factor(table$cluster_id_tool_interpreted),
                                     as.factor(table$cluster_id_process))
      
      nmi_raw = randnet::NMI(as.factor(table$cluster_id_tool),
                                     as.factor(table$cluster_id_process))
      
      ari_raw = aricode::ARI(as.factor(table$cluster_id_tool),
                                     as.factor(table$cluster_id_process))
      
      nmi_interpreted_driver = randnet::NMI(as.factor(table$cluster_id_tool_interpreted_driver),
                                     as.factor(table$cluster_id_process))
      
      # Now we need to only consider clusters != 'Subclonal' in tool
      table_no_tail = table %>% filter(cluster_id_tool_interpreted != 'Subclonal')
      
      if((table_no_tail %>% nrow()) >0){
        nmi_interpreted_no_tail = randnet::NMI(as.factor(table_no_tail$cluster_id_tool),
                               as.factor(table_no_tail$cluster_id_process))
        
        ari_interpreted_no_tail = aricode::ARI(as.factor(table_no_tail$cluster_id_tool_interpreted),
                                       as.factor(table_no_tail$cluster_id_process))
        
        n_interpreted_clones_tool_no_tail = length(unique(table_no_tail$cluster_id_tool_interpreted))
        n_true_driver_clones_process_no_tail = length(unique(table_no_tail$cluster_id_process))
      }
      else{
        nmi_interpreted_no_tail=0
        ari_interpreted_no_tail=0
        n_interpreted_clones_tool_no_tail=0
        n_true_driver_clones_process_no_tail=0
      }
      df = data.frame(spn = spn, purity = purity, 
                      coverage = coverage, vcf_caller = vcf_caller, 
                      cna_caller=cna_caller,
                      n_raw_tool=n_raw_clones_tool,
                      nmi_raw=nmi_raw,
                      nmi_interpreted=nmi_interpreted,
                      nmi_interpreted_driver=nmi_interpreted_driver,
                      ari_interpreted=ari_interpreted,
                      ari_raw=ari_raw,
                      nmi_interpreted_no_tail=nmi_interpreted_no_tail,
                      ari_interpreted_no_tail=ari_interpreted_no_tail,
                      n_true_driver_process = n_true_driver_clones_process,
                      n_interpreted_tool = n_interpreted_driver_clones_tool,
                      n_interpreted_driver_tool=n_interpreted_driver_clones_process,
                      n_interpreted_tool_no_tail=n_interpreted_clones_tool_no_tail,
                      n_true_driver_process_no_tail=n_true_driver_clones_process_no_tail,
                      tool=tool,
                      Kolmogorov_distance = ccf_ks,
                      wasserstein_raw=was_raw,
                      wasserstein_interpreted=was_interpreted)
    }
    metrics_table = bind_rows(metrics_table, df)
    # metrics_table = rbind(metrics_table, df)
  }
}

saveRDS(metrics_table, file.path(save_path, "metrics_tables/new_table_clusters_metrics_v2.rds"))