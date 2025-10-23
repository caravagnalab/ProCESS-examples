process GENERATE_SOMATIC_REPORT_ALL_CALLER {
    tag "${meta.spn}_${meta.coverage}_${meta.purity}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/process_validation:v2' :
        'docker.io/lvaleriani/process_validation:v2' }"

    input:
    tuple val(meta), val(rds_process_list), val(process_sample), val(process), val(rds_caller_list), val(caller_sample), val(caller)

    output:
    tuple val(meta), path("**/*.rds"),	emit:metrics_somatic_all_caller_rds
    tuple val(meta), path("**/*.png"),	emit:metrics_somatic_all_caller_report
    
    publishDir "${params.outdir}/${meta.spn}/somatic/${meta.coverage}x_${meta.purity}p/report/", mode: 'copy'

    script:
    """
  
    #!/usr/bin/env Rscript
    
    rm(list = ls())
    options(bitmapType='cairo')
    require(tidyverse)
    library(optparse)
    
    source("${projectDir}/bin/getters/sarek_getters.R")
    source("${projectDir}/bin/getters/process_getters.R")
    source("${projectDir}/bin/somatic/utils/utils.R")
    source("${projectDir}/bin/somatic/utils/plot_utils.R")

    spn_id = "${meta.spn}"
    coverage = "${meta.coverage}"
    purity = "${meta.purity}"
    gender = "${meta.sex}"

    mut_types = c("INDEL", "SNV")
    # INPUT PARAMATERS ####
    # min_vaf = .02
    min_ccf = 0.0
    min_vaf_caller = 0.0
    mut_types = c("INDEL", "SNV")
    
    process_samples = substr("$process_sample", 2, nchar("$process_sample")-1) 
    process_samples = strsplit(process_samples, ", ")[[1]]

    caller = substr("$caller", 2, nchar("$caller")-1) 
    caller = strsplit(caller, ", ")[[1]]

    caller_samples = strsplit("$caller_sample", '], ') %>% 
                      unlist() %>% 
                      stringr::str_replace_all('\\\\[', '') %>% 
                      stringr::str_replace_all('\\\\]', '') %>% 
                      strsplit(', ')
    names(caller_samples) = caller

    combination = paste0(coverage, "x_", purity, "p")
    
    process_rds_per_sample = strsplit("${rds_process_list}", '], ') %>% unlist()
    process_rds_sample = lapply(1:length(process_rds_per_sample), FUN = function(i){
      rds = process_rds_per_sample[[i]] %>% 
            stringr::str_replace_all('\\\\[', '') %>% 
            stringr::str_replace_all('\\\\]', '') %>% 
            strsplit(', ') %>% 
            unlist()
    })
    names(process_rds_sample) = process_samples

    caller_rds_per_caller = strsplit("${rds_caller_list}", ']], \\\\[\\\\[') %>% unlist()
    names(caller_rds_per_caller) = caller
    caller_rds_per_caller_sample = lapply(names(caller_rds_per_caller), FUN = function(call){
        df =  caller_rds_per_caller[[call]] %>% 
                strsplit('], ') %>% 
                unlist() %>%  
                stringr::str_replace_all('\\\\[', '') %>% 
                stringr::str_replace_all('\\\\]', '') %>% 
                strsplit(', ') 
        names(df) = caller_samples[[call]]
        return(df)
    })
    names(caller_rds_per_caller_sample) = caller

    for (mut_type in mut_types) {
      print(mut_type)

      results = lapply(process_samples, function(sample_id) { 
        gt_res = lapply(process_rds_sample[[sample_id]], function(rds){
          readRDS(rds)[[mut_type]] 
        }) %>% do.call("bind_rows", .) %>% 
          dplyr::mutate(sample = sample_id) %>% 
          dplyr::mutate(mutationID = paste0(mutationID,":",sample_id))
        
        caller_res_list = lapply(caller, function(call) {
          data = caller_rds_per_caller_sample[[call]][[sample_id]]
          caller_res = lapply(data, function(rds) {
            readRDS(rds)[[mut_type]] %>% 
              dplyr::filter(FILTER == "PASS")
            }) %>% 
            do.call("bind_rows", .) %>% 
            dplyr::filter(!is.na(VAF)) %>% 
            dplyr::mutate(sample = sample_id) %>% 
            dplyr::mutate(mutationID = paste0(mutationID,":",sample_id))
          caller_res
        })
        names(caller_res_list) = caller
        
        list(gt = gt_res, callers = caller_res_list)
      })
      names(results) = process_samples

      gt_res = lapply(results, function(x) x\$gt) %>% do.call("bind_rows", .)
      caller_res_list = lapply(caller, function(call) {
        lapply(results, function(x) x\$callers[[call]]) %>% do.call("bind_rows", .)
      })
      names(caller_res_list) = caller

      sample_info = list(mut_type=mut_type, spn=spn_id, purity=purity, coverage=coverage)
      report_plot = get_multi_caller_report(seq_res_long = gt_res, 
                                       caller_res_list = caller_res_list, 
                                       # pi = as.numeric(purity),
                                       sample_info = sample_info, 
                                       min_ccf = min_ccf,
                                       min_vaf_caller = min_vaf_caller,
                                       # min_vaf = min_vaf, 
                                       only_pass = TRUE)

      results_folder_path = file.path(mut_type)
      dir.create(results_folder_path, recursive = T)
      filename = paste(combination, "allCaller", mut_type, sep = '_')
      file_path = file.path(results_folder_path, filename)
      ggsave(paste0(file_path, '.png'), plot = report_plot, width = 15, height = 20, units = "in", dpi = 400)

      metrics = lapply(names(caller_res_list), function(nc) {
          caller_res = caller_res_list[[nc]]
          samples_results = lapply(process_samples, function(sample_id) {
            sample_ground_truth_res = gt_res %>% dplyr::filter(sample == sample_id)
            sample_caller_res = caller_res %>% dplyr::filter(sample == sample_id)
            
            # metrics_results = analyze_vaf_performance(sample_ground_truth_res, 
            #                                           sample_caller_res, 
            #                                           only_pass = TRUE, 
            #                                           min_vaf_threshold = min_vaf, 
            #                                           vaf_tolerance_pct = 5)  # 5% VAF tolerance
            
            metrics_results = analyze_ccf_performance(sample_ground_truth_res, 
                                                      sample_caller_res, 
                                                      tolerance_pct = 5, 
                                                      min_vaf_caller = min_vaf_caller,
                                                      min_ccf_threshold = min_ccf,
                                                      only_pass = TRUE)
            
            metrics_results\$raw_data = NULL
            metrics_results
          })
          names(samples_results) = process_samples
          samples_results
        })
        names(metrics) = names(caller_res_list)
        metrics_path = file.path(results_folder_path, "metrics.rds")
        saveRDS(metrics, metrics_path)
    }
    """
}




