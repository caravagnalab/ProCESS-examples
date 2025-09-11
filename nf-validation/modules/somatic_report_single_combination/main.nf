process GENERATE_SOMATIC_REPORT_SINGLE_COMBINATION {
    tag "${meta.spn}_${meta.coverage}_${meta.purity}_${caller}"

    input:
    tuple val(meta), val(rds_process_list), val(process_sample), val(process), val(rds_caller_list), val(caller_sample), val(caller)

    output:
    tuple val(meta), path("**/*.rds"),	emit:metrics_somatic_rds
    tuple val(meta), path("**/*.png"),	emit:metrics_somatic_report
    
    publishDir "${params.outdir}/${meta.spn}/somatic/${meta.coverage}x_${meta.purity}p/", mode: 'copy'

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
    
    process_samples = substr("$process_sample", 2, nchar("$process_sample")-1) 
    process_samples = strsplit(process_samples, ", ")[[1]]

    caller_samples = substr("$caller_sample", 2, nchar("$caller_sample")-1) 
    caller_samples = strsplit(caller_samples, ", ")[[1]]

    caller = "${caller}"

    combination = paste0(coverage, "x_", purity, "p")

    if (gender=="XX"){
      chromosomes = c(paste0('chr',1:22), 'chrX')
    } else {
      chromosomes = c(paste0('chr',1:22), 'chrX', 'chrY')
    }
    
    # INPUT PARAMATERS ####
    min_ccf = .02
    min_vaf_caller = 0.0
    mut_types = c("INDEL", "SNV")
    
    process_rds_per_sample = strsplit("${rds_process_list}", '], ') %>% unlist()
    process_rds_sample = lapply(1:length(process_rds_per_sample), FUN = function(i){
      rds = process_rds_per_sample[[i]] %>% 
            stringr::str_replace_all('\\\\[', '') %>% 
            stringr::str_replace_all('\\\\]', '') %>% 
            strsplit(', ') %>% 
            unlist()
    })
    names(process_rds_sample) = process_samples

    caller_rds_per_sample = strsplit("${rds_caller_list}", '], ') %>% unlist()
    caller_rds_sample = lapply(1:length(caller_rds_per_sample), FUN = function(i){
      rds = caller_rds_per_sample[[i]] %>% 
            stringr::str_replace_all('\\\\[', '') %>% 
            stringr::str_replace_all('\\\\]', '') %>% 
            strsplit(', ') %>% 
            unlist()
    })
    names(caller_rds_sample) = caller_samples

    for (mut_type in mut_types) {
      print(mut_type)
      
      gt_res = lapply(process_samples, function(sample_id) {
        gt_res = lapply(process_rds_sample[[sample_id]], function(rds){
          readRDS(rds)[[mut_type]]
        }) %>% do.call("bind_rows", .) %>% 
          dplyr::mutate(sample = sample_id) %>% 
          dplyr::mutate(mutationID = paste0(mutationID,":",sample_id))
        gt_res
      }) %>% do.call("bind_rows", .)

      caller_res = lapply(caller_samples, function(sample_id) {
        caller_res = lapply(caller_rds_sample[[sample_id]], function(rds) {
          readRDS(rds)[[mut_type]] %>% 
            dplyr::filter(FILTER == "PASS")
        }) %>% do.call("bind_rows", .) %>% 
          dplyr::filter(!is.na(VAF)) %>% 
          dplyr::mutate(sample = sample_id) %>% 
          dplyr::mutate(mutationID = paste0(mutationID,":",sample_id))
        caller_res
      }) %>% do.call("bind_rows", .)

      sample_info = list(caller_name=caller, mut_type=mut_type, spn=spn_id, purity=purity, coverage=coverage)
      report = get_report(seq_res_long = gt_res, 
                          caller_res = caller_res, 
                          sample_info = sample_info, 
                          min_ccf = min_ccf, min_vaf_caller = min_vaf_caller)
        
      results_folder_path = file.path(caller, "report", mut_type)
      dir.create(results_folder_path, recursive = T)
      metrics_path = file.path(results_folder_path, "metrics.rds")
      saveRDS(list(report_metrics=report\$report_metrics, vaf_comparison=report\$vaf_comparison), metrics_path)
      
      filename = paste(spn_id, combination, caller, mut_type, sep = '_')
      file_path = file.path(results_folder_path, filename)
      ggsave(paste0(file_path, '.png'), plot = report\$report_plot, width = 18, height = 18, units = "in", dpi = 400)
      message("        Report done!")

    }


    """
}