process CNAQC_VALIDATION_COMBINATION {
    tag "${meta.spn}_${meta.coverage}_${meta.purity}-cnaqc_validation"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/process_validation:v2' :
        'docker.io/lvaleriani/process_validation:v2' }"

    input:
    tuple val(meta), val(path), val(vcf_caller), val(cna_caller)

    output:
    tuple val(meta), path("cnaqc_validate.rds"),             emit:comb_cnaqc_validate
    tuple val(meta), path("cnaqc_stats.rds"),                emit:comb_cnaqc_stats
    tuple val(meta), path("report_cnaqc_validation.png"),    emit:comb_cnaqc_validate_report

    
    publishDir "${params.outdir}/${meta.spn}/cnaqc/${meta.coverage}x_${meta.purity}/${vcf_caller}_${cna_caller}", mode: 'copy'

    script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(patchwork)
    library(ProCESS)
    library(CNAqc)
    library(dplyr)
    



    source("${projectDir}/bin/getters/process_getters.R")
    source("${projectDir}/bin/getters/tumourevo_getters.R")
    source("${projectDir}/bin/qc/utils.R")




    data_dir = "${path}"
    spn = "${meta.spn}"
    coverage = as.numeric("${meta.coverage}")
    purity = as.numeric("${meta.purity}")

    vcf_caller <- "${vcf_caller}"
    cna_caller <- "${cna_caller}"
    callers <- paste0(vcf_caller,"_",cna_caller)
    
    samples = get_sample_names(spn = spn)
    
    cnaqc_obj <- lapply(samples, FUN = function(sample){
      file <- get_tumourevo_qc(spn = spn, coverage = coverage, purity = purity, tool = 'cnaqc', 
                               vcf_caller = vcf_caller, cna_caller = cna_caller, sample = sample)
      if (length(file) > 0){
        if (file.exists(file\$qc_rds)){
          readRDS(file\$qc_rds)
        }
      }
    })
    names(cnaqc_obj) <- samples
    if (!is.null(cnaqc_obj[[1]])){
        
        # get right QC segments
        true_cna <- lapply(samples, FUN = function(sample){ 
          readRDS(get_process_cna(spn, sample = sample)) %>% 
            mutate(sample = sample) %>% 
            filter(ratio > 0.1) %>% 
            mutate(chr = paste0('chr', chr)) %>% 
            dplyr::rename(true_major = major,
                          true_minor = minor)
    }) %>% bind_rows()   
    
    cna_cnaqc <- lapply(samples, FUN = function(sample){ 
        p <- cnaqc_obj[[sample]][['purity']]
        cnaqc_obj[[sample]][['cna']] %>% mutate(sample = sample, 
                                                cn_purity = p)
    }) %>% bind_rows()

    thr = 1e6
    compare <- true_cna %>%
        left_join(cna_cnaqc, 
                  by = join_by(chr,sample), 
                  relationship = "many-to-many") %>% 
        filter(from >= begin-thr & to <= end+thr | begin-thr >= from & end <= to+thr) %>% 
        mutate(true_purity = purity, delta_purity = abs(true_purity-cn_purity)) 
  
    validate <- compare %>% 
        mutate(CNAqc = 'NA') %>%
        mutate(ratio = ifelse(ratio > 0.9, 1, ratio)) %>% 
        mutate(CNAqc = ifelse(paste(Major, minor,sep = ':') == paste(true_major,true_minor,sep = ':')  & delta_purity < 0.1 & QC_PASS == TRUE & ratio > 0.9, 'CNAqc OK - Caller OK', CNAqc)) %>%
        mutate(CNAqc = ifelse(paste(Major, minor,sep = ':') == paste(true_major,true_minor,sep = ':')  & delta_purity < 0.1 & QC_PASS == FALSE & ratio > 0.9, 'CNAqc FAIL - Caller OK', CNAqc)) %>%
        mutate(CNAqc = ifelse(ratio < 0.9, 'Subclonal', CNAqc)) %>%
        
        mutate(CNAqc = ifelse(paste(Major, minor,sep = ':') != paste(true_major,true_minor,sep = ':')  & delta_purity < 0.1 & QC_PASS == FALSE & ratio > 0.9, 'CNAqc OK - Caller FAIL', CNAqc)) %>%
        mutate(CNAqc = ifelse(paste(Major, minor,sep = ':') == paste(true_major,true_minor,sep = ':')  & delta_purity > 0.1 & QC_PASS == FALSE & ratio > 0.9, 'CNAqc OK - Caller FAIL', CNAqc)) %>%
        mutate(CNAqc = ifelse(paste(Major, minor,sep = ':') != paste(true_major,true_minor,sep = ':')  & delta_purity > 0.1 & QC_PASS == FALSE & ratio > 0.9, 'CNAqc OK - Caller FAIL', CNAqc)) %>%
        
        mutate(CNAqc = ifelse(paste(Major, minor,sep = ':') != paste(true_major,true_minor,sep = ':')  & delta_purity < 0.1 & QC_PASS == TRUE & ratio > 0.9, 'CNAqc FAIL - Caller FAIL', CNAqc)) %>%
        mutate(CNAqc = ifelse(paste(Major, minor,sep = ':') != paste(true_major,true_minor,sep = ':')  & delta_purity > 0.1 & QC_PASS == TRUE & ratio > 0.9, 'CNAqc FAIL - Caller FAIL', CNAqc)) %>%
        mutate(CNAqc = ifelse(paste(Major, minor,sep = ':') == paste(true_major,true_minor,sep = ':')  & delta_purity > 0.1 & QC_PASS == TRUE & ratio > 0.9, 'CNAqc FAIL - Caller FAIL', CNAqc)) %>%
        select(-CCF, -length, -segment_id, -n, -begin, -end) %>%
        distinct() %>% 
        filter(chr != 'chrX') %>% 
        mutate(true_karyo = paste(true_major, true_minor, sep = ':'))
     
    stats_table <- get_statistics_qc(cnaqc_obj, purity = purity)
    table <- stats_table\$table %>% tibble() %>% select(-gg) %>% filter(info != 'sample') %>% 
        tidyr::pivot_wider(values_from = values, names_from = info, names_repair = 'minimal') %>%
        mutate(spn = spn, coverage = coverage, vcf_caller = vcf_caller, cna_caller = cna_caller)
    validate <- validate %>% mutate(spn = spn, coverage = coverage, 
                                    purity = purity, vcf_caller = vcf_caller, cna_caller = cna_caller)

 
    }
    validate_table_filt <- validate %>%
      mutate(comb_te=paste0(vcf_caller,"_",cna_caller)) %>% 
      mutate(CNAqc = ifelse(is.na(CNAqc), 'NA', CNAqc)) %>% 
      group_by(sample, coverage, purity, true_karyo,comb_te) %>%
      filter(!(all(CNAqc == "Subclonal", na.rm = T))) %>%
      ungroup()
    
    plt <- validate_table_filt %>% 
      group_by(true_karyo, CNAqc, sample) %>% 
      summarize(n = n())  %>% 
      ggplot() + 
      geom_col(aes(x = true_karyo, y=n, fill = CNAqc), position = "fill") + # position = "fill"
      ggplot2::scale_fill_manual('', values = c(
        `CNAqc OK - Caller OK` = 'seagreen',
        `CNAqc OK - Caller FAIL` = 'darkseagreen',
        `CNAqc FAIL - Caller FAIL` = 'indianred3', 
        `CNAqc FAIL - Caller OK` = 'tomato4', 
        `Subclonal` = 'gainsboro',
        `NA` = 'gray60'
      )) +     
      ylab('% of segments') +
      facet_wrap(~sample) +
      theme_light() + theme(legend.position = 'bottom', axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))
      
    saveRDS(object=validate,file="cnaqc_validate.rds")
    saveRDS(object=table,file="cnaqc_stats.rds")
    ggsave(filename="report_cnaqc_validation.png",plot=plt)
    """
}