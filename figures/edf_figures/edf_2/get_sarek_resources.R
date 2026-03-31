#rm(list=ls())
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(stringr)
library(lubridate)

get_resources_pipeline <- function(spn,coverage,purity,pipeline){
  comb <- paste0(coverage,"x_",purity,"p")
  if (pipeline=="sarek"){
    execution_dir <- file.path(spn,"sarek",comb,"pipeline_info")
    execution_dir_normal <- file.path(spn,"sarek","normal","pipeline_info")
    execution_traces <- list.files(path = execution_dir,pattern = "execution_trace",recursive = T,full.names = T)
    execution_traces_normal <- list.files(path = execution_dir_normal,pattern = "execution_trace",recursive = T,full.names = T)
    trace_df <- lapply(execution_traces, function(x) {
      df <- read.table(x, header = TRUE, sep = "\t")
      df <- df %>% 
        mutate(purity=purity,
               coverage=coverage)
      
      if (nrow(df) == 0) {
        message("Skipping empty file: ", x)
        return(NULL) 
      }
      return(df)
    })
    trace_df_normal <- lapply(execution_traces_normal, function(x) {
      df <- read.table(x, header = TRUE, sep = "\t")
      df <- df %>% 
        mutate(purity=purity,
               coverage=coverage,
               type="normal") 
      
      if (nrow(df) == 0) {
        message("Skipping empty file: ", x)
        return(NULL) 
      }
      return(df)
    })
    df_final <- do.call("rbind",trace_df)
    df_final_normal <- do.call("rbind",trace_df_normal)
    df_resources <- df_final %>% 
      dplyr::filter(status=="COMPLETED") %>% 
      tidyr::separate(col = "name", into = c("process","tag"),sep = " ") %>% 
      mutate(
        duration_seconds = parse_with_lubridate(duration),
        realtime_seconds = parse_with_lubridate(realtime))
    df_resources_normal <- df_final_normal %>% 
      dplyr::filter(status=="COMPLETED") %>% 
      tidyr::separate(col = "name", into = c("process","tag"),sep = " ") %>% 
      mutate(
        duration_seconds = parse_with_lubridate(duration),
        realtime_seconds = parse_with_lubridate(realtime)
      )
    
    df_resources_sarek <- df_resources %>% 
      filter(grepl("SAREK", process)) %>%
      tidyr::separate(col = "process", into = c("pipeline","subworkflow1","subworkflow2","subworkflow3","module"),sep = ":")
    
    
    df_resources_sequenza <- df_resources %>% 
      filter(grepl("SEQUENZA", process)) %>% 
      tidyr::separate(col = "process", into = c("subworkflow1","module"),sep = ":")
    
    df_resources_battengerg <- df_resources %>% 
      filter(grepl("BATTENBERG", process))
    
    selected_subworkflow3 <- c("GATK4_MARKDUPLICATES","BWAMEM2_MEM",
                               "GATK4_BASERECALIBRATOR","GATK4_APPLYBQSR"
    )
    
    selected_modules <- c("GETPILEUPSUMMARIES_TUMOR","GETPILEUPSUMMARIES_NORMAL",
                          "GATHERPILEUPSUMMARIES_TUMOR",
                          "FREEBAYES",
                          "CNVKIT_BATCH","CNVKIT_CALL","CNVKIT_GENEMETRICS","CNVKIT_EXPORT",
                          "ASCAT",
                          "STRELKA_SOMATIC","STRELKA_SINGLE","MUTECT2_PAIRED",
                          "CALCULATECONTAMINATION","LEARNREADORIENTATIONMODEL",
                          "FILTERMUTECTCALLS"
    )
    df_resources_sarek <- df_resources_sarek %>%
      mutate(
        peak_rss_gb = case_when(
          str_detect(peak_rss, "GB") ~ as.numeric(str_remove(peak_rss, " GB")),
          str_detect(peak_rss, "MB") ~ as.numeric(str_remove(peak_rss, " MB")) / 1024
        )
      ) %>% 
      mutate(
        peak_vmem_gb = case_when(
          str_detect(peak_vmem, "GB") ~ as.numeric(str_remove(peak_vmem, " GB")),
          str_detect(peak_vmem, "MB") ~ as.numeric(str_remove(peak_vmem, " MB")) / 1024
        )
      ) %>% 
      mutate(X.cpu = as.numeric(gsub("%", "", X.cpu)))
    
    df_resources_sequenza <- df_resources_sequenza %>%
      mutate(
        peak_rss_gb = case_when(
          str_detect(peak_rss, "GB") ~ as.numeric(str_remove(peak_rss, " GB")),
          str_detect(peak_rss, "MB") ~ as.numeric(str_remove(peak_rss, " MB")) / 1024
        )
      ) %>% 
      mutate(
        peak_vmem_gb = case_when(
          str_detect(peak_vmem, "GB") ~ as.numeric(str_remove(peak_vmem, " GB")),
          str_detect(peak_vmem, "MB") ~ as.numeric(str_remove(peak_vmem, " MB")) / 1024
        )
      ) %>% 
      mutate(X.cpu = as.numeric(gsub("%", "", X.cpu)))
    
    df_resources_battengerg <- df_resources_battengerg %>%
      mutate(
        peak_rss_gb = case_when(
          str_detect(peak_rss, "GB") ~ as.numeric(str_remove(peak_rss, " GB")),
          str_detect(peak_rss, "MB") ~ as.numeric(str_remove(peak_rss, " MB")) / 1024
        )
      ) %>% 
      mutate(
        peak_vmem_gb = case_when(
          str_detect(peak_vmem, "GB") ~ as.numeric(str_remove(peak_vmem, " GB")),
          str_detect(peak_vmem, "MB") ~ as.numeric(str_remove(peak_vmem, " MB")) / 1024
        )
      ) %>% 
      mutate(X.cpu = as.numeric(gsub("%", "", X.cpu)))
    
    df_sarek_preprocess <- df_resources_sarek %>%
      filter(subworkflow3%in%selected_subworkflow3) %>% 
      group_by(subworkflow3) %>% 
      summarise(
        mean_cpus = mean(X.cpu),
        sd_cpus = sd(X.cpu),
        mean_peak_rss = mean(peak_rss_gb),
        sd_peak_rss = sd(peak_rss_gb),
        mean_peak_vmem = mean(peak_vmem_gb),
        sd_peak_vmem = sd(peak_vmem_gb)
      ) %>% 
      mutate(process=subworkflow3,
             coverage = coverage,
             purity = purity,
             step = "preprocess"
      ) %>% 
      select(process,mean_cpus,sd_cpus,mean_peak_rss,sd_peak_rss,mean_peak_vmem,sd_peak_vmem,coverage,purity,step)
    
    df_sarek_variant_calling <- df_resources_sarek %>%
      filter(module%in%selected_modules) %>% 
      group_by(module) %>% 
      summarise(
        mean_cpus = mean(X.cpu),
        sd_cpus = sd(X.cpu),
        mean_peak_rss = mean(peak_rss_gb),
        sd_peak_rss = sd(peak_rss_gb),
        mean_peak_vmem = mean(peak_vmem_gb),
        sd_peak_vmem = sd(peak_vmem_gb)
      ) %>% 
      mutate(process=module,
             coverage = coverage,
             purity = purity,
             step = "variant_calling"
      )  %>% 
      select(process,mean_cpus,sd_cpus,mean_peak_rss,sd_peak_rss,mean_peak_vmem,sd_peak_vmem,coverage,purity,step)
    
    df_sarek_sequenza <- df_resources_sequenza %>%
      group_by(module) %>% 
      summarise(
        mean_cpus = mean(X.cpu),
        sd_cpus = sd(X.cpu),
        mean_peak_rss = mean(peak_rss_gb),
        sd_peak_rss = sd(peak_rss_gb),
        mean_peak_vmem = mean(peak_vmem_gb),
        sd_peak_vmem = sd(peak_vmem_gb)
      )  %>% 
      mutate(process=module,
             coverage = coverage,
             purity = purity,
             step = "variant_calling"
      ) %>% 
      select(process,mean_cpus,sd_cpus,mean_peak_rss,sd_peak_rss,mean_peak_vmem,sd_peak_vmem,coverage,purity,step)
    df_sarek_battenberg <- df_resources_battengerg %>%
      group_by(process) %>% 
      summarise(
        mean_cpus = mean(X.cpu),
        sd_cpus = sd(X.cpu),
        mean_peak_rss = mean(peak_rss_gb),
        sd_peak_rss = sd(peak_rss_gb),
        mean_peak_vmem = mean(peak_vmem_gb),
        sd_peak_vmem = sd(peak_vmem_gb)
      ) %>% 
      mutate(process=process,
             coverage = coverage,
             purity = purity,
             step = "variant_calling"
      ) %>% 
      select(process,mean_cpus,sd_cpus,mean_peak_rss,sd_peak_rss,mean_peak_vmem,sd_peak_vmem,coverage,purity,step)
    
    
    
    df_final <- rbind(df_sarek_preprocess,df_sarek_variant_calling,df_sarek_sequenza,df_sarek_battenberg)
    df_final <- df_final %>% 
      dplyr::mutate(SPN=spn)
  } else if (pipeline=="tumourevo"){
    vcf_caller = "mutect2"
    cna_caller = "ascat"
    comb <- paste0(coverage,"x_",purity,"p_",vcf_caller,"_",cna_caller)
    execution_dir <- file.path(spn,"tumourevo",comb,"pipeline_info")
    execution_traces <- list.files(path = execution_dir,pattern = "execution_trace",recursive = T,full.names = T)
    
    trace_df <- lapply(execution_traces, function(x) {
      df <- read.table(x, header = TRUE, sep = "\t")
      df <- df %>% 
        mutate(purity=purity,
               coverage=coverage)
      
      if (nrow(df) == 0) {
        message("Skipping empty file: ", x)
        return(NULL) 
      }
      return(df)
    })
    df_final <- do.call("rbind",trace_df)
    
    df_resources <- df_final %>% 
      dplyr::filter(status=="COMPLETED") %>% 
      tidyr::separate(col = "name", into = c("process","tag"),sep = " ") %>% 
      mutate(
        duration_seconds = parse_with_lubridate(duration),
        realtime_seconds = parse_with_lubridate(realtime)
      )
    df_resources_tumourevo <- df_resources %>% 
      tidyr::separate(col = "process", into = c("pipeline","subworkflow1","subworkflow2","subworkflow3","module"),sep = ":")
    
    selected_subworkflow3 <- c("ENSEMBLVEP_VEP","ANNOTATE_DRIVER",
                               "TINC","CNAQC","JOIN_CNAQC","VIBER","CTREE_VIBER" ,
                               "MOBSTERh","CTREE_MOBSTERh",
                               "PYCLONEVI","CTREE_PYCLONEVI",
                               "SIGPROFILER","SPARSE_SIGNATURES"
    )
    df_tumourevo <- df_resources_tumourevo %>%
      filter(subworkflow3%in%selected_subworkflow3) %>% 
      group_by(subworkflow2) %>% 
      summarise(
        mean_duration = mean(as.numeric(duration_seconds, units = "secs")),
        st_duration = sd(as.numeric(duration_seconds, units = "secs"))
      ) %>% 
      mutate(process=subworkflow2,
             coverage = coverage,
             purity = purity,
             step = "tumourevo"
      ) %>% 
      select(process,mean_duration,st_duration,coverage,purity,step)
    df_final <- df_tumourevo
    df_final <- df_final %>% 
      dplyr::mutate(SPN=spn)
  } 
  return(df_final)
}