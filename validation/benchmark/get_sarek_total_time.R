rm(list=ls())
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(stringr)
library(lubridate)
############# AUXILIARY FUNCTIONS ##############
parse_with_lubridate <- function(x) {
  # Extract numbers or zero if missing
  h  <- as.numeric(str_extract(x, "\\d+(?=h)"));  h[is.na(h)] <- 0
  m  <- as.numeric(str_extract(x, "\\d+(?=m(?!s))")); m[is.na(m)] <- 0
  s  <- as.numeric(str_extract(x, "\\d+(?=s)"));  s[is.na(s)] <- 0
  ms <- as.numeric(str_extract(x, "\\d+(?=ms)")); ms[is.na(ms)] <- 0
  
  # Create lubridate period
  period <- hours(h) + minutes(m) + seconds(s) + milliseconds(ms)
  
  # Convert to minutes as numeric
  as.numeric(period, "seconds")
}
#######################################################
setwd("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/")
dir <- "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/benchmark/"
spn <- "SPN01"

#purities <- c(0.3,0.6,0.9)
#coverages <- c(50,100)

purity <- 0.9
coverage <- 100


get_time_pipeline <- function(spn,coverage,purity,pipeline){
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

    df_sarek_preprocess <- df_resources_sarek %>%
      filter(subworkflow3%in%selected_subworkflow3) %>% 
      group_by(subworkflow3) %>% 
      summarise(
        mean_duration = mean(as.numeric(duration_seconds, units = "secs")),
        st_duration = sd(as.numeric(duration_seconds, units = "secs"))
      ) %>% 
      mutate(process=subworkflow3,
             coverage = coverage,
             purity = purity,
             step = "preprocess"
      ) %>% 
      select(process,mean_duration,st_duration,coverage,purity,step)
    
    df_sarek_variant_calling <- df_resources_sarek %>%
      filter(module%in%selected_modules) %>% 
      group_by(module) %>% 
      summarise(
        mean_duration = mean(as.numeric(duration_seconds, units = "secs")),
        st_duration = sd(as.numeric(duration_seconds, units = "secs"))
      ) %>% 
      mutate(process=module,
             coverage = coverage,
             purity = purity,
             step = "variant_calling"
      ) %>% 
      select(process,mean_duration,st_duration,coverage,purity,step)
      
    df_sarek_sequenza <- df_resources_sequenza %>%
      group_by(module) %>% 
      summarise(
        mean_duration = mean(as.numeric(duration_seconds, units = "secs")),
        st_duration = sd(as.numeric(duration_seconds, units = "secs"))
      ) %>% 
      mutate(process=module,
             coverage = coverage,
             purity = purity,
             step = "variant_calling"
      ) %>% 
      select(process,mean_duration,st_duration,coverage,purity,step)
    df_sarek_battenberg <- df_resources_battengerg %>%
      group_by(process) %>% 
      summarise(
        mean_duration = mean(as.numeric(duration_seconds, units = "secs")),
        st_duration = sd(as.numeric(duration_seconds, units = "secs"))
      ) %>% 
      mutate(process=process,
             coverage = coverage,
             purity = purity,
             step = "variant_calling"
      ) %>% 
      select(process,mean_duration,st_duration,coverage,purity,step)

    
    
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
spns <- list("SPN01")
df_nexflow <- lapply(spns, function(x){
  df_final_150_s <- get_time_pipeline(spn = x,purity = 0.9,coverage = 150,pipeline = "sarek")
  df_final_100_s <- get_time_pipeline(spn = x,purity = 0.9,coverage = 100,pipeline = "sarek")
  df_final_50_s <- get_time_pipeline(spn = x,purity = 0.9,coverage = 50,pipeline = "sarek")
  df_final_150_t <- get_time_pipeline(spn = x,purity = 0.9,coverage = 150,pipeline = "tumourevo")
  df_final_100_t <- get_time_pipeline(spn = x,purity = 0.9,coverage = 100,pipeline = "tumourevo")
  df_final_50_t <- get_time_pipeline(spn = x,purity = 0.9,coverage = 50,pipeline = "tumourevo")
  df_final <- rbind(df_final_150_s,df_final_100_s,df_final_50_s,df_final_150_t,df_final_100_t,df_final_50_t)
  return(df_final)
}) %>% bind_rows() %>% 
  mutate(substep=case_when(
    process%in%c("ASCAT","BATTENBERG")~"cna calling",
    str_detect(process, ("SEQUENZA|CNVKIT")) ~ "cna calling",
    step=="tumourevo" & process =="SIGNATURE_DECONVOLUTION" ~ "signature deconvolution",
    step=="tumourevo" & process =="SUBCLONAL_DECONVOLUTION" ~ "subclonal deconvolution",
    step=="tumourevo" & process%in%c("DRIVER_ANNOTATION","VCF_ANNOTATE_ENSEMBLVEP") ~ "variant/driver annotation",
    step=="tumourevo" & process=="QC" ~ "quality control",
    TRUE~"snv/indel calling"))


