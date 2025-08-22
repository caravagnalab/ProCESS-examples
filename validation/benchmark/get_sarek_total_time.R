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


get_time_pipeline <- function(spn,coverage,purity){
  comb <- paste0(coverage,"x_",purity,"p")
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
      realtime_seconds = parse_with_lubridate(realtime)
    ) %>% 
    group_by(process) %>% 
    mutate(duration_seconds_sum=sum(duration_seconds),
           realtime_seconds_sum=sum(realtime_seconds)) 
  df_resources_normal <- df_final_normal %>% 
    dplyr::filter(status=="COMPLETED") %>% 
    tidyr::separate(col = "name", into = c("process","tag"),sep = " ") %>% 
    mutate(
      duration_seconds = parse_with_lubridate(duration),
      realtime_seconds = parse_with_lubridate(realtime)
    ) %>% 
    group_by(process) %>% 
    mutate(duration_seconds_sum=sum(duration_seconds),
           realtime_seconds_sum=sum(realtime_seconds)) 
  
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
  
  df_sub <- data.frame(
    process = character(),
    mean_duration = period(),
    st_duration = period(),
    coverage = numeric(),
    purity = numeric(),
    n_jobs =numeric(),
    n_samples=numeric(),
    step = character(),
    stringsAsFactors = FALSE
  )
  
  df_mod <- data.frame(
    process = character(),
    mean_duration = period(),
    st_duration = period(),
    coverage = numeric(),
    purity = numeric(),
    n_jobs =numeric(),
    n_samples=numeric(),
    step = character(),
    stringsAsFactors = FALSE
  )
  
  df_seq <- data.frame(
    process = character(),
    mean_duration = period(),
    st_duration = period(),
    coverage = numeric(),
    purity = numeric(),
    n_jobs =numeric(),
    n_samples=numeric(),
    step = character(),
    stringsAsFactors = FALSE
  )
  
  df_batt <- data.frame(
    process = character(),
    mean_duration = period(),
    st_duration = period(),
    coverage = numeric(),
    purity = numeric(),
    n_jobs =numeric(),
    n_samples=numeric(),
    step = character(),
    stringsAsFactors = FALSE
  )
  
  for (proc in selected_subworkflow3){
    tmp <- df_resources_sarek %>%
      filter(subworkflow3 == proc)
    tot_jobs <- nrow(tmp)
    tot_samples <- length(unique(tmp$tag))
    tmp <- tmp %>% 
      mutate(
        # Parse the duration string into a Period
        duration_parsed = parse_with_lubridate(duration)
      ) %>%
      summarise(
        mean_duration = mean(as.numeric(duration_parsed, units = "secs")),
        st_duration = sd(as.numeric(duration_parsed, units = "secs"))
      ) %>%
      mutate(
        mean_duration = round(seconds_to_period(mean_duration),0),
        st_duration = seconds_to_period(st_duration)
      ) %>% 
      mutate(process=proc,
             coverage = coverage,
             purity = purity,
             n_jobs=tot_jobs,
             n_samples=tot_samples,
             step = "preprocess"
      )
    df_sub <- rbind(df_sub, tmp)
  }
  
  
  for (proc in selected_modules){
    tmp <- df_resources_sarek %>%
      filter(module == proc)
    tot_jobs <- nrow(tmp)
    tot_samples <- length(unique(tmp$tag))
    tmp <- tmp %>% 
      mutate(
        # Parse the duration string into a Period
        duration_parsed = parse_with_lubridate(duration)
      ) %>%
      summarise(
        mean_duration = mean(as.numeric(duration_parsed, units = "secs")),
        st_duration = sd(as.numeric(duration_parsed, units = "secs"))
      ) %>%
      mutate(
        mean_duration = round(seconds_to_period(mean_duration),0),
        st_duration = seconds_to_period(st_duration),
      ) %>% 
      mutate(process=proc,
             coverage = coverage,
             purity = purity,
             n_jobs=tot_jobs,
             n_samples=tot_samples,
             step = "variant_calling"
      )
    df_mod <- rbind(df_mod, tmp)
  }
  
  
  for (proc in unique(df_resources_sequenza$module)){
    tmp <- df_resources_sequenza %>%
      filter(module == proc)
    tot_jobs <- nrow(tmp)
    tot_samples <- length(unique(tmp$tag))
    tmp <- tmp %>% 
      mutate(
        # Parse the duration string into a Period
        duration_parsed = parse_with_lubridate(duration)
      ) %>%
      summarise(
        mean_duration = mean(as.numeric(duration_parsed, units = "secs")),
        st_duration = sd(as.numeric(duration_parsed, units = "secs"))
      ) %>%
      mutate(
        mean_duration = round(seconds_to_period(mean_duration),0),
        st_duration = seconds_to_period(st_duration),
      ) %>% 
      mutate(process=proc,
             coverage = coverage,
             purity = purity,
             n_jobs=tot_jobs,
             n_samples=tot_samples,
             step = "variant_calling"
      )
    df_seq <- rbind(df_seq, tmp)
  }
  
  
  for (proc in unique(df_resources_battenbgerg$process)){
    tmp <- df_resources_battenbgerg %>%
      filter(process == proc)
    tot_jobs <- nrow(tmp)
    tot_samples <- length(unique(tmp$tag))
    tmp <- tmp %>% 
      mutate(
        # Parse the duration string into a Period
        duration_parsed = parse_with_lubridate(duration)
      ) %>%
      summarise(
        mean_duration = mean(as.numeric(duration_parsed, units = "secs")),
        st_duration = sd(as.numeric(duration_parsed, units = "secs"))
      ) %>%
      mutate(
        mean_duration = round(seconds_to_period(mean_duration),0),
        st_duration = seconds_to_period(st_duration),
      ) %>% 
      mutate(process=proc,
             coverage = coverage,
             purity = purity,
             n_jobs=tot_jobs,
             n_samples=tot_samples,
             step = "variant_calling"
      )
    df_batt <- rbind(df_batt, tmp)
  }
  
  df_final <- rbind(df_mod,df_sub,df_batt,df_seq)
  df_final <- df_final %>% 
    dplyr::mutate(SPN=spn)
  return(df_final)
}

spns <- list("SPN01","SPN02","SPN03","SPN04","SPN06","SPN07")
df_tmp <- lapply(spns, function(x){
  df_final_100 <- get_time_pipeline(spn = x,purity = 0.9,coverage = 100)
  df_final_50 <- get_time_pipeline(spn = x,purity = 0.9,coverage = 50)
  df_final <- rbind(df_final_100,df_final_50)
  return(df_final)
}) %>% bind_rows()


df_plot <- df_tmp %>%
  mutate(
    mean_seconds = as.numeric(mean_duration, units = "secs"),
    sd_seconds   = as.numeric(st_duration, units = "secs"),
    mean_label   = as.character(mean_duration)   # keep the original period as string
  )

plot_time_all <- ggplot(df_plot, aes(x = reorder(process, mean_seconds), 
                    y = mean_seconds,
                    fill = factor(step))) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_seconds - sd_seconds,
                    ymax = mean_seconds + sd_seconds),
                width = 0.2, color = "black") +
  geom_text(aes(label = mean_label),
            hjust = -0.2,vjust = 1.5,
            size = 3) +
  coord_flip() +
  facet_grid(coverage~SPN)+
  labs(x = "Process", 
       y = "Mean duration (seconds)",
       fill = "# parallel jobs",
       title = "Mean ± SD of Workflow Step Durations")
ggsave(filename = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/benchmark/plot_time_sarek.png",plot = plot_time_all)
