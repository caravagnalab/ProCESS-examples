setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')
source("../validation/SCOUT/colors.R")

out = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/')
dir.create(out, recursive = T, showWarnings = F)

base = '/orfeo/cephfs/scratch/cdslab/shared/SCOUT/validation_subclonal_new/tables/'

color_palette_process = RColorBrewer::brewer.pal(n = 8, name = "Dark2") 
names(color_palette_process) <- c('Clonal', paste0('Clone ', 1:7))
color_palette_process['Subclonal'] = 'gray70'

df_score <- tibble(SPN = paste0('SPN0', 1:7),
                   signature = c(F,F,F,F,F,T,T),
                   driver = c(T,T,T,T,T,T,T),
                   tail = c(T,T,T,T,T,T,T))


get_exposure <- function(table){
  
  nmuts <- table %>% 
    group_by(patient_id, coverage, purity, cluster_id_process) %>% 
    summarise(
      tot_nmuts = n_distinct(mutation_id),
      .groups = "drop"
    )
  
  exposure_raw <- table %>% 
    group_by(patient_id, coverage, purity, causes, cluster_id_process) %>% 
    summarise(
      nmuts_cause = n_distinct(mutation_id),
      .groups = "drop"
    )
  
  exposure <- exposure_raw %>% 
    left_join(
      nmuts,
      by = c("patient_id", "coverage", "purity", "cluster_id_process")
    ) %>% 
    mutate(
      exposure = nmuts_cause / tot_nmuts
    )
  return(exposure)
}

cosine_similarity <- function(vec1, vec2) {
  sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
}


compare_signatures <- function(df1, df2) {
  
  df1_f <- df1 %>% filter(exposure>.1)
  df2_f <- df2 %>% filter(exposure>.1)
  
  # Check if signatures are identical
  sigs_match <- setequal(df1_f$causes, df2_f$causes)
  ndiff <- length(setdiff(df1_f$causes, df2_f$causes))
  ntot <- length(c(df1_f$causes, df2_f$causes) %>% unique())
  
  comparison <- full_join(
    df1 %>% select(causes, exposure, nmuts_cause),
    df2 %>% select(causes, exposure, nmuts_cause),
    by = "causes",
    suffix = c("_tmp", "_clonal")
  ) %>%
    mutate(
      diff_exposure = exposure_tmp - exposure_clonal,
      diff_n_sig = nmuts_cause_tmp - nmuts_cause_clonal
    )%>% 
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  comparison_data <- full_join(
    df1, 
    df2, 
    by = "causes", 
    suffix = c("_tmp", "_clonal")
  ) %>% 
    select(causes, exposure_tmp, exposure_clonal, nmuts_cause_tmp, nmuts_cause_clonal) %>% 
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  # 4. Calculate similarities
  cos_sim_exposure <- cosine_similarity(
    comparison_data$exposure_tmp, 
    comparison_data$exposure_clonal
  )
  
  cos_sim_n_sig <- cosine_similarity(
    comparison_data$nmuts_cause_tmp, 
    comparison_data$nmuts_cause_clonal
  )
  
  
  return(list('df' = comparison,
              'match' = sigs_match, 
              'cs_exp' = cos_sim_exposure, 
              'cs_n' = cos_sim_n_sig,
              'n_diff' = ndiff,
              'n_tot' = ntot))
}

spn='SPN04'
cov=50
pur=.9


df_all <- tibble()
for (spn in paste0('SPN0',1:7)){

  true_drivers_table = readRDS(file.path("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/drivers", spn, "process_drivers.rds")) %>% as_tibble()

  true_drivers_table = true_drivers_table  %>% rowwise() %>%
    mutate(code=replace(code, is.na(code), paste0(type, '_', CNA_type, '_', chr, '_', start, '_', end)),
           mutation_id=paste0(spn, ":", chr, ":", start,  ":", alt)) %>%
    ungroup() %>%
    select(mutation_id, code)

  print(spn)
  for (cov in c(50, 100, 150)){
    for (pur in c(0.9, 0.6, 0.3)){

      table_multi <- readRDS(paste0(base, 'table_process_', spn, '_', cov, 'x_', pur, 'p_mutect2_ascat.rds')) %>%
        filter(!str_detect(causes, "errors"))

      table_multi <- table_multi %>%
        left_join(true_drivers_table, by = 'mutation_id') %>%
        select(-driver_label_process) %>%
        mutate(driver_label_process=code)


      if(spn=='SPN03'){
        table_multi = table_multi %>% mutate(mutation_id = if_else(
          mutation_id=="SPN03:9:136496197:",
          'SPN03:9:136496196:C',
          mutation_id
        ))
      }


      n_samples = length(table_multi$sample_id %>% unique())

      c = table_multi %>%
        distinct(cluster_id_process, sample_id) %>%   # keep unique pairs
        dplyr::count(cluster_id_process) %>%
        filter(n==n_samples) %>%
        select(cluster_id_process)

      table_multi = table_multi %>%
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


      table_uni <- readRDS(paste0(base, 'table_process_univariate_w_private_', spn, '_', cov, 'x_', pur, 'p_mutect2_ascat.rds')) %>%
        filter(!str_detect(causes, "errors"))


      table_multi_plot <- table_multi %>% select(sample_id,
                                                 mutation_id,
                                                 is_driver_process,
                                                 cluster_id_process,
                                                 vaf_process,
                                                 driver_label_process)

      if (spn == 'SPN07'){
        dups <- table_multi_plot %>%
          ungroup() %>%
          dplyr::summarise(n = dplyr::n(), .by = c(mutation_id, is_driver_process, cluster_id_process,
                                                   driver_label_process, sample_id)) %>%
          dplyr::filter(n > 1L)
        table_multi_plot <- table_multi_plot %>% filter(!mutation_id%in%dups$mutation_id)
      }

      table_wide <- table_multi_plot %>%
        pivot_wider(
          names_from  = sample_id,
          values_from = vaf_process,
          values_fill = 0
        )


      table_uni <- table_uni %>%
        group_by(cluster_id_process, sample_id) %>%
        mutate(is_clonal_process=replace(FALSE, ccf_process > 0.9, TRUE)) %>%
        ungroup() %>%
        mutate(cluster_id_process_full = cluster_id_process) %>%
        mutate(cluster_id_process = replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal')) %>%
        separate(cluster_id_process, sep = '_', into = c('cluster_id_process', 'tmp')) %>%
        select(-tmp)

      # get driver
      driver <- table_multi %>%
        filter(is_driver_process==T) %>%
        select(cluster_id_process, driver_label_process) %>%
        distinct() %>%
        mutate(contains_driver = T)

      no_driver <- table_multi %>%
        filter(is_driver_process==F) %>%
        select(cluster_id_process) %>%
        distinct() %>%
        filter(!(cluster_id_process %in% driver$cluster_id_process)) %>%
        mutate(contains_driver = F)

      df_driver <- bind_rows(driver, no_driver)

      # get tail
      df_tail <- table_multi %>%
        left_join(table_uni %>% select(mutation_id, cluster_id_process) %>%
                    dplyr::rename(cluster_id_process_uni = cluster_id_process)) %>%
        group_by(mutation_id) %>%
        mutate(never_tail = all(cluster_id_process_uni != "Subclonal")) %>%
        group_by(cluster_id_process) %>%
        summarise(n_never_tail = sum(never_tail, na.rm = T)/n())

      # get exposure
      n_muts <- table_multi %>%
        group_by(cluster_id_process) %>%
        summarise(n = n())

      table <- table_multi #%>% filter(cluster_id_process %in% n_muts$cluster_id_process)

      table_sbs <- table %>%
        filter(str_detect(causes, "SBS"))

      table_id <- table %>%
        filter(str_detect(causes, "ID"))

      exposure_sbs <- get_exposure(table_sbs)
      exposure_id <- get_exposure(table_id)

      clonal_sig_sbs = exposure_sbs %>% filter(cluster_id_process == 'Clonal')
      clonal_sig_id = exposure_id %>% filter(cluster_id_process == 'Clonal')


      background_sbs = exposure_sbs %>%
        filter(cluster_id_process != 'Clonal') %>%
        group_by(causes) %>%
        summarise(nmuts_cause = sum(nmuts_cause)) %>%
        mutate(tot_nmuts = sum(nmuts_cause)) %>%
        mutate(exposure = nmuts_cause/tot_nmuts)

      background_id = exposure_id %>%
        filter(cluster_id_process != 'Clonal') %>%
        group_by(causes) %>%
        summarise(nmuts_cause = sum(nmuts_cause)) %>%
        mutate(tot_nmuts = sum(nmuts_cause)) %>%
        mutate(exposure = nmuts_cause/tot_nmuts)


      table_sbs <- lapply(unique(exposure_sbs$cluster_id_process), FUN = function(c){
        if (c != 'Clonal'){
          tmp <- exposure_sbs %>% filter(cluster_id_process == c)
          result_table <- compare_signatures(df1 = tmp, df2 = clonal_sig_sbs)
          result <- result_table$df %>% mutate(match = result_table$match,
                                               cs_exp = result_table$cs_exp,
                                               cs_nsig = result_table$cs_n,
                                               n_diff = result_table$n_diff,
                                               n_tot = result_table$n_tot) %>%
            mutate(cluster = as.character(c))
          return(result)
        }
      }) %>% bind_rows()


      table_sbs_background <- lapply(unique(exposure_sbs$cluster_id_process), FUN = function(c){
        if (c != 'Clonal'){
          tmp <- exposure_sbs %>% filter(cluster_id_process == c)
          result_table <- compare_signatures(df1 = tmp, df2 = background_sbs)
          result <- result_table$df %>% mutate(match = result_table$match,
                                               cs_exp = result_table$cs_exp,
                                               cs_nsig = result_table$cs_n,
                                               n_diff = result_table$n_diff,
                                               n_tot = result_table$n_tot) %>%
            mutate(cluster = as.character(c))
          return(result)
        }
      }) %>% bind_rows()


      table_id <- lapply(unique(exposure_id$cluster_id_process), FUN = function(c){
        if (c != 'Clonal'){
          tmp <- exposure_id %>% filter(cluster_id_process == c)
          result_table <- compare_signatures(df1 = tmp, df2 = clonal_sig_id)
          result <- result_table$df %>% mutate(match = result_table$match,
                                               cs_exp = result_table$cs_exp,
                                               cs_nsig = result_table$cs_n,
                                               n_diff = result_table$n_diff,
                                               n_tot = result_table$n_tot) %>%
            mutate(cluster = as.character(c))
          return(result)
        }
      }) %>% bind_rows()

      table_id_background <- lapply(unique(exposure_sbs$cluster_id_process), FUN = function(c){
        if (c != 'Clonal'){
          tmp <- exposure_id %>% filter(cluster_id_process == c)
          result_table <- compare_signatures(df1 = tmp, df2 = background_id)
          result <- result_table$df %>% mutate(match = result_table$match,
                                               cs_exp = result_table$cs_exp,
                                               cs_nsig = result_table$cs_n,
                                               n_diff = result_table$n_diff,
                                               n_tot = result_table$n_tot) %>%
            mutate(cluster = as.character(c))
          return(result)
        }
      }) %>% bind_rows()

      df_signature <- table_sbs %>%
        select(cluster, match, cs_exp, cs_nsig, n_diff, n_tot) %>%
        mutate(type = 'SBS') %>%
        distinct() %>%
        bind_rows(table_id %>%
                    select(cluster, match, cs_exp, cs_nsig, n_diff, n_tot) %>%
                    distinct() %>%
                    mutate(type = 'ID')) %>%
        mutate(n_diff_rel = n_diff/n_tot)

      df_signature_background <-  table_sbs_background %>%
        select(cluster, match, cs_exp, cs_nsig, n_diff, n_tot) %>%
        mutate(type = 'SBS') %>%
        distinct() %>%
        bind_rows(table_id_background %>%
                    select(cluster, match, cs_exp, cs_nsig, n_diff, n_tot) %>%
                    distinct() %>%
                    mutate(type = 'ID')) %>%
        mutate(n_diff_rel = n_diff/n_tot) %>%
        group_by(across(c(-type, -cs_nsig, -cs_exp, -match, -n_diff_rel, -n_diff, -n_tot))) %>%
        summarize(match = sum(match),
                  cs_exp = mean(cs_exp),
                  n_diff_rel = mean(n_diff_rel)) %>%
        dplyr::rename(bg_match = match, bg_cs = cs_exp, bg_diff = n_diff_rel)

      table <- full_join(df_tail, df_driver) %>%
        dplyr::rename(cluster = cluster_id_process) %>%
        full_join(df_signature) %>%
        full_join(df_signature_background)

      score <- df_score %>%
        filter(SPN == spn) %>%
        mutate(signature = as.numeric(signature),
               driver = as.numeric(driver),
               tail = as.numeric(tail))
      n <- sum(score$signature, score$driver, score$tail)

      interpret_table <- table %>%
        select(-match) %>%
        group_by(across(c(-type, -cs_nsig, -cs_exp, -n_diff_rel, -n_diff, -n_tot))) %>%
        summarize(cs_exp = min(cs_exp),
                  n_diff_rel = max(n_diff_rel),
                  bg_cs = min(bg_cs)) %>%
        ungroup() %>%
        mutate(driver = ifelse(contains_driver == F, 0, 1),
               cs_sign = 1-cs_exp,
               bg_cs = 1-bg_cs,
               cs_sign = ifelse(is.na(cs_exp) & cluster == 'Clonal', 1, cs_sign),
               bg_cs = ifelse(is.na(bg_cs) & cluster == 'Clonal', 1, bg_cs),
               n_diff_rel = ifelse(is.na(cs_exp) & cluster == 'Clonal', 1, n_diff_rel)) %>%
        mutate(score_spn = (score$driver * driver + score$tail * n_never_tail + score$signature* ((cs_sign+n_diff_rel+bg_cs)/3))/n,
               score_all = (driver + ((cs_sign+n_diff_rel+bg_cs)/3) + n_never_tail)/3,
               score_driver = driver,
               score_sign = (cs_sign+n_diff_rel+bg_cs)/3,
               score_tail = n_never_tail,
               score_no_tail = (driver + ((cs_sign+n_diff_rel+bg_cs)/3))/2,
               score_no_driver = (n_never_tail + ((cs_sign+n_diff_rel+bg_cs)/3))/2,
               score_no_sign = (driver + n_never_tail)/2) %>%
        mutate(coverage = cov, purity = pur, spn = spn) %>%
        left_join(score %>% dplyr::rename(spn = SPN, driver_v = driver))

      df_all <- bind_rows(df_all, interpret_table)
    }
  }
}
saveRDS(object = df_all, file = paste0(out, 'process.rds'))


color_palette_process = c(RColorBrewer::brewer.pal(n = 8, name = "Dark2"),
                          RColorBrewer::brewer.pal(n = 9, name = "Set1"),
                          RColorBrewer::brewer.pal(n = 8, name = "Set2"))
names(color_palette_process) <- df_all$cluster %>% unique()
color_palette_process['Subclonal'] = 'gray70'

df_all <- readRDS(paste0(out, 'process.rds'))
plt <- df_all %>%
  pivot_longer(cols = c(score_driver, score_all, score_tail, score_no_driver, score_no_tail, score_no_sign, score_sign)) %>%
  mutate(name = factor(name, levels = c('score_driver', 'score_tail', 'score_sign','score_no_driver', 'score_no_tail', 'score_no_sign', 'score_all'))) %>%
  ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.9, ymax = 1,
           fill = "palegreen4", alpha = 0.2) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.55, ymax = .9,
           fill = "goldenrod", alpha = 0.2) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.2, ymax = 0.55,
           fill = "salmon1", alpha = 0.2) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.2, ymax = 0,
           fill = "gainsboro", alpha = 0.2) +
  stat_summary(
    aes(x = name, y = value, color = cluster),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3),
    size = .4,
    show.legend = T
  ) +
  stat_summary(
    aes(x = name, y = value, color = cluster, group=cluster),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3),
    geom = 'line',
    linewidth=.6,
    show.legend = F
  ) +
  #geom_text(data = ~ filter(.x, contains_driver & name == 'score_all'), aes(x = name, y = value+0.03, col = cluster, group = cluster, label = driver_label_process)) +
  theme_minimal() +
  ylab('Score') +
  xlab('')+
  scale_color_manual('Cluster',
                     values = color_palette_process) +
  scale_x_discrete(labels = c('score_driver' = 'Driver',
                              'score_tail'   = 'Tail',
                              'score_sign'   = 'Signature',
                              'score_no_driver' = 'Signature\nTail',
                              'score_no_tail' = 'Driver\nSignature',
                              'score_no_sign' = 'Driver\nTail',
                              'score_all'    = 'All')) +
  facet_wrap(.~spn)

ggsave(plot = plt,
       filename = paste0(out, "/plot_process/process.png"),
       width = 12, height = 9, units = 'in', dpi = 200)
