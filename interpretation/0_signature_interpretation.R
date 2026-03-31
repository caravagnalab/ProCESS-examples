setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(tidyverse)
source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')

indir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/assign_signature/"

mut_caller = "mutect2"
cna_caller = "sequenza"
out_path = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/interpretation_', mut_caller, "_", cna_caller, '/')
dir.create(out_path, recursive = T, showWarnings = F)


cosine_similarity <- function(vec1, vec2) {
  sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
}

compare_signatures <- function(df1, df2) {
  
  df1_f <- df1 %>% filter(Exposure>.1)
  df2_f <- df2 %>% filter(Exposure>.1)
  
  # Check if signatures are identical
  sigs_match <- setequal(df1_f$Signature, df2_f$Signature)
  ndiff <- length(setdiff(df1_f$Signature, df2_f$Signature))
  ntot <- length(c(df1_f$Signature, df2_f$Signature) %>% unique())
  
  comparison <- full_join(
    df1 %>% select(Signature, Exposure, n_sig),
    df2 %>% select(Signature, Exposure, n_sig),
    by = "Signature",
    suffix = c("_tmp", "_clonal")
  ) %>%
    mutate(
      diff_exposure = Exposure_tmp - Exposure_clonal,
      diff_n_sig = n_sig_tmp - n_sig_clonal
    )%>% 
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  comparison_data <- full_join(
    df1, 
    df2, 
    by = "Signature", 
    suffix = c("_tmp", "_clonal")
  ) %>% 
    select(Signature, Exposure_tmp, Exposure_clonal, n_sig_tmp, n_sig_clonal) %>% 
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  # 4. Calculate similarities
  cos_sim_exposure <- cosine_similarity(
    comparison_data$Exposure_tmp, 
    comparison_data$Exposure_clonal
  )
  
  
  
  return(list('df' = comparison,
              'match' = sigs_match, 
              'cs_exp' = cos_sim_exposure, 
              'n_diff' = ndiff,
              'n_tot' = ntot))
}


spn = 'SPN02'
cov = 100
pur = 0.9
tool = 'pyclonevi'
signature_tool = 'BASCULE'
sign_type='SBS'


spn_list = paste0('SPN0', c(1,2,3,4,5,6,7))
summary_table <- tibble()
final_table <- tibble()
for (spn in spn_list){
  print(spn)
  for (cov in c(50, 100, 150)){
    for (pur in c(0.3, 0.6, 0.9)){
      for (tool in c('viber', 'pyclonevi')){
        for (signature_tool in c('BASCULE', 'SigProfiler')){
          for (sign_type in c('SBS', 'ID')){
            #print(c(spn, cov, pur, tool, signature_tool, sign_type))
            if (sign_type == 'SBS'){
              s_type = 'SBS96'
            } else {
              s_type = 'ID83'
            }
            
            type='raw'
            
            if (signature_tool == 'BASCULE'){
              name_file = paste0('bascule_fit_',type,'.rds')
            } else{
              name_file = paste0('Assignment_Solution/Activities/Assignment_Solution_Activities.txt')
            }
            
            
            if (signature_tool == 'BASCULE'){
              if (file.exists(paste0(indir, spn, '/', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '/', tool, '_', signature_tool, '_', type, '/', s_type, '/',name_file))){
                data = readRDS(paste0(indir, spn, '/', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '/', tool, '_', signature_tool, '_', type, '/', s_type, '/',name_file))
                
                df = data[['nmf']][[sign_type]][['exposure']] %>% 
                  dplyr::rename(Signature=sigs,
                                Exposure=value,
                                Samples=samples)
                
                counts = data[['input']][[sign_type]][['counts']] %>% 
                  group_by(samples) %>% 
                  summarise(n = sum(value)) %>% 
                  dplyr::rename(Samples = samples)
              } else {
                df = NULL
              }
              
            } else {
              if (file.exists(paste0(indir, spn, '/', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '/', tool, '_', signature_tool, '_', type, '/', s_type, '/',name_file))){
                data = read.table(file = paste0(indir, spn, '/', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '/', tool, '_', signature_tool, '_', type, '/', s_type, '/',name_file),header = T) 
                df = data %>% pivot_longer(
                  cols = starts_with(sign_type),
                  names_to = "Signature",
                  values_to = "Exposure")
                counts = data %>%  mutate(n = rowSums(across(-Samples))) %>% select(Samples, n)
              } else {
                df = NULL
              }
            }
            
            if (!is.null(df)){
            
              mutations <-  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/validation_subclonal_new/tables_interpreted/', tool, '_', spn, '_', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '.rds'))
              clonal <- mutations %>% filter(is_clonal_tool == T) %>% pull(cluster_id_tool) %>% unique()
              if (length(clonal) > 0){
              
                if (tool == 'pyclonevi' & signature_tool == 'BASCULE'){
                  clonal <- paste0('X', clonal)
                }
                
                df <- df %>% 
                  mutate(is_clonal = ifelse(Samples==clonal, T, F)) %>% 
                  left_join(counts) %>% 
                  mutate(n_sig = round(n*Exposure,0))
                clonal_sig <- df %>% filter(is_clonal == T) %>% filter(n_sig > 0)
                other_sig <- df %>% filter(is_clonal == F) %>% filter(n_sig > 0)
                
                background = other_sig %>% 
                  group_by(Signature) %>% 
                  summarise(n_sig = sum(n_sig)) %>% 
                  mutate(n = sum(n_sig)) %>% 
                  mutate(Exposure = n_sig/n)
                
                
                tmp_final_table <- lapply(unique(other_sig$Samples), FUN = function(c){
                  tmp <- other_sig %>% filter(Samples == c)
                  result_table <- compare_signatures(df1 = tmp, df2 = clonal_sig) 
                  result <- result_table$df %>% mutate(match = result_table$match,
                                                       cs_exp = result_table$cs_exp,
                                                       cs_nsig = result_table$cs_n,
                                                       n_diff = result_table$n_diff,
                                                       n_tot = result_table$n_tot) %>% 
                    mutate(cluster = as.character(c))
                  return(result)
                }) %>% bind_rows()
                
                
                tmp_final_table_background <- lapply(unique(other_sig$Samples), FUN = function(c){
                  tmp <- other_sig %>% filter(Samples == c)
                  result_table <- compare_signatures(df1 = tmp, df2 = background) 
                  result <- result_table$df %>% mutate(match = result_table$match,
                                                       cs_exp = result_table$cs_exp,
                                                       cs_nsig = result_table$cs_n,
                                                       n_diff = result_table$n_diff,
                                                       n_tot = result_table$n_tot) %>% 
                    mutate(cluster = as.character(c))
                  return(result)
                }) %>% bind_rows()
                
                
                if (nrow(tmp_final_table) > 0){
                  
                  f_background <- tmp_final_table_background %>% 
                    mutate(n_rel = n_diff/n_tot) %>% 
                    select(cluster, match, cs_exp, n_rel) %>% 
                    distinct() %>% 
                    dplyr::rename(bg_match = match, bg_cs_exp =cs_exp, bg_n_rel = n_rel)
                  
                  f_table <- tmp_final_table %>% 
                    mutate(n_rel = n_diff/n_tot) %>% 
                    select(cluster, match, cs_exp, n_rel) %>% 
                    distinct() %>% 
                    left_join(f_background) %>% 
                    mutate(sig = sign_type, coverage = cov, purity = pur, dec_tool = tool, sig_tool = signature_tool, spn = spn)
                  
                  summary_table <- summary_table %>% 
                    bind_rows(f_table)
                }
              }
            }
          } 
        }
      }
    }
  }
}
saveRDS(summary_table, paste0(out_path,'signature_deconvolution_summary.rds'))
