setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

indir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/assing_signature/"

# option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN03'),
#                     make_option(c("--cna_caller"), type = "character", default = 'ascat'),
#                     make_option(c("--vcf_caller"), type = "character", default = 'mutect2')
# )
# 
# opt_parser <- OptionParser(option_list = option_list)
# opt <- parse_args(opt_parser)
# 
# spn = opt$spn_id
# cna_caller = opt$cna_caller
# mut_caller = opt$vcf_caller


cosine_similarity <- function(vec1, vec2) {
  sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
}

compare_signatures <- function(df1, df2) {
  
  # Check if signatures are identical
  sigs_match <- setequal(df1$Signature, df2$Signature)
  
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
    select(Signature, Exposure_tmp, Exposure_clonal,n_sig_tmp, n_sig_clonal) %>% 
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  # 4. Calculate similarities
  cos_sim_exposure <- cosine_similarity(
    comparison_data$Exposure_tmp, 
    comparison_data$Exposure_clonal
  )
  
  cos_sim_n_sig <- cosine_similarity(
    comparison_data$n_sig_tmp, 
    comparison_data$n_sig_clonal
  )
  
  
  return(list('df' = comparison,
              'match' = sigs_match, 
              'cs_exp' = cos_sim_exposure, 
              'cs_n' = cos_sim_n_sig))
}



cov = 50
pur = 0.3
tool = 'pyclonevi'
signature_tool = 'BASCULE'
sign_type='ID'

mut_caller= 'mutect2'
cna_caller = 'ascat'


summary_table <- tibble()
final_table <- tibble()
for (spn in paste0('SPN0', 1:7)){
  print(spn)
  for (cov in c(50, 100, 150)){
    for (pur in c(0.3, 0.6, 0.9)){
      for (tool in c('viber', 'pyclonevi')){
        for (signature_tool in c('BASCULE', 'SigProfiler')){
          for (sign_type in c('SBS', 'ID')){
            #print(c(spn, cov, pur, tool, signature_tool, sign_type))
            if (sign_type == 'SBS'){
              s_type = 'SBS96'
              colors = sbs_colors
            } else {
              s_type = 'ID83'
              colors = id_colors
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
            
              mutations <-  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal/tables_interpreted/', tool, '_', spn, '_', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '.rds'))
              clonal <- mutations %>% filter(is_clonal_tool == T) %>% pull(cluster_id_tool) %>% unique()
              if (tool == 'pyclonevi' & signature_tool == 'BASCULE'){
                clonal <- paste0('X', clonal)
              }
              
              df <- df %>% 
                mutate(is_clonal = ifelse(Samples==clonal, T, F)) %>% 
                left_join(counts) %>% 
                mutate(n_sig = round(n*Exposure,0))
              clonal_sig <- df %>% filter(is_clonal == T) %>% filter(n_sig > 0)
              other_sig <- df %>% filter(is_clonal == F) %>% filter(n_sig > 0)
              
              tmp_final_table <- lapply(unique(other_sig$Samples), FUN = function(c){
                tmp <- other_sig %>% filter(Samples == c)
                if (unique(tmp$n) >= 50){ 
                  result_table <- compare_signatures(df1 = tmp, df2 = clonal_sig) 
                  result <- result_table$df %>% mutate(match = result_table$match,
                                                       cs_exp = result_table$cs_exp,
                                                       cs_nsig = result_table$cs_n) %>% 
                    mutate(cluster = as.character(c))
                  return(result)
                }
              }) %>% bind_rows()
              if (nrow(tmp_final_table) > 0){
                
                final_table <- final_table %>% 
                  bind_rows(tmp_final_table %>%
                  mutate(sig = sign_type, coverage = cov, purity = pur, dec_tool = tool, sig_tool = signature_tool,  spn = spn))
                
                summary_table <- summary_table %>% 
                  bind_rows(tmp_final_table %>% 
                  select(cluster, match, cs_exp, cs_nsig) %>% 
                    distinct() %>% 
                    mutate(sig = sign_type, coverage = cov, purity = pur, dec_tool = tool, sig_tool = signature_tool, spn = spn)) 
              }
            }
          } 
        }
      }
    }
  }
}


saveRDS(summary_table, '/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/signature_deconvolution_summary.rds')
saveRDS(final_table, '/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/signature_deconvolution_full.rds')
