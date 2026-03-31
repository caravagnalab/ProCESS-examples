library(dplyr)
library(ggplot2)
source('../../getters/tumourevo_getters.R')
source('../../getters/process_getters.R')
setwd("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/drivers/")

coverage_list = c(50, 100, 150)
purity_list = c(0.9, 0.6, 0.3)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04','SPN05', 'SPN06', 'SPN07')
combs = expand.grid(spn = spn_list,
                    coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller = cna_caller_list)


i=3
results <- lapply(1:nrow(combs), FUN = function(i){
  print(i)
  sp = combs[i, "spn"]
  cov = combs[i, "coverage"]
  pur = combs[i, "purity"]
  vcf_call = combs[i, 'vcf_caller']
  cna_call = combs[i, 'cna_caller']

  
  true_drivers = readRDS(file.path("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/drivers", sp, "process_drivers.rds")) %>% 
    as_tibble() %>% 
    filter(type == 'SID') %>% 
    select(chr, start, end, ref, alt, code)
  
  samples = get_sample_names(sp)
  
  inf_driver = lapply(samples, FUN = function(s){
    tmp_sample = paste(sp, s, sep = '_')
    # file = get_tumourevo_driver(spn = sp, 
    #                             coverage = cov, 
    #                             purity = pur,  
    #                             vcf_caller = vcf_call, 
    #                             cna_caller = cna_call, 
    #                             sample = s)
    if (sp == 'SPN02'){
      file = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/',sp,'/tumourevo/',cov, 'x_', pur, 'p_', vcf_call,'_',cna_call,
                    '/driver_annotation/annotate_driver/SCOUT/',sp,'/',tmp_sample,'/',paste('SCOUT', sp, tmp_sample, 'driver.rds', sep = '_'))
    } else {cna_call
      file = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/',sp,'/tumourevo/',cov, 'x_', pur, 'p_', vcf_call,'_',cna_call,
                    '/driver_annotation/annotate_driver/SCOUT/',sp,'/',tmp_sample,'/',paste('SCOUT', sp, tmp_sample, 'driver.rds', sep = '_'))
    }
    if (length(file) > 0 & file.exists(file)){
      readRDS(file)[[tmp_sample]][['mutations']] %>% 
        filter(is_driver == T) %>% 
        select(chr, from, to, ref, alt, SYMBOL, driver_label)
    }
    
  }) %>% bind_rows() %>% distinct()

  if (nrow(inf_driver) > 0){
    true <- true_drivers %>%
      transmute(chr = paste0("chr", chr),
                pos = start,
                ref, alt)
    
    pred <- inf_driver %>%
      transmute(chr,
                pos = from,
                ref, alt)
    
    if (sp == 'SPN03'){
      pred <- pred %>% 
        mutate(ref = ifelse(pos == 136496196, 'AG', ref)) %>% 
        mutate(alt = ifelse(pos == 136496196, '', alt)) %>% 
        mutate(pos = ifelse(pos == 136496196, 136496197, pos))
    }
    
    TP <- inner_join(true, pred, by = c("chr", "pos", "ref", "alt")) %>% nrow()
    FP <- anti_join(pred, true, by = c("chr", "pos", "ref", "alt")) %>% nrow()
    FN <- anti_join(true, pred, by = c("chr", "pos", "ref", "alt")) %>% nrow()
    
    Precision <- TP / (TP + FP)
    Recall <- TP / (TP + FN)
    F1 <- 2 * Precision * Recall / (Precision + Recall)
    
    tibble(spn = sp, coverage = cov, purity = pur,
           vcf_caller = vcf_call, cna_caller = cna_call,
           TP=TP, FP=FP, FN=FN,
           Precision=Precision,Recall=Recall,F1=F1)
  }
}) %>% bind_rows()


saveRDS(results, file = 'performance_driver.rds')
