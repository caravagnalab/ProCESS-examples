rm(list = ls())
options(bitmapType='cairo')
library(optparse)
library(tidyverse)
library(gt)

source('utils.R')

source('../../getters/process_getters.R')
source('../../getters/tumourevo_getters.R')

option_list <- list( 
  make_option(c("-v", "--variantcaller"), type="character", default='mutect2', help="variantcaller"),
  make_option(c("-c", "--cnacaller"), type="character", default='ascat', help="cna caller")
)

param <- parse_args(OptionParser(option_list=option_list))

SPN <- paste0('SPN0', c(7,1,2,3,4,6))
COV <- c(50, 100)
PUR <- c(0.3, 0.6, 0.9)

full_table <- tibble()
validate_table <- tibble()

for (spn in SPN){
  print(spn)
  
  samples <- get_sample_names(spn)
  for (cov in COV){
    for (pur in PUR){
      print(paste0('Running for: ', spn, '-', cov, 'x-', pur, 'p'))
      
      cnaqc_obj <- lapply(samples, FUN = function(sample){
        file <- get_tumourevo_qc(spn = spn, coverage = cov, purity = pur, tool = 'CNAqc', vcf_caller = param$variantcaller, cna_caller = param$cnacaller, sample = sample)
        if (length(file) > 0){
          if (file.exists(file$qc_rds)){
          readRDS(file$qc_rds)
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
          mutate(true_purity = pur, delta_purity = abs(true_purity-cn_purity)) 
        
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
        
        stats_table <- get_statistics_qc(cnaqc_obj, purity = pur)
        table <- stats_table$table %>% tibble() %>% select(-gg) %>% filter(info != 'sample') %>% tidyr::pivot_wider(values_from = values, names_from = info, names_repair = 'minimal')
        full_table <- bind_rows(full_table, table %>% mutate(spn = spn, coverage = cov, vcf_caller = param$variantcaller, cna_caller = param$cnacaller))
        validate_table <- bind_rows(validate_table, validate %>% mutate(spn = spn, coverage = cov, purity = pur, vcf_caller = param$variantcaller, cna_caller = param$cnacaller))
      }
    }
  }
}

saveRDS(full_table, file = 'table_CNAqc.rds')
saveRDS(validate_table, file = 'plot_CNAqc.rds')

validate_table_filt <- validate_table %>%
  mutate(CNAqc = ifelse(is.na(CNAqc), 'NA', CNAqc)) %>% 
  group_by(sample, coverage, purity, true_karyo) %>%
  filter(!(all(CNAqc == "Subclonal", na.rm = T))) %>%
  ungroup()

plt <- validate_table_filt %>% 
  group_by(true_karyo, CNAqc, sample, coverage, purity, spn) %>% 
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
  ggh4x::facet_nested(coverage+purity ~  spn + sample , scales = 'free') +
  theme_light() + theme(legend.position = 'bottom', axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))
plt

ggsave(filename = 'validation_CNAqc.png', plot = plt, width = 18, height = 7, units = 'in', dpi = 400)
