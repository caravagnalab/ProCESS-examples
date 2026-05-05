rm(list = ls())
options(bitmapType='cairo')
library(tidyverse)
library(optparse)
library(caret)
library(dplyr)
library(patchwork)

source('../../getters/process_getters.R')
source('../../getters/tumourevo_getters.R')
df_all_combs_SPN_cna<- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/cna/cna_all.rds")
df_all_combs_SPN_cna <- df_all_combs_SPN_cna %>% 
  dplyr::rename(cna_caller=tool) %>% 
  mutate(cna_caller=case_when(cna_caller=="Sequenza"~"sequenza",
                              cna_caller=="ASCAT" ~ "ascat",
                              TRUE~cna_caller)) %>% 
  
  mutate(purity=as.numeric(true_purity),
         cn_purity=as.numeric(purity),
         coverage=as.numeric(coverage)) %>% 
  select(sample,spn,coverage,cn_purity,class,cna_caller,purity)

fga_df <-readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/fga_df.rds") %>% 
  dplyr::rename(sample=spn) %>% 
  mutate(fga_class=case_when(sample%in%c("SPN01_1.1","SPN01_1.3","SPN06_3.2")~"WGD",
                             TRUE ~ fga_class))

SPN <- paste0('SPN0', c(1,3,4,5,6,7))
coverages <- c(50, 100, 150)
purities <- c(0.3, 0.6, 0.9)
cna_callers <- c("sequenza","ascat")
variant_callers <- c("mutect2","strelka")
params_grid = expand.grid(coverages, purities,cna_callers,variant_callers)
colnames(params_grid) = c("coverage", "purity","cna_caller","snv_caller")


validation_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION/"
validate_spns <- list()
for (spn in SPN){
  print(spn)
  validate_tables <- list()
  for (i in 1:nrow(params_grid)){
    print(i)
    cov <- params_grid$coverage[i]
    pur <- params_grid$purity[i]
    cna_caller <- params_grid$cna_caller[i]
    snv_caller <- params_grid$snv_caller[i]
    
    samples <- get_sample_names(spn)
    
    table <- lapply(samples, FUN = function(sample){
      file <- get_tumourevo_qc(spn = spn, coverage = cov, purity = pur, tool = 'tinc', 
                               vcf_caller =snv_caller, cna_caller =cna_caller, sample = sample)
      
      if (length(file) > 0){
        if (file.exists(file$fit_rds)){
          data <- readRDS(file$fit_rds)
          tmp <- tibble(TIT = data$TIT, TIN = data$TIN, sample = sample)
          return(tmp)
        }
      }
    }) %>% bind_rows()
    table$coverage = as.numeric(cov)
    table$purity = as.numeric(pur)
    table$vcf_caller = snv_caller
    table$cna_caller = cna_caller
    table$spn = spn
    table$N = length(samples)
    table <- table %>% mutate(error = abs(purity-TIT))
    validate_tables[[i]] <- table
  }
  validate_spns[[spn]] <- do.call("bind_rows",validate_tables)
}
final_table_tinc <- do.call("bind_rows",validate_spns)
final_table_tinc <- final_table_tinc %>% 
  rename(spn=SPN) %>% 
  left_join(df_all_combs_SPN_cna) %>% 
  left_join(fga_df)

v1 <- final_table_tinc %>% 
  filter(fga_class!="WGD") %>% 
  ggplot() +
  geom_abline(linewidth = 0.5, col = 'gray') +
  geom_point(aes(y = TIT, x = purity, col=spn, shape=fga_class), size = 3) +
  geom_point(aes(x = purity, y = purity), size = 3, shape = 8) +
  facet_grid(. ~ class) + 
  ylab('TIT (TINC)') +
  xlab('purity (ProCESS)') + 
  scale_color_manual(values = c('steelblue', 'seagreen', 'goldenrod', 'coral', 'palevioletred', 'indianred3')) +
  xlim(0,1) +
  ylim(0,1) +
  theme_bw()
v1
v2 <- final_table_tinc %>% 
  filter(fga_class!="WGD") %>% 
  ggplot( aes(x = fga_class, y = error, col = class)) +
  geom_boxplot() +
#  geom_jitter() +
  ylab('|ProCESS - TINC|') +
  # scale_color_manual(values = c('steelblue', 'seagreen', 'goldenrod', 'coral', 'palevioletred', 'indianred3'))+
  facet_grid(.~cna_caller) +
  my_ggplot_theme()

v2

v2 <- ggplot(final_table_tinc, aes(x = spn, y = error, col = spn)) +
  geom_boxplot() +
  geom_jitter() +
  ylab('|ProCESS - TINC|') +
  scale_color_manual(values = c('steelblue', 'seagreen', 'goldenrod', 'coral', 'palevioletred', 'indianred3')) +
  facet_grid(coverage ~ purity) +
  theme_bw()

# per spn, cov, purity
metrics_1 <- final_table %>%
  mutate(
    abs_err = abs(TIT - purity),
    sq_err  = (TIT - purity)^2,
    rel_err = if_else(purity > 0, abs_err / purity, NA_real_)  # avoid div-by-zero
  ) %>%
  group_by(SPN, coverage, purity) %>%
  summarise(
    n        = n(),
    MAE      = mean(abs_err, na.rm = TRUE),
    RMSE     = sqrt(mean(sq_err, na.rm = TRUE)),
    MedianAE = median(abs_err, na.rm = TRUE),
    MAPE     = mean(rel_err, na.rm = TRUE),   # % error (scale-free)
    Bias     = mean(TIT - purity, na.rm = TRUE),  # signed error
    .groups  = "drop"
  ) %>%
  arrange(SPN, coverage, purity)
  
# per spn, cov
metrics_2 <- final_table %>%
  group_by(SPN, coverage) %>%
  summarise(
    n           = n(),
    MAE         = mean(abs(TIT - purity), na.rm = TRUE),
    RMSE        = sqrt(mean((TIT - purity)^2, na.rm = TRUE)),
    Bias        = mean(TIT - purity, na.rm = TRUE),
    Pearson_r   = cor(TIT, purity, use = "complete.obs", method = "pearson"),
    Spearman_rho= cor(TIT, purity, use = "complete.obs", method = "spearman"),
    .groups     = "drop"
  ) %>%
  arrange(SPN, coverage)


rmse1 <- metrics_1 %>% 
  ggplot() +
  geom_point(aes(x = as.factor(coverage), y = RMSE, col = SPN), size = 2) +
  scale_color_manual(values = c('steelblue', 'seagreen', 'goldenrod', 'coral', 'palevioletred', 'indianred3')) +
  facet_grid(.~purity)+
  xlab('Coverage') + 
  theme_bw()

rmse2 <- ggplot(metrics_1, aes(x = purity, y = RMSE,
                               color = SPN, linetype = factor(coverage))) +
  geom_line(aes(group = interaction(SPN, coverage)), size = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = unique(metrics_1$purity)) +
  scale_color_manual(values = c('steelblue', 'seagreen', 'goldenrod', 'coral', 'palevioletred', 'indianred3')) +
  labs(x = "Tumor purity (truth)",
       y = "RMSE (TIT vs purity)",
       color = "SPN",
       linetype = "Coverage") +
  theme_minimal(base_size = 12)

ggsave(filename = 'validation_TINC_v1.png', plot = v1, width = 7, height = 3, units = 'in', dpi = 600)
ggsave(filename = 'validation_TINC_v2.png', plot = v2, width = 11, height = 4, units = 'in', dpi = 600)
ggsave(filename = 'RMSE_TINC.png', plot = rmse2, width = 6, height = 6, units = 'in', dpi = 600)
saveRDS(object = final_table, file = 'table_TINC.rds')
saveRDS(object = metrics_1, file = 'metric_spn_cov_pur_TINC.rds')
saveRDS(object = metrics_2, file = 'metric_spn_cov_TINC.rds')