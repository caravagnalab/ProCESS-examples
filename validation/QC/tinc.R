rm(list = ls())
options(bitmapType='cairo')
library(tidyverse)
library(optparse)
library(caret)
library(dplyr)
library(patchwork)

source('../../getters/process_getters.R')
source('../../getters/tumourevo_getters.R')

option_list <- list( 
  make_option(c("-v", "--variantcaller"), type="character", default='mutect2', help="variant caller"),
  make_option(c("-c", "--cnacaller"), type="character", default='ascat', help="cna caller")
)
param <- parse_args(OptionParser(option_list=option_list))

SPN <- paste0('SPN0', seq(1,7))
COV <- c(50, 100, 150, 200)
PUR <- c(0.3, 0.6, 0.9)

final_table <- tibble()
for (spn in SPN){
  print(spn)
  for (cov in COV){
    for (pur in PUR){
      samples <- get_sample_names(spn)
      
      table <- lapply(samples, FUN = function(sample){
        file <- get_tumourevo_qc(spn = spn, coverage = cov, purity = pur, tool = 'tinc', vcf_caller = param$variantcaller, cna_caller = param$cnacaller, sample = sample)

        if (length(file) > 0){
          if (file.exists(file$fit_rds)){
            data <- readRDS(file$fit_rds)
            tmp <- tibble(TIT = data$TIT, TIN = data$TIN, sample = sample)
            return(tmp)
          }
        }
      }) %>% bind_rows()
      
      if (nrow(table) > 1){
        table$coverage = as.numeric(cov)
        table$purity = as.numeric(pur)
        table$vc = param$variantcaller
        table$cnc = 'ascat'
        table$SPN = spn
        table$N = length(samples)
        table <- table %>% mutate(error = abs(purity-TIT))
      }
      final_table <- bind_rows(final_table, table)
    }
  }
}
#saveRDS(object = table, file = paste0(out, param$cov,'x_',param$pur,'p_',param$variantcaller,'_ascat.rds'))

#sp <- ggpubr::ggscatter(final_table, x = "TIT", y = "purity",
#                color = "SPN", palette = "jco", position = 'jitter',
#                add = "reg.line") 
#corr <- sp + ggpubr::stat_cor(aes(color = SPN), label.x = 0.6, label.y.npc = c(0.27, 0.3))


v1 <- final_table %>% 
  ggplot() +
  geom_abline(linewidth = 0.5, col = 'gray') +
  geom_point(aes(y = TIT, x = purity, col = SPN), size = 3) +
  geom_point(aes(x = purity, y = purity), size = 3, shape = 8) +
  facet_grid(. ~ coverage) + 
  ylab('TIT (TINC)') +
  xlab('purity (ProCESS)') + 
  scale_color_manual(values = c('steelblue', 'seagreen', 'goldenrod', 'coral', 'palevioletred', 'indianred3')) +
  xlim(0,1) +
  ylim(0,1) +
  theme_bw()
v1

v2 <- ggplot(final_table, aes(x = SPN, y = error, col = SPN)) +
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