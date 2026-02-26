library(tidyverse)
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/SCOUT/colors.R')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils.R")


# spn = 'SPN03'
# cov = 50
# pur = .3

# df = tibble()
# for (spn in paste0('SPN0', 1:7)){
#   for (cov in c(50, 100, 150)){
#     for (pur in c(0.3, 0.6, 0.9)){
# 
#       r_mobster <- readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/validation_subclonal/tables/table_mobster_',spn,'_',cov, 'x_', pur, 'p_mutect2_ascat.rds')) %>%
#         select(sample_id,SNV_rate_mobster) %>% distinct()
# 
#       r_process <- readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/validation_subclonal/tables/table_process_univariate_w_private_',spn,'_',cov, 'x_', pur, 'p_mutect2_ascat.rds')) %>%
#         select(sample_id, SNV_rate_process) %>% distinct() %>% filter(!is.na(SNV_rate_process))
# 
#       join <- left_join(r_mobster, r_process) %>% mutate(coverage = cov, purity = pur, spn = spn)
#       df = bind_rows(df, join)
#     }
#   }
# }
# #saveRDS(df, '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/muts_rate.rds')

df <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/muts_rate.rds')
df_f <- df %>% 
  filter(!is.na(SNV_rate_mobster), !is.na(SNV_rate_process))  %>% 
  filter(!(spn %in% c('SPN02', 'SPN05') & SNV_rate_process == 1e-8),
         !(sample_id %in% c('SPN07_SPN07_2.1','SPN07_SPN07_2.2') & SNV_rate_process == 1e-8))

# SPN01
clone3 <- c(growth = 0.01, death = 0.01)
clone4 <- c(growth_rates = 2.5 , death_rates = 0.01)
df_rates <- tibble(spn = 'SPN01', sample = 'SPN01_SPN01_1.1', )

plt_mut_rate <- ggplot() +
  geom_rect(aes(xmin = -9.5, xmax = -8.5, ymin = -9.5, ymax = -8.5), fill = 'steelblue', alpha = .1) + 
  geom_rect(aes(xmin = -8.5, xmax = -7.5, ymin = -8.5, ymax = -7.5), fill = 'goldenrod', alpha = .1) + 
  geom_rect(aes(xmin = -7.5, xmax = -6, ymin = -7.5, ymax = -6), fill = 'coral3', alpha = .1) + 
  geom_abline(color = 'gray', linetype = 2, linewidth = .3) +
  geom_jitter(data = df_f, aes(y= log10(SNV_rate_mobster/2), x = log10(SNV_rate_process), col = spn), width = .1) +
  scale_color_manual('SPN', values = SPN_colors) + 
  my_ggplot_theme() +
  xlab('log10(mu ProCESS)') +
  ylab('log10(mu MOBSTER)') +
  xlim(-9.5, -6) + 
  ylim(-9.5, -6) +
  theme(panel.grid.minor = element_blank())
  
 