library(tidyverse)
library(patchwork)

base = "/orfeo/scratch/cdslab/shared/SCOUT/VALIDATION/"
SPNS <- paste0('SPN0', c(1,2,3,4,5,6,7))
CALLERS <- c('haplotypecaller', 'strelka', 'freebayes')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
df_all_SPN_germline = lapply(SPNS, function(i) {
  file =  paste0('/orfeo/scratch/cdslab/shared/SCOUT/VALIDATION/', i ,'/germline/report/normal_metrics.rds')
  if (file.exists(file)) {
    metrics = readRDS(file)[['report_metrics']]
  }
}) %>% do.call("bind_rows", .)
#saveRDS(object = df, file = 'germline.rds')

#df <- readRDS('germline.rds')


metrics <- c("Accuracy","Sensitivity","Precision","F1_Score")
df_all_SPN_germline_orig <- df_all_SPN_germline

df_all_SPN_germline  <- df_all_SPN_germline %>%
  pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") %>%
  group_by(tool, metric) %>%
  summarise(
    n     = n(),
    mean  = mean(value, na.rm = TRUE),
    sd    = sd(value,  na.rm = TRUE),
    se    = sd / sqrt(n),
    tcrit = qt(0.975, df = n - 1),           # 95% CI
    ci_lo = mean - tcrit * se,
    ci_hi = mean + tcrit * se,
    .groups = "drop"
  )

caller_line_germline <- ggplot(df_all_SPN_germline, aes(x = metric, y = mean,
                              col = as.factor(tool), fill = as.factor(tool), group = tool)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = as.factor(tool)), alpha = 0.5, color = NA, position = position_dodge(width = 0.3)) +
  geom_line(size = .7, position = position_dodge(width = 0.3), linetype = 3) +
  geom_point(size = 2, position = position_dodge(width = 0.3)) +
  scale_color_manual('Caller', values = col_germline_tools) +
  scale_fill_manual('Caller', values = col_germline_tools) +
  ylab('Value') +
  xlab('Metric') +
  ylim(0.9, 1) + 
  my_ggplot_theme()


caller_col_germline <- ggplot(df_all_SPN_germline, aes(x = metric, y = mean,
                                       col = as.factor(tool), fill = as.factor(tool))) +
  geom_col(position = position_dodge())+
  scale_color_manual('Caller', values = col_germline_tools) +
  scale_fill_manual('Caller', values = col_germline_tools) +
  ylab('Value') +
  xlab('Metric') +
  my_ggplot_theme()
  # ylim(0.9, 1)


df_mut_germline = lapply(SPNS, function(i) {
  file =  paste0('/orfeo/scratch/cdslab/shared/SCOUT/VALIDATION/', i ,'/germline/report/normal_metrics.rds')
  if (file.exists(file)) {
    metrics = readRDS(file)[['baf_metric']]
  }
}) %>% do.call("bind_rows", .)


caller_box_germline <- df_mut_germline %>% 
  pivot_longer(c(cor_coeff, RMSE)) %>% 
  mutate(name = ifelse(name == 'RMSE', 'RMSE', 'Correlation')) %>% 
  ggplot(aes(x = name, y = value, fill = as.factor(tool),colour = as.factor(tool))) +
  geom_boxplot(outliers = F,alpha = 0.3) + 
  scale_color_manual('Germline caller', values = col_germline_tools) +
  scale_fill_manual('Germline caller', values = col_germline_tools) +
  ylab('Value') +
  xlab('Metric') +
  my_ggplot_theme()+
  facet_wrap(.~name, scales = 'free') +
  my_ggplot_theme()+
  theme(axis.text.x = element_text(size = 0))

#plt_germline <- caller_line_germline + caller_box_germline + plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom')
#ggsave(filename = '/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/germline/germline.pdf', 
#       width = 8, height = 3.5, units = 'in', dpi = 200)
