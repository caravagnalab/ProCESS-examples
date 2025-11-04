library(tidyverse)
library(patchwork)

base = "/orfeo/scratch/cdslab/shared/SCOUT/VALIDATION/"
SPNS <- paste0('SPN0', c(1,2,3,4,6,7))
CALLERS <- c('haplotypecaller', 'strelka', 'freebayes')

# df = lapply(SPNS, function(i) {
#   file =  paste0('/orfeo/scratch/cdslab/shared/SCOUT/VALIDATION/', i ,'/germline/report/normal_metrics.rds')
#   if (file.exists(file)) {
#     metrics = readRDS(file)[['report_metrics']]
#   }
# }) %>% do.call("bind_rows", .)
# saveRDS(object = df, file = 'germline.rds')

df <- readRDS('germline.rds')
method_colors = c(
  "ProCESS" = "gray80",
  "mutect2" = "lightsteelblue",
  "mutect2 (all)" = "steelblue4",  # darker steelblue
  "strelka" = "coral",
  "strelka (all)" = "coral4",     # darker coral
  "freebayes" = "#8FBC8B",
  "freebayes (all)" = "#228B22",   # darker green (ForestGreen )
  'haplotypecaller' = 'palevioletred'
)

metrics <- c("Accuracy","Sensitivity","Precision","F1_Score")

df  <- df %>%
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

caller_line <- ggplot(df, aes(x = metric, y = mean,
                              col = as.factor(tool), fill = as.factor(tool), group = tool)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = as.factor(tool)), alpha = 0.5, color = NA) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual('Caller', values = method_colors) +
  scale_fill_manual('Caller', values = method_colors) +
  ylab('Value') +
  xlab('Metric') +
  theme_minimal() +
  ylim(0.9, 1) #+ 
#theme(axis.text.x = element_text(angle = 30, vjust = 1))

caller_line


df_mut = lapply(SPNS, function(i) {
  file =  paste0('/orfeo/scratch/cdslab/shared/SCOUT/VALIDATION/', i ,'/germline/report/normal_metrics.rds')
  if (file.exists(file)) {
    metrics = readRDS(file)[['baf_metric']]
  }
}) %>% do.call("bind_rows", .)


caller_box <- df_mut %>% 
  pivot_longer(c(cor_coeff, RMSE)) %>% 
  mutate(name = ifelse(name == 'RMSE', 'RMSE', 'Correlation')) %>% 
  ggplot(aes(x = name, y = value, col = as.factor(tool))) +
  geom_boxplot(outliers = F) + 
  scale_color_manual('Caller', values = method_colors) +
  scale_fill_manual('Caller', values = method_colors) +
  ylab('Value') +
  xlab('Metric') +
  theme_minimal() +
  facet_wrap(.~name, scales = 'free') +
  theme(axis.text.x = element_text(size = 0)) 

plt <- caller_line + caller_box + plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom')
ggsave(filename = 'germline_all.pdf', 
       width = 8, height = 3.5, units = 'in', dpi = 200)
