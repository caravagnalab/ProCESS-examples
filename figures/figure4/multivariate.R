setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')
source('/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/utils_plot.R')


# coverage_list = c(50, 100, 150)
# purity_list = c(0.9, 0.6, 0.3)
# vcf_caller_list = c("strelka", "mutect2")
# cna_caller_list = c("sequenza", "ascat")
# spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04','SPN05', 'SPN06', 'SPN07')
# combs = expand.grid(spn = spn_list,
#                     coverage = coverage_list,
#                     purity = purity_list,
#                     vcf_caller = vcf_caller_list,
#                     cna_caller = cna_caller_list)
# 
# 
class <- c('No Bias',
           'Sampling Bias',
           'Sampling Bias',
           'Sampling Bias',
           'No Bias',
           'No Bias',
           'No Bias')
class = tibble(spn = spn_list, type = class)

# i=1
# results <- lapply(1:nrow(combs), FUN = function(i){
#   print(i)
#   sp = combs[i, "spn"]
#   cov = combs[i, "coverage"]
#   pur = combs[i, "purity"]
#   vcf_call = combs[i, 'vcf_caller']
#   cna_call = combs[i, 'cna_caller']
#   
#   samples = get_sample_names(sp)
#   
#   file = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/validation_subclonal_new/tables_interpreted/viber_', sp, '_', cov, 'x_', pur, 'p_', vcf_call, '_', cna_call, '.rds')
#   if (file.exists(file)){
#     viber_cluster = readRDS(file) %>% pull(cluster_id_tool) %>% unique() %>% length()
#   }
#   
#   file = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/validation_subclonal_new/tables_interpreted/pyclonevi_', sp, '_', cov, 'x_', pur, 'p_', vcf_call, '_', cna_call, '.rds')
#   if (file.exists(file)){
#     pyclone_cluster = readRDS(file) %>% pull(cluster_id_tool) %>% unique() %>% length()
#   }
#   
#   tibble(spn = sp,
#          coverage = cov,
#          purity = pur, 
#          vcf_caller = vcf_call,
#          cna_caller = cna_call,
#          cluster_viber = viber_cluster,
#          cluster_pyclone = pyclone_cluster,
#          samples = length(samples))
#   
# }) %>% bind_rows()
# 
# saveRDS(results, file = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/performance_multivariate.rds')
results = readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/performance_multivariate.rds')
results = results %>% 
  select(-type) %>% 
  left_join(class, by = join_by(spn)) %>% 
  mutate(viber = cluster_viber/samples,
         pyclone = cluster_pyclone/samples)

plt_multivariate <- results %>% 
  pivot_longer(cols = c(viber, pyclone)) %>% 
  ggplot() +
  geom_boxplot(aes(x = as.factor(purity), y = value, col = type, fill = type), alpha =.2) +
  facet_grid(.~name) +
  my_ggplot_theme() + 
  xlab('Purity') +
  scale_fill_manual('SPN', values = c('darkseagreen4', 'orangered3')) + 
  scale_color_manual('SPN', values = c('darkseagreen4', 'orangered3')) + 
  ylab('Relative number of inferred cluster\nClusters / Samples')



results_long <- results %>%
  pivot_longer(cols = c(viber, pyclone), names_to = "tool", values_to = "value")

# Run Wilcoxon test for each tool x purity combination
significance <- results_long %>%
  group_by(tool, purity) %>%
  summarise(
    p_value = wilcox.test(
      value[type == "No Bias"],
      value[type == "Sampling Bias"]
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),  # Benjamini-Hochberg correction
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

print(significance)
