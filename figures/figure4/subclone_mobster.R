setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')
source('/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/utils_plot.R')


coverage_list = c(50, 100, 150)
purity_list = c(0.9, 0.6, 0.3)
vcf_caller_list = c("strelka")
cna_caller_list = c("ascat")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04','SPN05', 'SPN06', 'SPN07')
combs = expand.grid(spn = spn_list,
                    coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller = cna_caller_list)

# i=1
# results <- lapply(1:nrow(combs), FUN = function(i){
#   print(i)
#   sp = combs[i, "spn"]
#   cov = combs[i, "coverage"]
#   pur = combs[i, "purity"]
#   vcf_call = combs[i, 'vcf_caller']
#   cna_call = combs[i, 'cna_caller']
#   
#   
#   samples <- get_sample_names(sp)
#   inf_data = lapply(samples, FUN = function(s){
#     
#     tmp_sample = paste(sp, s, sep = '_')
#     tmp_class = table %>% filter(sample_id == tmp_sample)
#     if (sp == 'SPN02'){
#       file = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/',sp,'/tumourevo/',cov, 'x_', pur, 'p_', vcf_call,'_',cna_call,
#                     '/subclonal_deconvolution/mobster/SCOUT/',sp,'/',tmp_sample,'/',paste('SCOUT', sp, tmp_sample, 'mobster_best_fit.rds', sep = '_'))
#     } else {cna_call
#       file = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/',sp,'/tumourevo/',cov, 'x_', pur, 'p_', vcf_call,'_',cna_call,
#                     '/subclonal_deconvolution/mobster/SCOUT/',sp,'/',tmp_sample,'/',paste('SCOUT', sp, tmp_sample, 'mobster_best_fit.rds', sep = '_'))
#     }
#     if (file.exists(file)){
#       fit = readRDS(file)
#       data = tibble(sample_id = tmp_sample, 
#                     n_subclones = sum(!unique(fit$Clusters$cluster) %in% c('C1', 'Tail')))
#       return(data)
#     }
# 
#   }) %>% bind_rows()
#   
#   if (nrow(inf_data) > 0){
#     inf_data = inf_data %>% 
#       mutate(with_subclone = ifelse(n_subclones>0, T, F)) %>% 
#       mutate(coverage = cov, purity = pur,
#              vcf_caller = vcf_call, cna_caller = cna_call)
#     
#   }
# }) %>% bind_rows()
# 
# saveRDS(results, file = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/performance_subclones_strelka_ascat.rds')
results_mutect_ascat = readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/performance_subclones_mutect_ascat.rds')
results_mutect_sequenza = readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/performance_subclones_mutect_sequenza.rds')
results_strelka_ascat = readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/performance_subclones_strelka_ascat.rds')
results = bind_rows(results_mutect_ascat, results_mutect_sequenza, results_strelka_ascat)

table = readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/class_table.rds')

all_samples <- c()
for (sp in spn_list){
  samples <- get_sample_names(sp)
  all_samples <- c(all_samples, paste0(sp, '_', samples))
}

results = results %>% 
  left_join(table) %>% 
  mutate(gt = ifelse(class == 'Monoclonal', F, T)) %>% 
  dplyr::rename(pred = with_subclone)

plt_subclone <- results %>% 
  filter(!is.na(gt)) %>% 
  mutate(TP = gt & pred,
           FP = !gt & pred,
           FN = gt & !pred,
           TN = !gt & !pred) %>% 
  group_by(subclass, purity, cna_caller, vcf_caller) %>% 
  summarise(Polyclonal = sum(TP) / (sum(TP) + sum(FN)),  #recall
          Monoclonal = sum(TN) / (sum(TN) + sum(FP))) %>%  #Specificity
  pivot_longer(cols = c(Polyclonal, Monoclonal),
                 names_to = "metric",
                 values_to = "value") %>%
  ggplot(aes(x = as.factor(purity), 
             y = value, 
             color = subclass, 
             group = subclass)) +
  geom_line() +
  geom_point(size = 2) +
  scale_color_manual('Subclass', values = c('Polyclonal' = 'goldenrod', 
                                            'Monoclonal'='#645394', 
                                            'Weak evidence' = 'khaki3',
                                            'Sampling Bias' = 'plum')) +
  ylab("Sample composition detection\nTPR") +
  xlab("Purity") +
  my_ggplot_theme() +
  ylim(0,1) + 
  ggh4x::facet_grid2(vcf_caller + cna_caller~metric)
