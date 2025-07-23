rm(list = ls())
options(bitmapType='cairo')
library(tidyverse)
library(vcfR)
library(optparse)
library(caret)
library(dplyr)
library(patchwork)

source('vcf_parser.R')
source('utils.R')
source('../../getters/sarek_getters.R')
source('../../getters/process_getters.R')

option_list <- list( 
  make_option(c("-s", "--SPN"), type="character", default='SPN03', help="SPN name")
)


param <- parse_args(OptionParser(option_list=option_list))
dir <- '/orfeo/scratch/cdslab/shared/SCOUT/'
spn <- param$SPN
vcf <- paste0(dir, spn, "/validation/germline/vcf")
output <- paste0(dir, spn, "/validation/germline/rds")
report <- paste0(dir, spn, "/validation/germline/report")
dir.create(output, recursive = T, showWarnings = F)
dir.create(report, recursive = T, showWarnings = F)

tool <- c('haplotypecaller', 'freebayes', 'strelka')
list_tool <- list()
list_tool_pass <- list()

for (t in tool){
  print(t)
  vcf_list = lapply(1:22, FUN = function(chr){
    print(chr)
    rds_file = paste0(output, '/chr', chr, '_normal_sample.', t, '.rds')
    
    if (!file.exists(rds_file)){
      if (t == 'haplotypecaller'){
        vcf_file = get_sarek_vcf_file(spn = spn, type = 'normal', caller = 'haplotypecaller', coverage = 50, purity = .3)$vcf
        rds_vcf = parse_HaplotypeCaller(file = vcf_file, out_file = rds_file, save = T)[[paste0(spn, '_normal_sample')]]$mutations
      } else if (t == 'freebayes') {
        vcf_file = get_sarek_vcf_file(spn = spn, type = 'normal', caller = 'freebayes', coverage = 50, purity = .3)$vcf
        rds_vcf = parse_freebayes(file = vcf_file, out_file = rds_file, save = T, cutoff = 0.3)[[paste0(spn, '_normal_sample')]]$mutations
      } else if (t == 'strelka'){
        vcf_file = get_sarek_vcf_file(spn = spn, type = 'normal', caller = 'strelka', coverage = 50, purity = .3)$variants_vcf
        rds_vcf = parse_strelka(file = vcf_file, out_file = rds_file, save = T)[[paste0(spn, '_normal_sample')]]$mutations
      }
    } else  {
      rds_vcf = readRDS(rds_file)[[paste0(spn, '_normal_sample')]]$mutations    
    }
    
    rds_vcf = rds_vcf %>% 
      dplyr::select(chr, from, to, ref, alt, NV, DP, BAF, QUAL, FILTER) %>% 
      dplyr::mutate(mutationID = paste(chr,from, sep = ':'))
    return(rds_vcf)
  })
  
  list_tool[[t]] =  bind_rows(vcf_list)
  list_tool_pass[[t]] = bind_rows(vcf_list) %>% dplyr::filter(FILTER == "PASS" & !is.na(FILTER))
}

process_normal <- readRDS(get_mutations(spn = spn, type = 'normal')) %>% 
  filter(classes =='germinal') %>% 
  dplyr::mutate(chr = paste0('chr', chr),
                mutationID = paste(chr,chr_pos, sep = ':')) %>% 
  dplyr::rename(BAF = normal_sample.VAF,
         DP = normal_sample.coverage,
         NV = normal_sample.occurrences)

N = 1e5
DP_list = list(process_normal$DP %>% sample(1e5), list_tool_pass$freebayes$DP %>% sample(1e5), list_tool_pass$haplotypecaller$DP %>% sample(1e5), list_tool_pass$strelka$DP%>% sample(1e5))
BAF_list = list(process_normal$BAF%>% sample(1e5), list_tool_pass$freebayes$BAF%>% sample(1e5), list_tool_pass$haplotypecaller$BAF%>% sample(1e5), list_tool_pass$strelka$BAF%>% sample(1e5))

labels = c("ProCESS", "haplotypecaller", "strelka", "freebayes")
colors = c("ProCESS"="deepskyblue3", "haplotypecaller"="coral3", "strelka"="palegreen4", "freebayes"="maroon")

DP_density = plot_density_comparison_multi(DP_list, labels, colors, x_label = "Read Depth", n=N)
DP_ecdf = plot_ecdf_comparison_multi(DP_list, labels, colors, x_label = "Read Depth", n=N)

BAF_density = plot_density_comparison_multi(BAF_list, labels, colors, x_label = "BAF", n=N)
BAF_ecdf = plot_ecdf_comparison_multi(BAF_list, labels, colors, x_label = "BAF", n=N)

#freebayes
merged_freebayes = merge_datasets(snp_caller = list_tool$freebayes, ground_truth = process_normal)
y_true_freebayes = as.numeric(factor((merged_freebayes$BAF.races > 0), levels=c(FALSE, TRUE))) - 1
y_pred_freebayes = as.numeric(factor((merged_freebayes$BAF.caller > 0), levels=c(FALSE, TRUE))) - 1
metrics_freebayes = compute_metrics(y_true_freebayes, y_pred_freebayes) %>% mutate(tool = "freebayes")

merged_freebayes_pass = merge_datasets(snp_caller = list_tool_pass$freebayes, ground_truth = process_normal)
baf_comparison_freebayes <- merged_freebayes_pass %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(true_positive) %>%  # Only analyze true positive variants
  dplyr::group_by(chr.caller) %>%   # Group by chromosome
  dplyr::summarise(
    n_variants = dplyr::n(),                             # Count of variants per chromosome
    cor_coeff = safe_cor_test(BAF.caller, BAF.races),    # Correlation coefficient
    RMSE = safe_rmse(BAF.races, BAF.caller),             # Root mean square error
    .groups = 'drop'
  )


merged_strelka = merge_datasets(snp_caller = list_tool$strelka, ground_truth = process_normal)
y_true_strelka = as.numeric(factor((merged_strelka$BAF.races > 0), levels=c(FALSE, TRUE))) - 1
y_pred_strelka = as.numeric(factor((merged_strelka$BAF.caller > 0), levels=c(FALSE, TRUE))) - 1
metrics_strelka = compute_metrics(y_true_strelka, y_pred_strelka) %>% mutate(tool = "strelka")

merged_strelka_pass = merge_datasets(snp_caller = list_tool_pass$strelka, ground_truth = process_normal)
baf_comparison_strelka <- merged_strelka_pass %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(true_positive) %>%  # Only analyze true positive variants
  dplyr::group_by(chr.caller) %>%   # Group by chromosome
  dplyr::summarise(
    n_variants = dplyr::n(),                             # Count of variants per chromosome
    cor_coeff = safe_cor_test(BAF.caller, BAF.races),    # Correlation coefficient
    RMSE = safe_rmse(BAF.races, BAF.caller),             # Root mean square error
    .groups = 'drop'
  )


merged_haplotypecaller = merge_datasets(snp_caller = list_tool$haplotypecaller, ground_truth = process_normal)
y_true_haplotypecaller = as.numeric(factor((merged_haplotypecaller$BAF.races > 0), levels=c(FALSE, TRUE))) - 1
y_pred_haplotypecaller = as.numeric(factor((merged_haplotypecaller$BAF.caller > 0), levels=c(FALSE, TRUE))) - 1
metrics_haplotypecaller = compute_metrics(y_true_haplotypecaller, y_pred_haplotypecaller) %>% mutate(tool = "haplotypecaller")

merged_haplotypecaller_pass = merge_datasets(snp_caller = list_tool_pass$haplotypecaller, ground_truth = process_normal)
baf_comparison_haplotypecaller <- merged_haplotypecaller_pass %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(true_positive) %>%  # Only analyze true positive variants
  dplyr::group_by(chr.caller) %>%   # Group by chromosome
  dplyr::summarise(
    n_variants = dplyr::n(),                             # Count of variants per chromosome
    cor_coeff = safe_cor_test(BAF.caller, BAF.races),    # Correlation coefficient
    RMSE = safe_rmse(BAF.races, BAF.caller),             # Root mean square error
    .groups = 'drop'
  )

all_metric <- bind_rows(metrics_freebayes, metrics_strelka, metrics_haplotypecaller) %>% mutate(spn = spn)
all_baf <- bind_rows(baf_comparison_haplotypecaller %>% mutate(tool = 'haplotypecaller'), 
                     baf_comparison_strelka %>% mutate(tool = 'strelka'), 
                     baf_comparison_freebayes %>% mutate(tool = 'freebayes')) %>% mutate(spn = spn)

saveRDS(object = list(report_metrics=all_metric, baf_metric = all_baf), file = paste0(report,'/normal_metrics.rds'))

print('start plotting')
x = list(
  "ProCESS" = process_normal %>% dplyr::pull(mutationID),
  "haplotypecaller" = list_tool_pass$haplotypecaller  %>% dplyr::pull(mutationID),
  "strelka" = list_tool_pass$strelka  %>% dplyr::pull(mutationID),
  "freebayes" = list_tool_pass$freebayes  %>% dplyr::pull(mutationID)
)
upset_plot = ggVennDiagram::ggVennDiagram(x, force_upset = TRUE)
upset_plot = ggplotify::as.ggplot(upset_plot)



metric_plot <- all_metric %>% 
  pivot_longer(cols = c(Accuracy, Sensitivity, Precision, Recall, F1_Score)) %>% 
  ggplot() +
  geom_point(aes(x = name, y = value, col = tool), size = 3) +
  geom_line(aes(x = name, y = value, col = tool, group = tool), linetype = 2) +
  scale_color_manual(values = colors) +
  xlab('metric') + 
  theme_bw()


baf_corr_plot <- all_baf %>% 
  ggplot() +
  geom_boxplot(aes(x = tool, y = cor_coeff, fill = tool)) +
  scale_fill_manual(values = colors) +
  theme_bw()

baf_rmse_plot <- all_baf %>% 
  ggplot() +
  geom_boxplot(aes(x = tool, y = RMSE, fill = tool)) +
  scale_fill_manual(values = colors) +
  theme_bw()

title = spn
design <- 'AABB\nCCDD\nEEEE\nEEEE\nFFGH'
report_plot <- DP_density + theme(legend.position = 'none') + DP_ecdf + theme(legend.position = 'none') +
  BAF_density + theme(legend.position = 'none') + BAF_ecdf +
  upset_plot + metric_plot + theme(legend.position = 'none') +
  baf_corr_plot + theme(legend.position = 'none') +
  baf_rmse_plot + theme(legend.position = 'none') +
  plot_layout(design = design, guides = 'collect') +
  plot_annotation(title)

ggsave(plot = report_plot, filename = paste0(report,'/normal.png'), dpi = 500, width = 11, height = 12, units = 'in')
