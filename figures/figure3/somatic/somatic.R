library(ggplot2)
library(tidyverse)
library(patchwork)

base = "/orfeo/scratch/cdslab/shared/SCOUT/VALIDATION/"
SPNS <- paste0('SPN0', c(1,2,3,4,6,7))
CALLERS <- c('mutect2', 'strelka', 'freebayes')
COVERAGES <- c(50,100,150)
PURITY <- c(0.3, 0.6, 0.9)
MUTS <- c('SNV', 'INDEL')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
params_grid = expand.grid(COVERAGES, PURITY, MUTS, SPNS)
names(params_grid) <- c('cov', 'pur', 'muts', 'spn')

df_all_SPN_somatic = lapply(1:nrow(params_grid), function(i) {
  mut_type = params_grid[i,]$muts
  coverage = params_grid[i,]$cov
  purity = params_grid[i,]$pur
  spn = params_grid[i,]$spn

  file =  paste0('/orfeo/scratch/cdslab/shared/SCOUT/VALIDATION/', spn ,'/somatic/', coverage,'x_',purity,'p/report/',mut_type,'/metrics_new_binning.rds')

  if (file.exists(file)) {
    metrics = readRDS(file)

    parsed_metrics = lapply(names(metrics), function(caller_name) {
      caller_metrics = metrics[[caller_name]]
      sample_names = names(caller_metrics)
      lapply(sample_names, function(sample_id) {
        FP = caller_metrics[[sample_id]]$detection_summary["False Positive"]
        d = caller_metrics[[sample_id]]$performance_table %>%
          dplyr::mutate(sample_id = sample_id, SPN = spn) %>%
          dplyr::mutate(purity = as.numeric(purity), coverage = as.numeric(coverage))
        d$false_positive[1] = FP
        d
      }) %>% do.call("bind_rows", .) %>%
        dplyr::mutate(caller = caller_name)
    }) %>% do.call("bind_rows", .)

    dplyr::bind_cols(parsed_metrics, params_grid[i,])
  }
}) %>% do.call("bind_rows", .)

df_all_SPN_somatic <- df_all_SPN_somatic %>% 
  mutate(CCF_bin=case_when(CCF_bin=="Clonal"~"95-100%",
                           TRUE ~ CCF_bin))
#saveRDS(object = df, file = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/somatic/somatic.rds')
# 
# df <- readRDS(
#   '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/somatic/somatic.rds')
# 



# coverage_line <- df_all_SPN_somatic %>% 
#   dplyr::group_by(sample_id, caller, spn, pur, muts) %>% 
#   dplyr::select(sensitivity, cov, CCF_bin) %>% 
#   dplyr::mutate(normalized_sens = sensitivity / max(sensitivity)) %>% 
#   ggplot(mapping = aes(x = as.numeric(as.factor(CCF_bin)), y = normalized_sens, col = as.factor(cov), fill = as.factor(cov))) +
#   scale_color_manual('Coverage', values = coverage_colors) + 
#   scale_fill_manual('Coverage', values = coverage_colors) + 
#   ylab('Normalized sensitivity') + 
#   xlab('CCF_bin') + 
#   geom_smooth(method = "loess") +
#   scale_x_continuous(labels = unique(df$CCF_bin), breaks=seq(2,length(unique(df$CCF_bin))+1,1)) +   facet_wrap(~muts) +
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+xlab("CCF bin")
coverage_line_boxes <- df_all_SPN_somatic %>% 
  # dplyr::group_by(sample_id, pur, spn, cov, muts) %>% 
  # dplyr::select(sensitivity, cov, CCF_bin) %>% 
  # dplyr::mutate(normalized_sens = sensitivity / max(sensitivity)) %>% 
  dplyr::mutate(normalized_sens = sensitivity) %>% 
  ggplot(mapping = aes(x = as.factor(CCF_bin), y = normalized_sens, col = as.factor(cov), fill = as.factor(cov))) +
  stat_summary(fun = "median", fun.min=function(x) boxplot.stats(x)$stats[2],fun.max=function(x) boxplot.stats(x)$stats[4],aes(colour = as.factor(cov)), 
               linewidth = 1, size=0.5,position = position_dodge(width = 0.5))+ 
  # geom_smooth(aes(group = cov), method = "loess", se = FALSE,linewidth = 0.5,alpha = 1) +  # trend line
  scale_color_manual('Coverage', values = coverage_colors) + 
  scale_fill_manual('Coverage', values = coverage_colors) + 
  ylab('') +
  xlab('CCF_bin') +
  geom_hline(yintercept = 0.5,
             linetype = "dashed", color = "grey40")+
  facet_wrap(~muts) +
  my_ggplot_theme()+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  theme(axis.text.x = element_text(angle = 30, vjust = 1),
        axis.title.x = element_blank())


# purity_line <- df_all_SPN_somatic %>% 
#   dplyr::group_by(sample_id, caller, spn, cov, muts) %>% 
#   dplyr::select(sensitivity, pur, CCF_bin) %>% 
#   dplyr::mutate(normalized_sens = sensitivity / max(sensitivity)) %>% 
#   ggplot(mapping = aes(x = as.numeric(as.factor(CCF_bin)), y = normalized_sens, col = as.factor(pur), fill = as.factor(pur))) +
#   scale_color_manual('Purity', values = purity_colors) + 
#   scale_fill_manual('Purity', values = purity_colors) + 
#   ylab('Normalized sensitivity') + 
#   xlab('CCF_bin') + 
#   scale_x_continuous(labels = unique(df$CCF_bin), breaks=seq(2,length(unique(df$CCF_bin))+1,1)) +   geom_smooth(method = "loess") +
#   facet_wrap(~muts) +
#   theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),axis.title.x = element_blank())
purity_line_boxes <- df_all_SPN_somatic %>% 
  # dplyr::group_by(sample_id, pur, spn, cov, muts) %>% 
  # dplyr::select(sensitivity, pur, CCF_bin) %>% 
  # dplyr::mutate(normalized_sens = sensitivity / max(sensitivity)) %>% 
  dplyr::mutate(normalized_sens = sensitivity) %>% 
  ggplot(mapping = aes(x = as.factor(CCF_bin), y = normalized_sens, col = as.factor(pur), fill = as.factor(pur))) +
  stat_summary(fun = "median", fun.min=function(x) boxplot.stats(x)$stats[2],fun.max=function(x) boxplot.stats(x)$stats[4],aes(colour =  as.factor(pur)), 
               linewidth = 1, size=0.5,position = position_dodge(width = 0.5))+
  # geom_smooth(aes(group = pur), method = "loess", se = FALSE,linewidth = 0.5,alpha = 1) +  # trend line
  scale_color_manual('Purity', values = purity_colors) + 
  scale_fill_manual('Purity', values = purity_colors) + 
  ylab('Sensitivity') + 
  geom_hline(yintercept = 0.5,
             linetype = "dashed", color = "grey40")+
  # xlab('CCF_bin') + 
  facet_wrap(~muts) +
  my_ggplot_theme()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x =element_blank(), 
        axis.title.x = element_blank())+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
# caller_line <- df_all_SPN_somatic %>% 
#   dplyr::group_by(sample_id, pur, spn, cov, muts) %>% 
#   dplyr::select(sensitivity, caller, CCF_bin) %>% 
#   dplyr::mutate(normalized_sens = sensitivity / max(sensitivity)) %>% 
#   ggplot(mapping = aes(x = as.numeric(as.factor(CCF_bin)), y = normalized_sens, col = as.factor(caller), fill = as.factor(caller))) +
#   scale_color_manual('Caller', values = col_somatic_tools) + 
#   scale_fill_manual('Caller', values = col_somatic_tools) + 
#   ylab('Normalized sensitivity') + 
#   xlab('CCF_bin') + 
#   geom_smooth(method = "loess") +
#   scale_x_continuous(labels = unique(df$CCF_bin), breaks=seq(2,length(unique(df$CCF_bin))+1,1)) + 
#   facet_wrap(~muts) +
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 30, vjust = 1),axis.title.x = element_blank())
caller_line_boxes <- df_all_SPN_somatic %>% 
  # dplyr::group_by(sample_id, pur, spn, cov, muts) %>% 
  # dplyr::select(sensitivity, caller, CCF_bin) %>% 
  dplyr::mutate(normalized_sens = sensitivity) %>% 
  ggplot(mapping = aes(x = as.factor(CCF_bin), y = normalized_sens, col = as.factor(caller), fill = as.factor(caller))) +
  stat_summary(fun = "median", fun.min=function(x) boxplot.stats(x)$stats[2],fun.max=function(x) boxplot.stats(x)$stats[4],aes(colour = caller), 
               linewidth = 1, size=0.5,position = position_dodge(width = 0.5))+
  # geom_smooth(aes(group = caller), method = "loess", se = FALSE,linewidth = 0.5,alpha = 1) +  # trend line
  scale_color_manual('Somatic caller', values = col_somatic_tools) + 
  scale_fill_manual('Somatic caller', values = col_somatic_tools) + 
  ylab('') +
  geom_hline(yintercept = 0.5,
             linetype = "dashed", color = "grey40")+
  facet_wrap(~muts) +
  my_ggplot_theme()+
  theme(axis.text.x =element_blank(),axis.ticks.x =element_blank(), 
        axis.title.x = element_blank())


# somatic_panel <- caller_line + purity_line + coverage_line + plot_layout(nrow = 3) + plot_annotation(tag_levels = 'A')

somatic_panel <- caller_line_boxes + purity_line_boxes + coverage_line_boxes + plot_layout(nrow = 3,heights = 1) &
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) 

ggsave(plot = somatic_panel, filename = '/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/somatic/somatic.pdf',width = 8, height = 10, units = 'in')
#ggsave(plot = line, filename = '/orfeo/cephfs/scratch/area/lvaleriani/tesi/sarek/all_sarek.pdf',width = 6, height = 8, units = 'in')


tmp <- df %>% 
  dplyr::mutate(CCF_bin = case_when(
    CCF_bin == "Clonal" ~ "Clonal",
    CCF_bin %in% c("10-25%", "25-50%", "50-99%") ~ "10-99%",
    TRUE ~ "0-10%"
  )) %>% 
  select(sample_id, CCF_bin, pur, spn, cov, muts, sensitivity, caller) %>% 
  group_by(sample_id, CCF_bin, pur, spn, cov, muts, caller) %>%
  summarise(sensitivity = mean(sensitivity))
  
caller <- tmp %>% 
  dplyr::group_by(sample_id, CCF_bin, pur, spn, cov, muts) %>% 
  dplyr::select(sensitivity, caller) %>% 
  dplyr::mutate(normalized_sens = sensitivity / max(sensitivity)) %>% 
  ggplot(mapping = aes(x = CCF_bin, y = normalized_sens, fill = caller)) +
  geom_boxplot(outliers = F) +
  ylab('Normalized sensitivity') + 
  scale_fill_manual('Caller', values = method_colors) + 
  facet_wrap(~muts) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1.2, hjust=1),axis.title.x = element_blank())


purity <- tmp %>% 
  dplyr::group_by(sample_id, CCF_bin, caller, spn, cov, muts) %>% 
  dplyr::select(sensitivity, pur) %>% 
  dplyr::mutate(normalized_sens = sensitivity / max(sensitivity)) %>% 
  ggplot(mapping = aes(x = CCF_bin, y = normalized_sens, fill = as.factor(pur))) +
  geom_boxplot(outliers = F) +
  scale_fill_manual('Purity', values = purity_colors) + 
  ylab('Normalized sensitivity') + 
  facet_wrap(~muts) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 30, vjust = 1.2, hjust=1),axis.title.x = element_blank())

coverage <- tmp %>% 
  dplyr::group_by(sample_id, CCF_bin, caller, spn, pur, muts) %>% 
  dplyr::select(sensitivity, cov) %>% 
  dplyr::mutate(normalized_sens = sensitivity / max(sensitivity)) %>% 
  ggplot(mapping = aes(x = CCF_bin, y = normalized_sens, fill = as.factor(cov))) +
  geom_boxplot(outliers = F) +
  scale_fill_manual('Coverage', values = coverage_colors) + 
  ylab('Normalized sensitivity') + 
  facet_wrap(~muts) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 30, vjust = 1.2, hjust=1),axis.title.x="CCF bin")

boxplot <- caller + purity + coverage + plot_layout(nrow = 3)
# 
# tmp_fp <- df %>% 
#   dplyr::mutate(CCF_bin = case_when(
#     CCF_bin == "Clonal" ~ "Clonal",
#     CCF_bin %in% c("10-25%", "25-50%", "50-99%") ~ "10-99%",
#     TRUE ~ "0-10%"
#   )) %>% 
#   select(sample_id, CCF_bin, pur, spn, cov, muts, false_positive, true_positives, caller) %>% 
#   group_by(sample_id, CCF_bin, pur, spn, cov, muts, caller) %>%
#   dplyr::summarise(FP = sum(false_positive, na.rm = TRUE), TP = sum(true_positives, na.rm = TRUE)) %>% 
#   dplyr::mutate(FDR = FP / (TP + FP))
# 
# fp <- tmp_fp %>% 
#   dplyr::group_by(caller, muts, cov, pur, spn) %>% 
#   dplyr::select(FDR, caller) %>% 
#   dplyr::mutate(normalized_fdr = FDR / max(FDR)) %>% 
#   ggplot(mapping = aes(x = CCF_bin, y = normalized_fdr, fill = caller)) +
#   geom_boxplot() + 
#   facet_wrap(~muts) +
#   theme_minimal() +
#   scale_fill_manual(values = method_colors) +
#   labs(x = "Caller", y = "False Discovery Rate", fill = "Caller")



