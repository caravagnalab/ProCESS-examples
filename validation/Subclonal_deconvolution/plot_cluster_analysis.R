library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

metrics_drivers= readRDS(file.path(save_path, "metrics_tables/metrics_drivers_clonal_vs_subclonal.rds"))

## Import themes ####
source('/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/figures/figure3/utils_plot.R')

SPN_colors = c("01"='steelblue', "02"='seagreen', "03"='goldenrod', 
               "04"='coral', "05"="magenta4","06"='palevioletred', "07"='indianred3')


t1 = metrics_drivers %>%
  select(spn, purity, coverage, cna_caller, vcf_caller,
         tool,
         TP_c, FP_c, FN_c,
         TP_c_blind, FP_c_blind, FN_c_blind,
         precision_c, recall_c, 
         precision_c_blind, recall_c_blind) %>%
  mutate(spn = as.character(spn),
         spn = str_extract(spn, "\\d+"),
         spn = str_sub(spn, -2, -1),
         spn=factor(spn)) %>%
  pivot_longer(
    cols = matches("_(c_blind|c)$"),
    names_to = c("metric", "type"),
    names_pattern = "^(.*)_(c_blind|c)$",
    values_to = "value"
  ) %>%
  mutate(type = ifelse(type=='c', 'Interpreted', 'Blind'))

### Precision blind vs interpreted all tools ####
p = t1 %>%
  filter(metric %in% c('precision')) %>%
  mutate(metric = factor(metric,
                         levels = c("precision")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("viber", "pyclonevi", 'mobster')),
         type=factor(type,
                     levels = c('Blind','Interpreted'))) %>%
  ggplot(aes(x = type, y = value)) +
  geom_boxplot(outliers = F) +
  stat_summary(
    aes(x = type, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3)
  ) +
  stat_summary(
    aes(x = type, color = spn, group=spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3),
    geom = 'line'
  ) +
  # ggh4x::facet_nested(~"Analysis type" + ~factor(dim, levels=c('Univariate', 'Multivariate'))) +
  xlab("Cluster Filtering")+
  ylab("Mutation Clustering Precision")+
  scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()

p
ggsave(filename = paste0(save_path, "plots/metrics/cluster_analysis/precision_blind_vs_interpreted_all_tools.pdf"), 
       plot = p, device="pdf", width=4, height=3, units="in")

saveRDS(p, file = paste0(save_path, "plots/metrics/cluster_analysis/rds/precision_blind_vs_interpreted_all_tools.rds"))

### Precision blind vs interpreted univariate vs multivariate####

p = t1 %>%
  filter(metric %in% c('precision')) %>%
  mutate(metric = factor(metric,
                         levels = c("precision")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("viber", "pyclonevi", 'mobster')),
         type=factor(type,
                     levels = c('Blind','Interpreted')),
         dim = ifelse(tool == 'mobster', 'Univariate', 'Multivariate')) %>%
  ggplot(aes(x = type, y = value)) +
  geom_boxplot(outliers = F) +
  stat_summary(
    aes(x = type, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3)
  ) +
  stat_summary(
    aes(x = type, color = spn, group=spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3),
    geom = 'line'
  ) +
  ggh4x::facet_nested(~"Analysis type" + ~factor(dim, levels=c('Univariate', 'Multivariate'))) +
  xlab("Cluster Filtering")+
  ylab("Mutation Clustering Precision")+
  scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()

p
ggsave(filename = paste0(save_path, "plots/metrics/cluster_analysis/precision_blind_vs_interpreted_univ_vs_multi.pdf"), 
       plot = p, device="pdf", width=4, height=3, units="in")

saveRDS(p, file = paste0(save_path, "plots/metrics/cluster_analysis/rds/precision_blind_vs_interpreted_univ_vs_multi.rds"))

### Precision blind vs interpreted with tool ####
p = t1 %>%
  filter(metric %in% c('precision')) %>%
  mutate(metric = factor(metric,
                         levels = c("precision")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("viber", "pyclonevi", 'mobster')),
         type=factor(type,
                     levels = c('Blind','Interpreted'))) %>%
  ggplot(aes(x = type, y = value)) +
  geom_boxplot(outliers = F) +
  stat_summary(
    aes(x = type, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3)
  ) +
  stat_summary(
    aes(x = type, color = spn, group=spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3),
    geom = 'line'
  ) +
  ggh4x::facet_nested( ~factor(tool, levels=c('mobster','viber', 'pyclonevi'))) +
  xlab("Cluster Filtering")+
  ylab("Mutation Clustering Precision")+
  scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()
p

ggsave(filename = paste0(save_path, "plots/metrics/cluster_analysis/precision_blind_vs_interpreted_w_tool.pdf"), 
       plot = p, device="pdf", width=4, height=3, units="in")
saveRDS(p, file = paste0(save_path, "plots/metrics/cluster_analysis/rds/precision_blind_vs_interpreted_w_tool.rds"))


### Precision and recall blind vs interpreted ####

p = t1 %>%
  filter(metric %in% c('precision', 'recall')) %>%
  mutate(metric = factor(metric,
                         levels = c("precision", 'recall')),
         purity = factor(purity),
         coverage = factor(coverage),
         type=factor(type,
                     levels = c('Blind','Interpreted'))) %>%
  ggplot(aes(x = type, y = value)) +
  geom_boxplot(outliers = F) +
  stat_summary(
    aes(x = type, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3)
  ) +
  stat_summary(
    aes(x = type, color = spn, group=spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3),
    geom = 'line'
  ) +
  ggh4x::facet_nested(~"Metric" + ~factor(metric, levels=c('precision', 'recall'))) +
  xlab("Cluster Filtering")+
  ylab("Metric value")+
  scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()

p

ggsave(filename = paste0(save_path, "plots/metrics/cluster_analysis/precision_recall_blind_vs_interpreted_all_tools.pdf"), 
       plot = p, device="pdf", width=4, height=3, units="in")

saveRDS(p, file = paste0(save_path, "plots/metrics/cluster_analysis/rds/precision_recall_blind_vs_interpreted_all_tools.rds"))

### Precision blind vs interpreted with tool ####
p = t1 %>%
  filter(metric %in% c('precision', 'recall')) %>%
  mutate(metric = factor(metric,
                         levels = c("precision", "recall")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c('mobster',"viber", "pyclonevi")),
         type=factor(type,
                     levels = c('Blind','Interpreted'))) %>%
  ggplot(aes(x = type, y = value)) +
  geom_boxplot(outliers = F) +
  stat_summary(
    aes(x = type, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3)
  ) +
  stat_summary(
    aes(x = type, color = spn, group=spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3),
    geom = 'line'
  ) +
  ggh4x::facet_nested(
    ~"Metric" + ~factor(metric, levels=c('precision', 'recall'))
    ~"Tool" + ~factor(tool, levels=c('mobster','viber', 'pyclonevi'))) +
  xlab("Cluster Filtering")+
  # ylab("Mutation Clustering Precision")+
  scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()
p

ggsave(filename = paste0(save_path, "plots/metrics/cluster_analysis/precision_recall_blind_vs_interpreted_w_tool.pdf"), 
       plot = p, device="pdf", width=4, height=3, units="in")
saveRDS(p, file = paste0(save_path, "plots/metrics/cluster_analysis/rds/precision_recall_blind_vs_interpreted_w_tool.rds"))
