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


# Plots driver analysis with only TP ####
t = metrics_drivers %>%
  select(spn, purity, coverage, cna_caller, vcf_caller,
         tool,
         clonal_TP_d, clonal_FP_d, clonal_FN_d,
         subclonal_TP_d, subclonal_FP_d, subclonal_FN_d) %>%
  mutate(spn = as.character(spn),
         spn = str_extract(spn, "\\d+"),
         spn = str_sub(spn, -2, -1),
         spn=factor(spn)) %>%
  pivot_longer(
    cols = matches("^(clonal|subclonal)_.*_d$"),
    names_to = c("type", "metric"),
    names_pattern = "^(clonal|subclonal)_(.*)_d$",
    values_to = "value"
  )%>% 
  mutate(type = ifelse(type=='clonal', 'Clonal', 'Subclonal'))

pr = t %>%
  filter(metric %in% c("TP", "FP", "FN")) %>%
  pivot_wider(
    names_from = metric,
    values_from = value
  ) %>%
  mutate(
    precision = TP / (TP + FP),
    recall    = TP / (TP + FN),
    f1 = 2*((precision*recall)/(precision+recall))
  ) %>% 
  group_by(spn, purity, coverage, type) %>% 
  mutate(mean_precision = mean(precision, na.rm = TRUE),
         mean_recall = mean(recall, na.rm = TRUE),
         mean_f1 = mean(f1, na.rm = TRUE)) %>%
  ungroup()%>%
  distinct(
    spn, purity, tool, type, cna_caller, vcf_caller, # remove differences per coverage because I only want 1 point in the scatter
    mean_precision, mean_recall,mean_f1,
    .keep_all = TRUE
  )




#### Only Driver facet Clonal-Subclonal on TP clusters ####
p_driver = pr %>%
  mutate(purity = factor(purity),
         coverage = factor(coverage)) %>% 
  ggplot(aes(x = type, y = f1)) +
  geom_boxplot(outliers = F, col = 'gray80', fill = 'gainsboro', alpha = .4) +
  stat_summary(
    aes(x = type, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3),
    size = .2, show.legend = F
  ) +
  xlab("Cluster Type")+
  ylab("Driver assignment F1 Score")+
  scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()


# ggsave(filename = paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/driver.pdf"), 
#        plot = p_driver, device="pdf", width=4, height=3, units="in")

#saveRDS(p, file = paste0(save_path, "plots/metrics/driver_analysis/rds/precision_recall_clonal_vs_subclonal_only_TPc.rds"))
