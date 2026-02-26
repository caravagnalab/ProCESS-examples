library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

metrics_table = readRDS(file.path(save_path, "metrics_tables/metrics_drivers_clonal_vs_subclonal.rds"))

## Import themes ####
source('/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/figures/figure3/utils_plot.R')

SPN_colors = c("01"='steelblue', "02"='seagreen', "03"='goldenrod', 
               "04"='coral', "05"="magenta4","06"='palevioletred', "07"='indianred3')

palette_blind_int = RColorBrewer::brewer.pal(n=4, name="Dark2")

metrics_table_nmi_long = metrics_table %>%
  select(spn, purity, coverage, cna_caller, vcf_caller,
         tool,
         sample,
         nmi_blind,
         nmi_interpreted) %>%
  mutate(spn = as.character(spn),
         spn = str_extract(spn, "\\d+"),
         spn = str_sub(spn, -2, -1),
         spn=factor(spn)) %>% 
  pivot_longer(
    cols = c(nmi_blind, nmi_interpreted),
    names_to = "metric",
    values_to = "value"
  ) %>% 
  mutate(metric_label=case_when(
    metric=='nmi_blind'~"NMI blind",
    metric=='nmi_interpreted'~"NMI interpreted"
  )) 

mean_NMI = metrics_table_nmi_long %>%
  mutate(
    metric_label = factor(metric_label, levels = c("NMI blind", "NMI interpreted")),
    tool = factor(tool, levels = c("viber", "pyclonevi", "mobster"))
  ) %>%
  group_by(tool, spn, metric_label) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")


nmi_plot = metrics_table_nmi_long %>%
  mutate(tool = case_when(
    tool == 'viber' ~ 'VIBER',
    tool == 'mobster' ~ 'MOBSTER',
    tool == 'pyclonevi' ~ 'PyClone-VI',
  )) %>% 
  mutate(metric_label = factor(metric_label,
                               levels = c("NMI blind", "NMI interpreted")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("MOBSTER", "VIBER", 'PyClone-VI'))) %>%
  ggplot(aes(x = spn, y = value)) +
  geom_boxplot(aes(fill = metric_label), position = position_dodge(width = 0.8), outliers = T, outlier.size = .3, alpha=.9) +
  # geom_point(
  #   aes(group = metric_label),
  #   position = position_dodge(width = 0.8),
  #   size = 1.5
  # ) +
  ggh4x::facet_nested(~factor(tool, levels=c("MOBSTER", "VIBER", 'PyClone-VI'))) +
  xlab("SPN")+ 
  ylab("NMI")+
  scale_fill_manual('Metric', values = c('coral2', 'cadetblue')) +  
  scale_color_manual('Metric', values = c('coral2', 'cadetblue')) + 
  my_ggplot_theme()

nmi_plot

