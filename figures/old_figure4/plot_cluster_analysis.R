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
         tool,sample,
         TP_c, FP_c, FN_c,
         TP_c_blind, FP_c_blind, FN_c_blind,
         precision_c, recall_c, 
         precision_c_blind, recall_c_blind) %>%
  mutate(f1_c_blind = 2*((precision_c_blind*recall_c_blind)/(precision_c_blind+recall_c_blind))) %>% 
  mutate(f1_c = 2*((precision_c*recall_c)/(precision_c+recall_c))) %>% 
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

table_cluster <- t1 %>% 
  filter(metric == 'f1') %>% 
  group_by(type, tool) %>% 
  summarise(mean = mean(value, na.rm = T),
            median = median(value, na.rm = T),
            sd = sd(value, na.rm = T))
# write.table(x = table_cluster, quote = F, sep = '\t', row.names = F,
#             file = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/cluster_table.tsv')


### Precision blind vs interpreted with tool ####
p_cluster = t1 %>%
  mutate(type = ifelse(type == 'Interpreted', 'Driver\nInformed', 'Driver\nBlind')) %>% 
  mutate(tool = case_when(
    tool == 'viber' ~ 'VIBER',
    tool == 'mobster' ~ 'MOBSTER',
    tool == 'pyclonevi' ~ 'PyClone-VI',
  )) %>% 
  filter(metric %in% c('f1')) %>%
  mutate(metric = factor(metric,
                         levels = c("f1")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("MOBSTER", "VIBER", 'PyClone-VI')),
         type=factor(type,
                     levels = c('Driver\nBlind','Driver\nInformed'))) %>%
  ggplot(aes(x = type, y = value)) +
  geom_boxplot(outliers = F, col = 'gray40', fill = 'gainsboro', alpha = .4) +
  stat_summary(
    aes(x = type, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3),
    size = .2, show.legend = T
  ) +
  stat_summary(
    aes(x = type, color = spn, group=spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3),
    geom = 'line',
    linewidth=.3, show.legend = F
  ) +
  ggh4x::facet_nested(~ factor(tool, levels=c("MOBSTER", "VIBER", 'PyClone-VI'))) +
  xlab("Cluster")+
  ylab("Mutation Clustering F1 Score")+
  scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()
p_cluster

# ggsave(filename = paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/cluster.pdf"), 
#        plot = p_cluster, device="pdf", width=4, height=3, units="in")
#saveRDS(p, file = paste0(save_path, "plots/metrics/cluster_analysis/rds/precision_blind_vs_interpreted_w_tool.rds"))

# EDF: With purity and coverage ####

p_cluster_purity = t1 %>%
  mutate(type = ifelse(type == 'Interpreted', 'Driver\nInformed', 'Driver\nBlind')) %>% 
  mutate(tool = case_when(
    tool == 'viber' ~ 'VIBER',
    tool == 'mobster' ~ 'MOBSTER',
    tool == 'pyclonevi' ~ 'PyClone-VI',
  )) %>% 
  filter(metric %in% c('f1')) %>%
  mutate(metric = factor(metric,
                         levels = c("f1")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("MOBSTER", "VIBER", 'PyClone-VI')),
         type=factor(type,
                     levels = c('Driver\nBlind','Driver\nInformed'))) %>%
  ggplot(aes(x = type, y = value)) +
  geom_boxplot(outliers = F, col = 'gray40', fill = 'gainsboro', alpha = .4) +
  stat_summary(
    aes(x = type, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3),
    size = .2, show.legend = T
  ) +
  stat_summary(
    aes(x = type, color = spn, group=spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3),
    geom = 'line',
    linewidth=.3, show.legend = F
  ) +
  ggh4x::facet_nested(~ factor(tool, levels=c("MOBSTER", "VIBER", 'PyClone-VI'))
                      ~ factor(purity, levels=c(0.3, 0.6, 0.9))) +
  xlab("Cluster")+
  ylab("Mutation Clustering F1 Score")+
  scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()

p_cluster_purity


p_cluster_purity = t1 %>%
  mutate(type = ifelse(type == 'Interpreted', 'Driver Informed', 'Driver Blind')) %>% 
  mutate(tool = case_when(
    tool == 'viber' ~ 'VIBER',
    tool == 'mobster' ~ 'MOBSTER',
    tool == 'pyclonevi' ~ 'PyClone-VI',
  )) %>% 
  filter(metric %in% c('f1')) %>%
  mutate(metric = factor(metric,
                         levels = c("f1")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("MOBSTER", "VIBER", 'PyClone-VI')),
         type=factor(type,
                     levels = c('Driver Blind','Driver Informed'))) %>%
  ggplot(aes(x = purity, y = value, fill = type)) +
  geom_boxplot(outliers = F,alpha = .9) +
  ggh4x::facet_nested(~ factor(tool, levels=c("MOBSTER", "VIBER", 'PyClone-VI'))) +
  scale_fill_manual('Cluster', values = c('coral2', 'cadetblue')) +  
  scale_color_manual('Cluster', values = c('coral2', 'cadetblue')) + 
  xlab("Purity")+
  ylab("Mutation Clustering F1 Score")+
  # scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()

p_cluster_purity


# EDF: number of clusters ####

t2 = metrics_drivers %>%
  select(spn, purity, coverage, cna_caller, vcf_caller,
         tool,sample,
         n_cluster_process,
         n_cluster_blind,
         n_cluster_interpreted,
         n_cluster_interpreted_driver) %>%
  mutate(spn = as.character(spn),
         spn = str_extract(spn, "\\d+"),
         spn = str_sub(spn, -2, -1),
         spn=factor(spn)) %>% 
  mutate(relative_error_blind = abs((n_cluster_blind- n_cluster_interpreted_driver)/n_cluster_interpreted_driver),
         relative_error_interpreted = abs((n_cluster_interpreted- n_cluster_interpreted_driver)/n_cluster_interpreted_driver)) %>% 
  filter(relative_error_blind != Inf, relative_error_interpreted != Inf)


t3 = t2 %>% 
  pivot_longer(
    cols = c(relative_error_blind, relative_error_interpreted),
    names_to = "type",
    values_to = "relative_error"
  ) %>% 
  mutate(type = ifelse(type == 'relative_error_blind', 'Driver blind', 'Driver interpreted')) %>% 
  mutate(tool = case_when(
    tool == 'viber' ~ 'VIBER',
    tool == 'mobster' ~ 'MOBSTER',
    tool == 'pyclonevi' ~ 'PyClone-VI',
  )) %>% 
  group_by(purity, type, tool) %>% 
  summarise(mean = mean(relative_error),
            std = sd(relative_error)) %>% 
  ungroup()

p = ggplot(t3, aes(x=as.factor(purity), y=mean, fill = type)) + 
  geom_bar(stat="identity", color="gray40", 
           position=position_dodge(), alpha = .9) +
  # geom_errorbar(aes(ymin=pmax(mean - std, 0), ymax=mean+std), width=.2,
  #               position=position_dodge(.9)) +
  ggh4x::facet_nested(~ factor(tool, levels=c("MOBSTER", "VIBER", 'PyClone-VI'))) +
  scale_fill_manual('Cluster', values = c('coral2', 'cadetblue')) + 
  ylab('Relative Error N of Cluster')+
  xlab('Purity') +
  guides(alpha = "none")+
  my_ggplot_theme()

p
p_cluster_purity
