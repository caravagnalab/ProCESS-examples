library(tidyverse)
library(randnet)
library(scales)

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

# Simple final plot ####

metrics_table = readRDS(file.path(save_path, "metrics_tables/new_table_clusters_metrics_v2.rds"))
metrics_table = metrics_table %>% filter(tool!='mobster')

## Theme ####
my_theme = theme_light(base_size=12) +
  theme(legend.position="bottom",
        legend.key.size=unit(0.3, "cm"),
        panel.background=element_rect(fill="white"),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=12),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        text=element_text(size=12))

palette_raw_int = RColorBrewer::brewer.pal(n=3, name="Dark2")
palette_coverage = RColorBrewer::brewer.pal(n=3, name="Set1")
palette_spn = RColorBrewer::brewer.pal(n=length(spn_list), name="Dark2")

## NMI ####
metrics_table_nmi_long = metrics_table %>%
  select(spn, purity, coverage, cna_caller, vcf_caller,
         tool,
         nmi_raw,
         nmi_interpreted) %>%
  mutate(spn = as.character(spn),
         spn = str_extract(spn, "\\d+"),
         spn = str_sub(spn, -2, -1),
         spn=factor(spn)) %>% 
  pivot_longer(
    cols = c(nmi_raw, nmi_interpreted),
    names_to = "metric",
    values_to = "value"
  ) %>% 
  mutate(metric_label=case_when(
    metric=='nmi_raw'~"NMI raw",
    metric=='nmi_interpreted'~"NMI interpreted"
  )) 
# 
# mean_NMI = metrics_table_nmi_long %>%
#   group_by(metric_label,tool) %>%
#   summarize(mean_value = mean(value))

mean_NMI = metrics_table_nmi_long %>%
  mutate(
    metric_label = factor(metric_label, levels = c("NMI raw", "NMI interpreted")),
    tool = factor(tool, levels = c("viber", "pyclonevi"))
  ) %>%
  group_by(tool, spn, metric_label) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

metrics_table_nmi_long %>%
  mutate(
    metric_label = factor(metric_label, levels = c("NMI raw", "NMI interpreted")),
    purity = factor(purity),
    coverage = factor(coverage),
    tool = factor(tool, levels = c("viber", "pyclonevi"))
  ) %>%
  ggplot(aes(x = metric_label, y = value)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_line(
    data = mean_NMI,
    aes(color = spn, group = spn),
    linewidth = 0.6,
    alpha = 0.8
  ) +
  geom_point(
    data = mean_NMI,
    aes(color = spn),
    size = 2.5
  ) +
  ggh4x::facet_nested(~"Tool" + ~tool) +
  xlab("") +
  ylab("NMI") +
  scale_color_manual(values = palette_spn, name = "SPN") +
  my_theme

metrics_table_nmi_long %>%
  mutate(metric_label = factor(metric_label,
                               levels = c("NMI raw", "NMI interpreted")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("viber", "pyclonevi"))) %>%
  # group_by(tool, metric_label) %>%
  # mutate(spn_gpr=paste0(spn, tool, metric_label)) %>%
  # mutate(spn_gpr = fct_reorder(spn_gpr, value, .fun=mean, .desc=T)) %>%
  
  ggplot(aes(x = metric_label, y = value)) +
  geom_boxplot(outliers = F) +
  stat_summary(
    aes(x = metric_label, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3)
  ) +
  stat_summary(
    aes(x = metric_label, color = spn, group=spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3),
    geom = 'line'
  ) +
  ggh4x::facet_nested(~"Tool" + ~factor(tool, levels=c('viber', 'pyclonevi'))) +
  xlab("")+
  ylab("NMI")+
  scale_color_manual(values=palette_spn, name='SPN')+
  my_theme


## Relative error number of clusters ####
metrics_table_long_relative_error = metrics_table %>%
  select(spn, purity, coverage,tool,
         n_raw_tool,
         n_interpreted_tool,
         n_true_driver_process) %>%
  mutate(relative_raw_error = abs((n_raw_tool - n_true_driver_process)/n_true_driver_process),
         relative_interpreted_error = abs((n_interpreted_tool - n_true_driver_process)/n_true_driver_process)) %>% 
  pivot_longer(
    cols = c(relative_raw_error,
             relative_interpreted_error),
    names_to = "error",
    values_to = "value"
  ) %>% 
  mutate(error_label=case_when(
    error=='relative_raw_error'~"Raw",
    error=='relative_interpreted_error'~"Interpreted"
  ))

### Boxplots ####

# mean_relative_err = metrics_table_long_relative_error %>%
#   group_by(error_label, spn,tool) %>%
#   summarize(value = mean(value)) %>%
#   mutate(
#     mean_label = paste0(error_label, " (", round(value, 3), ")")
#   )

mean_relative_err = metrics_table_long_relative_error %>%
  mutate(
    error_label = factor(error_label, levels = c("Raw", "Interpreted")),
    tool = factor(tool, levels = c("viber", "pyclonevi"))
  ) %>%
  group_by(tool, spn, error_label) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")


metrics_table_long_relative_error %>%
  mutate(
    error_label = factor(error_label, levels = c("Raw", "Interpreted")),
    purity = factor(purity),
    coverage = factor(coverage),
    tool = factor(tool, levels = c("viber", "pyclonevi"))
  ) %>%
  ggplot(aes(x = error_label, y = value)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_line(
    data = mean_relative_err,
    aes(color = spn, group = spn),
    linewidth = 0.6,
    alpha = 0.8
  ) +
  geom_point(
    data = mean_relative_err,
    aes(color = spn),
    size = 2.5
  ) +
  ggh4x::facet_nested(~"Tool" + ~tool) +
  xlab("") +
  ylab("Relative error") +
  scale_color_manual(values = palette_spn, name = "SPN") +
  my_theme

# New plots ####

## NMI ####
metrics_table_nmi_long = metrics_table %>%
  select(spn, purity, coverage, cna_caller, vcf_caller,
         tool,
         nmi_raw,
         nmi_interpreted) %>%
  mutate(spn = as.character(spn),
         spn = str_extract(spn, "\\d+"),
         spn = str_sub(spn, -2, -1),
         spn=factor(spn)) %>% 
  pivot_longer(
    cols = c(nmi_raw, nmi_interpreted),
    names_to = "metric",
    values_to = "value"
  ) %>% 
  mutate(metric_label=case_when(
    metric=='nmi_raw'~"NMI raw",
    metric=='nmi_interpreted'~"NMI interpreted"
  )) 

mean_NMI = metrics_table_nmi_long %>%
  group_by(metric_label,tool) %>%
  summarize(mean_value = mean(value))


nmi_plot = metrics_table_nmi_long %>%
  mutate(metric_label = factor(metric_label,
                               levels = c("NMI raw", "NMI interpreted")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                levels = c("viber", "pyclonevi"))) %>%
  ggplot(aes(x = spn, y = value)) +
  geom_boxplot(aes(color = metric_label), position = position_dodge(width = 0.8)) +
  geom_point(
    aes(shape=purity, group = metric_label, alpha = coverage),
    position = position_dodge(width = 0.8),
    size = 2.5
  ) +
  ggh4x::facet_nested(~"Tool" + ~factor(tool, levels=c('viber', 'pyclonevi'))) +
  geom_hline(
    data = mean_NMI,
    aes(yintercept = mean_value, color = metric_label),
    linewidth = 0.5,
    linetype = 'dashed'
  ) +
  xlab("SPN")+ 
  ylab("NMI")+
  ylim(0,1)+
  scale_color_manual(values=palette_raw_int, name='Metric')+
  scale_shape_manual(values=c(16,17,3), name='Purity')+
  scale_alpha_manual(values=c(0.2,0.5,1), name = 'Coverage')+
  my_theme

nmi_plot
# ggsave(file.path(save_path, "plots/metrics/nmi_plot_tools.png"), nmi_plot)

## ARI ####
metrics_table_ari_long = metrics_table %>%
  select(spn, purity, coverage, tool,
         cna_caller, vcf_caller,
         ari_raw,
         ari_interpreted) %>%
  pivot_longer(
    cols = c(ari_raw, ari_interpreted),
    names_to = "metric",
    values_to = "value"
  ) %>% 
  mutate(metric_label=case_when(
    metric=='ari_raw'~"ARI raw",
    metric=='ari_interpreted'~"ARI interpreted"
  ))

mean_ARI = metrics_table_ari_long %>%
  group_by(metric_label, tool) %>%
  summarize(mean_value = mean(value))

ari_plot=metrics_table_ari_long %>%
  mutate(metric_label = factor(metric_label,
                               levels = c("ARI raw", "ARI interpreted")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("viber", "pyclonevi"))) %>%
  ggplot(aes(x = spn, y = value)) +
  geom_boxplot(aes(color = metric_label), position = position_dodge(width = 0.8)) +
  geom_point(
    aes(shape=purity, group = metric_label, alpha = coverage),
    position = position_dodge(width = 0.8),
    size = 2.5
  ) +
  ggh4x::facet_nested(~"Tool" + ~factor(tool, levels=c('viber', 'pyclonevi'))) +
  geom_hline(
    data = mean_ARI,
    aes(yintercept = mean_value, color = metric_label),
    linewidth = 0.5,
    linetype = 'dashed'
  ) +
  xlab("SPN")+ 
  ylab("NMI")+
  scale_color_manual(values=palette_raw_int, name='Metric')+
  scale_shape_manual(values=c(16,17,3), name='Purity')+
  scale_alpha_manual(values=c(0.2,0.5,1), name = 'Coverage')+
  my_theme


ari_plot

ggsave(file.path(save_path, "plots/metrics/ari_plot.png"), ari_plot)


## Relative error number of clusters ####
metrics_table_long_relative_error = metrics_table %>%
  select(spn, purity, coverage,tool,
         n_raw_tool,
         n_interpreted_tool,
         n_true_driver_process) %>%
  mutate(relative_raw_error = abs((n_raw_tool - n_true_driver_process)/n_true_driver_process),
         relative_interpreted_error = abs((n_interpreted_tool - n_true_driver_process)/n_true_driver_process)) %>% 
  pivot_longer(
    cols = c(relative_raw_error,
             relative_interpreted_error),
    names_to = "error",
    values_to = "value"
  ) %>% 
  mutate(error_label=case_when(
    error=='relative_raw_error'~"Raw",
    error=='relative_interpreted_error'~"Interpreted"
  ))
### Scatter ####

mean_relative_err = metrics_table_long_relative_error %>%
  group_by(error_label, tool) %>%
  summarize(mean_value = mean(value)) %>%
  mutate(
    mean_label = paste0(error_label, " (", round(mean_value, 3), ")")
  )

clusters_relative_error_scatter = metrics_table_long_relative_error %>%
  mutate(
    error_label = factor(error_label, levels = c("Raw", "Interpreted")),
    purity = factor(purity),
    coverage=factor(coverage)
  ) %>%
  ggplot(aes(x = n_true_driver_process, y = value, color=spn, shape=purity)) +
  ggh4x::facet_nested("Tool"+~factor(tool, levels=c('viber', 'pyclonevi'))~"Error" + 
                        ~factor(error_label, levels=c("Raw", "Interpreted"))) +
  geom_point(aes(alpha=coverage), size = 3) +
  geom_hline(
    data = mean_relative_err,
    aes(yintercept = mean_value, linetype = "Mean Error"),
    color = 'black',
    linetype = 'dashed',
    linewidth = 0.5
  ) +
  xlab('Number of simulated clones') +
  # ylab('Relative error |(#tool - #ProCESS)/#ProCESS|') +
  ylab('Relative error') +
  scale_color_manual(values=palette_spn, name='SPN')+
  scale_shape_manual(values=c(16,17,3), name='Purity')+
  scale_alpha_manual(values=c(0.2,0.5,1), name='Coverage')+
  scale_x_continuous(breaks = scales::pretty_breaks())+
  # scale_y_continuous(breaks = scales::pretty_breaks())+
  my_theme+
  theme(legend.position = 'right')

clusters_relative_error_scatter

ggsave(file.path(save_path, "plots/metrics/clusters_relative_error_scatter.png"), clusters_relative_error_scatter)

## Wasserstein distance ####
new_table = metrics_table %>%
  select(spn, purity, coverage,tool,
         n_true_driver_process,
         wasserstein_raw,
         wasserstein_interpreted) %>% 
  pivot_longer(
    cols = c(wasserstein_raw,
             wasserstein_interpreted),
    names_to = "distance",
    values_to = "value"
  ) %>% 
  mutate(distance_label=case_when(
    distance=='wasserstein_raw'~"Raw",
    distance=='wasserstein_interpreted'~"Interpreted"
  ))
### Scatter ####

# mean_relative_err = metrics_table_long_relative_error %>%
#   group_by(error_label, tool) %>%
#   summarize(mean_value = mean(value)) %>%
#   mutate(
#     mean_label = paste0(error_label, " (", round(mean_value, 3), ")")
#   )

new_table %>%
  mutate(
    error_label = factor(distance_label, levels = c("Raw", "Interpreted")),
    purity = factor(purity),
    coverage=factor(coverage)
  ) %>%
  ggplot(aes(x = n_true_driver_process, y = value, color=spn, shape=purity)) +
  ggh4x::facet_nested("Tool"+~factor(tool, levels=c('viber', 'pyclonevi'))~"Error" + 
                        ~factor(distance_label, levels=c("Raw", "Interpreted"))) +
  geom_point(aes(alpha=coverage), size = 3) +
  xlab('Number of simulated clones') +
  ylab('Wasserstein distance') +
  ylim(0,1)+
  scale_color_manual(values=palette_spn, name='SPN')+
  scale_shape_manual(values=c(16,17,3), name='Purity')+
  scale_alpha_manual(values=c(0.2,0.5,1), name='Coverage')+
  scale_x_continuous(breaks = scales::pretty_breaks())+
  my_theme+
  theme(legend.position = 'right')



# Patchwork ####

p = patchwork::wrap_plots(
  clusters_relative_error_scatter,
  ggplot() + labs(title="CCF"),
  nmi_plot + theme(legend.position = 'right'), 
  # ari_plot + labs(title="Clustering ARI"),
  ggplot() + labs(title="Signatures"),
  design="aabb\nccdd"
) &
  patchwork::plot_annotation(tag_levels="A") &
  theme(plot.tag=element_text(size=18, face="bold"))
p

ggsave(file.path(save_path, "plots/metrics/fig5.png"),p,
       width=35, height=25, units="cm")
ggsave(file.path(save_path, "plots/metrics/fig5.pdf"),p,
       device = 'pdf', width=35, height=25, units="cm")

# NMI vs relative error ####
base = metrics_table %>%
  select(spn, purity, coverage, tool,
         n_raw_tool, n_interpreted_tool, n_true_driver_process,
         nmi_raw, nmi_interpreted) %>%
  mutate(
    relative_raw_error = abs((n_raw_tool - n_true_driver_process) / n_true_driver_process),
    relative_interpreted_error = abs((n_interpreted_tool - n_true_driver_process) / n_true_driver_process)
  )

error_long = base %>%
  pivot_longer(
    cols = c(relative_raw_error, relative_interpreted_error),
    names_to = "metric",
    values_to = "error"
  ) %>%
  mutate(
    metric = recode(metric,
                    relative_raw_error = "Blind",
                    relative_interpreted_error = "Interpreted"
    )
  ) %>%
  select(spn, purity, coverage, tool, metric, error)

nmi_long = base %>%
  pivot_longer(
    cols = c(nmi_raw, nmi_interpreted),
    names_to = "metric",
    values_to = "NMI"
  ) %>%
  mutate(
    metric = recode(metric,
                    nmi_raw = "Blind",
                    nmi_interpreted = "Interpreted"
    )
  ) %>%
  select(spn, purity, coverage, tool, metric, NMI)

metrics_table_long_relative_error_nmi = left_join(
  error_long, nmi_long,
  by = c("spn", "purity", "coverage", "tool", "metric")
)

error_vs_NMI_scatter_plot = metrics_table_long_relative_error_nmi %>%
  mutate(
    error_complete = factor(metric, levels = c("Blind", "Interpreted")),
    purity = factor(purity),
    coverage=factor(coverage)
  ) %>%
  ggplot(aes(x = error, y = NMI, color=spn, shape=purity)) +
  ggh4x::facet_nested("Tool"+~factor(tool, levels=c('viber', 'pyclonevi'))~"Error" + 
                        ~factor(metric, levels=c("Blind", "Interpreted"))) +
  geom_point(aes(alpha=coverage), size = 3) +
  xlab('Relative error') +
  ylab('NMI') +
  ylim(0,1)+
  scale_color_manual(values=palette_spn, name='SPN')+
  scale_shape_manual(values=c(16,17,3), name='Purity')+
  scale_alpha_manual(values=c(0.2,0.5,1), name='Coverage')+
  my_theme+
  theme(legend.position = 'right')

error_vs_NMI_scatter_plot


metrics_table_long_relative_error_nmi %>%
  mutate(
    error_complete = factor(metric, levels = c("Blind", "Interpreted")),
    purity = factor(purity),
    coverage=factor(coverage)
  ) %>%
  ggplot(aes(x = error, y = NMI, color=spn)) +
  ggh4x::facet_nested("Tool"+~factor(tool, levels=c('viber', 'pyclonevi'))~"Error" + 
                        ~factor(metric, levels=c("Blind", "Interpreted"))) +
  geom_point(size = 1) +
  xlab('Relative error') +
  ylab('NMI') +
  ylim(0,1)+
  scale_color_manual(values=palette_spn, name='SPN')+
  # scale_shape_manual(values=c(16,17,3), name='Purity')+
  # scale_alpha_manual(values=c(0.2,0.5,1), name='Coverage')+
  my_theme+
  theme(legend.position = 'right')


ggsave(file.path(save_path, "plots/metrics/error_vs_NMI.png"),error_vs_NMI_scatter_plot,
       width=20, height=15, units="cm")


# Relative error vs Kolmogorov distance ####
base = metrics_table %>%
  select(spn, purity, coverage, tool,
         n_raw_tool, n_interpreted_tool, n_true_driver_process,
         Kolmogorov_distance) %>%
  mutate(
    relative_raw_error = abs((n_raw_tool - n_true_driver_process) / n_true_driver_process),
    relative_interpreted_error = abs((n_interpreted_tool - n_true_driver_process) / n_true_driver_process)
  )

error_long = base %>%
  pivot_longer(
    cols = c(relative_raw_error, relative_interpreted_error),
    names_to = "metric",
    values_to = "error"
  ) %>%
  mutate(
    metric = recode(metric,
                    relative_raw_error = "Raw",
                    relative_interpreted_error = "Interpreted"
    )
  ) %>%
  select(spn, purity, coverage, tool, metric, error)

distance_long = base %>%
  pivot_longer(
    cols = c(Kolmogorov_distance),
    names_to = "name",
    values_to = "distance"
  )  %>% 
  # distance_long has no metric, so create it by duplicating each row
  tidyr::crossing(metric = c("Raw", "Interpreted")) %>% 
  # alternatively, use mutate + bind_rows; crossing is cleaner
  select(spn, purity, coverage, tool, metric, distance)

metrics_table_long_relative_error_nmi = left_join(
  error_long, distance_long,
  by = c("spn", "purity", "coverage", "tool", "metric")
)

error_vs_Kolmogorv_scatter_plot = metrics_table_long_relative_error_nmi %>%
  mutate(
    error_complete = factor(metric, levels = c("Raw", "Interpreted")),
    purity = factor(purity),
    coverage=factor(coverage)
  ) %>%
  ggplot(aes(x = error, y = distance, color=spn, shape=purity)) +
  ggh4x::facet_nested("Tool"+~factor(tool, levels=c('viber', 'pyclonevi'))~"Error" + 
                        ~factor(metric, levels=c("Raw", "Interpreted"))) +
  geom_point(aes(alpha=coverage), size = 3) +
  xlab('Relative error') +
  # ylab('Relative error |(#tool - #ProCESS)/#ProCESS|') +
  ylab('Kolmogorov distance') +
  scale_color_manual(values=palette_spn, name='SPN')+
  scale_shape_manual(values=c(16,17,3), name='Purity')+
  scale_alpha_manual(values=c(0.2,0.5,1), name='Coverage')+
  my_theme+
  theme(legend.position = 'right')

error_vs_Kolmogorv_scatter_plot

ggsave(file.path(save_path, "plots/metrics/error_vs_Kdistance.png"),error_vs_Kolmogorv_scatter_plot,
       width=20, height=15, units="cm")

# Relative error vs Wasserstein distance ####
base = metrics_table %>%
  select(spn, purity, coverage, tool,
         n_raw_tool, n_interpreted_tool, n_true_driver_process,
         wasserstein_interpreted,
         wasserstein_raw) %>%
  mutate(
    relative_raw_error = abs((n_raw_tool - n_true_driver_process) / n_true_driver_process),
    relative_interpreted_error = abs((n_interpreted_tool - n_true_driver_process) / n_true_driver_process)
  )

error_long = base %>%
  pivot_longer(
    cols = c(relative_raw_error, relative_interpreted_error),
    names_to = "metric",
    values_to = "error"
  ) %>%
  mutate(
    metric = recode(metric,
                    relative_raw_error = "Raw",
                    relative_interpreted_error = "Interpreted"
    )
  ) %>%
  select(spn, purity, coverage, tool, metric, error)

distance_long = base %>%
  pivot_longer(
    cols = c(wasserstein_raw, wasserstein_interpreted),
    names_to = "metric",
    values_to = "distance"
  ) %>%
  mutate(
    metric = recode(metric,
                    wasserstein_raw = "Raw",
                    wasserstein_interpreted = "Interpreted"
    )
  ) %>%
  select(spn, purity, coverage, tool, metric, distance)

metrics_table_long_relative_error_wasserstein = left_join(
  error_long, distance_long,
  by = c("spn", "purity", "coverage", "tool", "metric")
)

metrics_table_long_relative_error_wasserstein %>%
  mutate(
    error_complete = factor(metric, levels = c("Raw", "Interpreted")),
    purity = factor(purity),
    coverage=factor(coverage)
  ) %>%
  ggplot(aes(x = error, y = distance)) +
  ggh4x::facet_nested("Tool"+~factor(tool, levels=c('viber', 'pyclonevi'))~"Error" + 
                        ~factor(metric, levels=c("Raw", "Interpreted"))) +
  geom_point(aes(color=spn),size = 2) +
  geom_smooth(color='black', method="lm")+
  ggpubr::stat_cor(aes(label=after_stat(rr.label)))+
  xlab('Relative error') +
  ylab('Wasserstein distance') +
  ylim(0,1)+
  scale_color_manual(values=palette_spn, name='SPN')+
  my_theme+
  theme(legend.position = 'right')

error_vs_wasserstein_scatter_plot

ggsave(file.path(save_path, "plots/metrics/error_vs_Kdistance.png"),error_vs_Kolmogorv_scatter_plot,
       width=20, height=15, units="cm")

# Relative error vs Kolmogorov distance ####
base = metrics_table %>%
  select(spn, purity, coverage, tool,
         n_raw_tool, n_interpreted_tool, n_true_driver_process,
         Kolmogorov_distance) %>%
  mutate(
    relative_raw_error = abs((n_raw_tool - n_true_driver_process) / n_true_driver_process),
    relative_interpreted_error = abs((n_interpreted_tool - n_true_driver_process) / n_true_driver_process)
  )

error_long = base %>%
  pivot_longer(
    cols = c(relative_raw_error, relative_interpreted_error),
    names_to = "metric",
    values_to = "error"
  ) %>%
  mutate(
    metric = recode(metric,
                    relative_raw_error = "Raw",
                    relative_interpreted_error = "Interpreted"
    )
  ) %>%
  select(spn, purity, coverage, tool, metric, error)

distance_long = base %>%
  pivot_longer(
    cols = c(Kolmogorov_distance),
    names_to = "name",
    values_to = "distance"
  )  %>% 
  # distance_long has no metric, so create it by duplicating each row
  tidyr::crossing(metric = c("Raw", "Interpreted")) %>% 
  # alternatively, use mutate + bind_rows; crossing is cleaner
  select(spn, purity, coverage, tool, metric, distance)

metrics_table_long_relative_error_nmi = left_join(
  error_long, distance_long,
  by = c("spn", "purity", "coverage", "tool", "metric")
)

error_vs_Kolmogorv_scatter_plot = metrics_table_long_relative_error_nmi %>%
  mutate(
    error_complete = factor(metric, levels = c("Raw", "Interpreted")),
    purity = factor(purity),
    coverage=factor(coverage)
  ) %>%
  ggplot(aes(x = error, y = distance, color=spn, shape=purity)) +
  ggh4x::facet_nested("Tool"+~factor(tool, levels=c('viber', 'pyclonevi'))~"Error" + 
                        ~factor(metric, levels=c("Raw", "Interpreted"))) +
  geom_point(aes(alpha=coverage), size = 3) +
  xlab('Relative error') +
  # ylab('Relative error |(#tool - #ProCESS)/#ProCESS|') +
  ylab('Kolmogorov distance') +
  scale_color_manual(values=palette_spn, name='SPN')+
  scale_shape_manual(values=c(16,17,3), name='Purity')+
  scale_alpha_manual(values=c(0.2,0.5,1), name='Coverage')+
  my_theme+
  theme(legend.position = 'right')

error_vs_Kolmogorv_scatter_plot

ggsave(file.path(save_path, "plots/metrics/error_vs_Kdistance.png"),error_vs_Kolmogorv_scatter_plot,
       width=20, height=15, units="cm")

# NMI vs Kolmogorov distance ####
base = metrics_table %>%
  select(spn, purity, coverage, tool,
         nmi_raw, nmi_interpreted,
         Kolmogorov_distance) 

NMI_long = base %>%
  pivot_longer(
    cols = c(nmi_raw, nmi_interpreted),
    names_to = "metric",
    values_to = "NMI"
  ) %>%
  mutate(
    metric = recode(metric,
                    nmi_raw = "Raw",
                    nmi_interpreted = "Interpreted"
    )
  ) %>%
  select(spn, purity, coverage, tool, metric, NMI)

distance_long = base %>%
  pivot_longer(
    cols = c(Kolmogorov_distance),
    names_to = "name",
    values_to = "distance"
  )  %>% 
  # distance_long has no metric, so create it by duplicating each row
  tidyr::crossing(metric = c("Raw", "Interpreted")) %>% 
  # alternatively, use mutate + bind_rows; crossing is cleaner
  select(spn, purity, coverage, tool, metric, distance)

metrics_table_long_distance_nmi = left_join(
  NMI_long, distance_long,
  by = c("spn", "purity", "coverage", "tool", "metric")
)

NMI_vs_Kolmogorv_scatter_plot = metrics_table_long_distance_nmi %>%
  mutate(
    error_complete = factor(metric, levels = c("Raw", "Interpreted")),
    purity = factor(purity),
    coverage=factor(coverage)
  ) %>%
  ggplot(aes(x = NMI, y = distance, color=spn, shape=purity)) +
  ggh4x::facet_nested("Tool"+~factor(tool, levels=c('viber', 'pyclonevi'))~"Error" + 
                        ~factor(metric, levels=c("Raw", "Interpreted"))) +
  geom_point(aes(alpha=coverage), size = 3) +
  xlab('NMI') +
  # ylab('Relative error |(#tool - #ProCESS)/#ProCESS|') +
  ylab('Kolmogorov distance') +
  scale_color_manual(values=palette_spn, name='SPN')+
  scale_shape_manual(values=c(16,17,3), name='Purity')+
  scale_alpha_manual(values=c(0.2,0.5,1), name='Coverage')+
  my_theme+
  theme(legend.position = 'right')

NMI_vs_Kolmogorv_scatter_plot

ggsave(file.path(save_path, "plots/metrics/NMI_vs_Kdistance.png"),NMI_vs_Kolmogorv_scatter_plot,
       width=20, height=15, units="cm")

# Plots with no tail mutations in interpreted ####
metrics_table = metrics_table %>% filter(tool !='mobster')
## NMI ####
metrics_table_nmi_long = metrics_table %>%
  select(spn, purity, coverage, cna_caller, vcf_caller,
         tool,
         nmi_raw,
         nmi_interpreted_no_tail) %>%
  mutate(spn = as.character(spn),
         spn = str_extract(spn, "\\d+"),
         spn = str_sub(spn, -2, -1),
         spn=factor(spn)) %>% 
  pivot_longer(
    cols = c(nmi_raw, nmi_interpreted_no_tail),
    names_to = "metric",
    values_to = "value"
  ) %>% 
  mutate(metric_label=case_when(
    metric=='nmi_raw'~"NMI raw",
    metric=='nmi_interpreted_no_tail'~"NMI interpreted"
  )) 

mean_NMI = metrics_table_nmi_long %>%
  group_by(metric_label,tool) %>%
  summarize(mean_value = mean(value))


nmi_plot = metrics_table_nmi_long %>%
  mutate(metric_label = factor(metric_label,
                               levels = c("NMI raw", "NMI interpreted")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("viber", "pyclonevi"))) %>%
  ggplot(aes(x = spn, y = value)) +
  geom_boxplot(aes(color = metric_label), position = position_dodge(width = 0.8)) +
  geom_point(
    aes(shape=purity, group = metric_label, alpha = coverage),
    position = position_dodge(width = 0.8),
    size = 2.5
  ) +
  ggh4x::facet_nested(~"Tool" + ~factor(tool, levels=c('viber', 'pyclonevi'))) +
  geom_hline(
    data = mean_NMI,
    aes(yintercept = mean_value, color = metric_label),
    linewidth = 0.5,
    linetype = 'dashed'
  ) +
  xlab("SPN")+ 
  ylab("NMI")+
  scale_color_manual(values=palette_raw_int, name='Metric')+
  scale_shape_manual(values=c(16,17,3), name='Purity')+
  scale_alpha_manual(values=c(0.2,0.5,1), name = 'Coverage')+
  my_theme

nmi_plot
# ggsave(file.path(save_path, "plots/metrics/nmi_plot_no_tail.png"), nmi_plot)



## Relative error number of clusters ####
metrics_table_long_relative_error = metrics_table %>%
  select(spn, purity, coverage,tool,
         n_raw_tool,
         n_interpreted_tool_no_tail,
         n_true_driver_process) %>%
  mutate(relative_raw_error = abs((n_raw_tool - n_true_driver_process)/n_true_driver_process),
         relative_interpreted_error = abs((n_interpreted_tool_no_tail - n_true_driver_process)/n_true_driver_process)) %>% 
  pivot_longer(
    cols = c(relative_raw_error,
             relative_interpreted_error),
    names_to = "error",
    values_to = "value"
  ) %>% 
  mutate(error_label=case_when(
    error=='relative_raw_error'~"Raw",
    error=='relative_interpreted_error'~"Interpreted"
  ))

### Scatter ####
mean_relative_err = metrics_table_long_relative_error %>%
  group_by(error_label, tool) %>%
  summarize(mean_value = mean(value)) %>%
  mutate(
    mean_label = paste0(error_label, " (", round(mean_value, 3), ")")
  )

clusters_relative_error_scatter = metrics_table_long_relative_error %>%
  mutate(
    error_label = factor(error_label, levels = c("Raw", "Interpreted")),
    purity = factor(purity),
    coverage=factor(coverage)
  ) %>%
  ggplot(aes(x = n_true_driver_process, y = value, color=spn, shape=purity)) +
  ggh4x::facet_nested("Tool"+~factor(tool, levels=c('viber', 'pyclonevi'))~"Error" + 
                        ~factor(error_label, levels=c("Raw", "Interpreted"))) +
  geom_point(aes(alpha=coverage), size = 3) +
  geom_hline(
    data = mean_relative_err,
    aes(yintercept = mean_value, linetype = "Mean Error"),
    color = 'black',
    linetype = 'dashed',
    linewidth = 0.5
  ) +
  xlab('Number of simulated clones') +
  ylab('Relative error') +
  scale_color_manual(values=palette_spn, name='SPN')+
  scale_shape_manual(values=c(16,17,3), name='Purity')+
  scale_alpha_manual(values=c(0.2,0.5,1), name='Coverage')+
  scale_x_continuous(breaks = scales::pretty_breaks())+
  # scale_y_continuous(breaks = scales::pretty_breaks())+
  my_theme+
  theme(legend.position = 'right')

clusters_relative_error_scatter

metrics_table_long_relative_error %>%
  mutate(
    error_label = factor(error_label, levels = c("Raw", "Interpreted")),
    purity = factor(purity),
    coverage=factor(coverage)
  ) %>%
  ggplot(aes(x = n_true_driver_process, y = value)) +
  ggh4x::facet_nested("Tool"+~factor(tool, levels=c('viber', 'pyclonevi'))~"Error" +
                        ~factor(error_label, levels=c("Raw", "Interpreted"))) +
  # ggh4x::facet_nested(~factor(error_label, levels=c("Raw", "Interpreted")))  +
  geom_hline(
    data = mean_relative_err,
    aes(yintercept = mean_value, linetype = "Mean Error"),
    color = 'gainsboro',
    linetype = 'dashed',
    linewidth = 0.5
  ) +
  geom_point(aes(color=spn), size = 2)+
  geom_smooth(color='black', linewidth=0.5)+
  xlab('Number of simulated clones') +
  ylab('Relative error') +
  scale_color_manual(values=palette_spn, name='SPN')+
  # scale_shape_manual(values=c(16,17,3), name='Purity')+
  # scale_alpha_manual(values=c(0.2,0.5,1), name='Coverage')+
  scale_x_continuous(breaks = scales::pretty_breaks())+
  # scale_y_continuous(breaks = scales::pretty_breaks())+
  my_theme+
  theme(legend.position = 'right')


p = patchwork::wrap_plots(
  clusters_relative_error_scatter,
  ggplot() + labs(title="CCF"),
  nmi_plot + theme(legend.position = 'right'), 
  # ari_plot + labs(title="Clustering ARI"),
  ggplot() + labs(title="Signatures"),
  design="aabb\nccdd"
) &
  patchwork::plot_annotation(tag_levels="A") &
  theme(plot.tag=element_text(size=18, face="bold"))
p

ggsave(file.path(save_path, "plots/metrics/fig5_no_tail.png"),p,
       width=35, height=25, units="cm")
ggsave(file.path(save_path, "plots/metrics/fig5_no_tail.pdf"),p,
       device = 'pdf', width=35, height=25, units="cm")


# Interpreted drivers ####
palette_raw_int = RColorBrewer::brewer.pal(n=4, name="Dark2")

## NMI ####
metrics_table_nmi_long = metrics_table %>%
  select(spn, purity, coverage, cna_caller, vcf_caller,
         tool,
         nmi_raw,
         nmi_interpreted,
         nmi_interpreted_driver) %>%
  mutate(spn = as.character(spn),
         spn = str_extract(spn, "\\d+"),
         spn = str_sub(spn, -2, -1),
         spn=factor(spn)) %>% 
  pivot_longer(
    cols = c(nmi_raw, nmi_interpreted,nmi_interpreted_driver),
    names_to = "metric",
    values_to = "value"
  ) %>% 
  mutate(metric_label=case_when(
    metric=='nmi_raw'~"NMI raw",
    metric=='nmi_interpreted'~"NMI interpreted",
    metric=='nmi_interpreted_driver'~"NMI interpreted driver"
  )) 

mean_NMI = metrics_table_nmi_long %>%
  group_by(metric_label,tool) %>%
  summarize(mean_value = mean(value))


nmi_plot = metrics_table_nmi_long %>%
  mutate(metric_label = factor(metric_label,
                               levels = c("NMI raw", "NMI interpreted","NMI interpreted driver")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("viber", "pyclonevi"))) %>%
  ggplot(aes(x = spn, y = value)) +
  geom_boxplot(aes(color = metric_label), position = position_dodge(width = 0.8)) +
  geom_point(
    aes(shape=purity, group = metric_label, alpha = coverage),
    position = position_dodge(width = 0.8),
    size = 2.5
  ) +
  ggh4x::facet_nested(~"Tool" + ~factor(tool, levels=c('viber', 'pyclonevi'))) +
  geom_hline(
    data = mean_NMI,
    aes(yintercept = mean_value, color = metric_label),
    linewidth = 0.5,
    linetype = 'dashed'
  ) +
  xlab("SPN")+ 
  ylab("NMI")+
  ylim(0,1)+
  scale_color_manual(values=palette_raw_int, name='Metric')+
  scale_shape_manual(values=c(16,17,3), name='Purity')+
  scale_alpha_manual(values=c(0.2,0.5,1), name = 'Coverage')+
  my_theme

nmi_plot 

## Relative error number of clusters ####

metrics_table_long_relative_error = metrics_table %>%
  select(spn, purity, coverage,tool,
         n_raw_tool,
         n_interpreted_tool,
         n_interpreted_driver_tool,
         n_true_driver_process) %>%
  mutate(relative_raw_error = abs((n_raw_tool - n_true_driver_process)/n_true_driver_process),
         relative_interpreted_error = abs((n_interpreted_tool - n_true_driver_process)/n_true_driver_process),
         relative_interpreted_driver_error = abs((n_interpreted_driver_tool - n_true_driver_process)/n_true_driver_process)) %>% 
  pivot_longer(
    cols = c(relative_raw_error,
             relative_interpreted_error,
             relative_interpreted_driver_error),
    names_to = "error",
    values_to = "value"
  ) %>% 
  mutate(error_label=case_when(
    error=='relative_raw_error'~"Raw",
    error=='relative_interpreted_error'~"Interpreted",
    error=='relative_interpreted_driver_error'~"Interpreted driver"
  ))
### Scatter ####

mean_relative_err = metrics_table_long_relative_error %>%
  group_by(error_label, tool) %>%
  summarize(mean_value = mean(value)) %>%
  mutate(
    mean_label = paste0(error_label, " (", round(mean_value, 3), ")")
  )

clusters_relative_error_scatter = metrics_table_long_relative_error %>%
  mutate(
    error_label = factor(error_label, levels = c("Raw", "Interpreted", "Interpreted driver")),
    purity = factor(purity),
    coverage=factor(coverage)
  ) %>%
  ggplot(aes(x = n_true_driver_process, y = value, color=spn, shape=purity)) +
  ggh4x::facet_nested("Tool"+~factor(tool, levels=c('viber', 'pyclonevi'))~"Error" + 
                        ~factor(error_label, levels=c("Raw", "Interpreted", "Interpreted driver"))) +
  geom_point(aes(alpha=coverage), size = 3) +
  geom_hline(
    data = mean_relative_err,
    aes(yintercept = mean_value, linetype = "Mean Error"),
    color = 'black',
    linetype = 'dashed',
    linewidth = 0.5
  ) +
  xlab('Number of simulated clones') +
  ylab('Relative error') +
  scale_color_manual(values=palette_spn, name='SPN')+
  scale_shape_manual(values=c(16,17,3), name='Purity')+
  scale_alpha_manual(values=c(0.2,0.5,1), name='Coverage')+
  scale_x_continuous(breaks = scales::pretty_breaks())+
  # scale_y_continuous(breaks = scales::pretty_breaks())+
  my_theme+
  theme(legend.position = 'right')

clusters_relative_error_scatter
