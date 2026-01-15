library(tidyverse)
library(randnet)
library(scales)

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

source(file.path(save_path, "utils_plots_final.R"))


coverage_list = c(50,100,150)
purity_list = c(0.3,0.6,0.9)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN06', 'SPN07')
tools = c('viber', 'pyclonevi', 'mobster')

combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller=cna_caller_list,
                    spn=spn_list)

spn = 'SPN03'
coverage=100
purity=0.9
vcf_caller = "mutect2"
cna_caller = "ascat"
simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)

columns = c("spn","purity","coverage",
            "vcf_caller","cna_caller", 
            'n_true_driver_process',
            'n_raw_tool',
            'n_pipeline_interpreted_tool',
            # 'n_true_driver_tool',
            "nmi_raw",
            "nmi_interpreted",
            "ari_raw",
            "ari_interpreted",
            "tool",
            "n_pipeline_interpreted_tool_no_tail",
            "n_true_driver_process_no_tail")

#,
            # "fraction_raw_clusters", 
            # "fraction_true_driver_clusters",
            # "fraction_false_driver_clusters",
            # "fraction_all_driver_clusters")

metrics_table = data.frame(matrix(nrow = 0, ncol = length(columns)))

colnames(metrics_table) = columns

# if (!is.na(i)) {
for(tool in tools){
  for(i in 1:nrow(combs)){
    coverage = combs[i, "coverage"]
    purity = combs[i, "purity"]
    vcf_caller = combs[i, "vcf_caller"]
    cna_caller = combs[i, "cna_caller"]
    spn = combs[i, "spn"]
    
    simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
    if(tool == 'mobster' & spn=='SPN02' & ((coverage==50 & purity==0.6)|
                       (coverage==50 & purity==0.9)|
                       (coverage==100 & purity==0.6)|
                       (coverage==100 & purity==0.9)|
                       (coverage==150 & purity==0.6)|
                       (coverage==150 & purity==0.9))){
      next
    }
    
    table = readRDS(file.path(main_path, "subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds")))
    
      # add_count(cluster_id_tool_interpreted_driver, name="n_mutations_tool_interpreted_driver")
    if(tool == 'mobster'){
      
      table = table %>% 
        group_by(cluster_id_process, sample_id) %>%
        mutate(ccf_process=mean(ccf_process, na.rm=TRUE)) %>%
        ungroup()
      
      table = table %>% group_by(cluster_id_process, sample_id) %>% 
        mutate(is_clonal_process=replace(FALSE, ccf_process > 0.95, TRUE)) %>% ungroup() %>% 
        mutate(cluster_id_process_full = cluster_id_process) %>% 
        mutate(cluster_id_process=replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal'))
      
      sample_names = table$sample_id %>% unique()
      
      # Create the df to fill sample per sample
      df = data.frame(matrix(nrow = 0, ncol = length(columns)))
      colnames(df) = columns
      
      for(sample_name in sample_names){
        
        table_sample = table %>% filter(sample_id == sample_name)
        
        # now I can compute the NMI and the metrics and put them in the dataframe
        # Number of clusters tool and process
        n_raw_clones_tool = length(unique(table_sample$cluster_id_tool))
        n_true_driver_clones_process = length(unique(table_sample$cluster_id_process))
        
        # NMI and ARI
        nmi_raw = randnet::NMI(as.factor(table$cluster_id_tool),
                               as.factor(table$cluster_id_process))
        
        ari_raw = aricode::ARI(as.factor(table$cluster_id_tool),
                               as.factor(table$cluster_id_process))
        
        # Absolute error
        # Relative error
        
        df_samples = data.frame(spn = spn, purity = purity, 
                        coverage = coverage, vcf_caller = vcf_caller, 
                        cna_caller=cna_caller,
                        nmi_raw=nmi_raw,
                        nmi_interpreted=NA,
                        ari_interpreted=NA,
                        ari_raw=ari_raw,
                        n_raw_tool=n_raw_clones_tool,
                        n_true_driver_process = n_true_driver_clones_process,
                        n_pipeline_interpreted_tool = NA,
                        tool=tool,
                        n_pipeline_interpreted_tool_no_tail=NA,
                        n_true_driver_process_no_tail=NA)
        
        df = rbind(df, df_samples)
      }
    }
    else{
      table = table %>% 
        group_by(cluster_id_process, sample_id) %>%
        mutate(ccf_process=mean(ccf_process, na.rm=TRUE)) %>%
        ungroup() %>% 
        add_count(cluster_id_process, name="n_mutations_process") %>%
        add_count(cluster_id_tool, name="n_mutations_tool") %>%
        add_count(cluster_id_tool_interpreted, name="n_mutations_tool_interpreted")
      
      table = table %>% group_by(cluster_id_process) %>% 
        mutate(is_clonal_process=replace(FALSE, all(ccf_process > 0.95), TRUE)) %>% ungroup() %>% 
        mutate(cluster_id_process_full = cluster_id_process) %>% 
        mutate(cluster_id_process=replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal'))
      
      # unique(table$cluster_id_tool_interpreted)
    
      n_raw_clones_tool = length(unique(table$cluster_id_tool))
      # n_true_driver_clones_tool = length(unique(table$cluster_id_tool_interpreted_driver))
      n_interpreted_driver_clones_tool = length(unique(table$cluster_id_tool_interpreted))
      n_true_driver_clones_process = length(unique(table$cluster_id_process))
        
      fraction_raw_clusters = n_raw_clones_tool/n_true_driver_clones_process
      
      # fraction_true_driver_clusters = n_true_driver_clones_tool/n_true_driver_clones_process
      
      # fraction_false_raw_driver_clusters = (n_raw_clones_tool - n_true_driver_clones_tool)/n_raw_clones_tool
      
      # fraction_all_driver_clusters = n_interpreted_driver_clones_tool/n_true_driver_clones_process
      # fraction_false_driver_clusters = (n_interpreted_driver_clones_tool - n_true_driver_clones_tool)/n_interpreted_driver_clones_tool
    
      nmi_interpreted = randnet::NMI(as.factor(table$cluster_id_tool_interpreted),
                                     as.factor(table$cluster_id_process))
      
      ari_interpreted = aricode::ARI(as.factor(table$cluster_id_tool_interpreted),
                                     as.factor(table$cluster_id_process))
      
      nmi_raw = randnet::NMI(as.factor(table$cluster_id_tool),
                                     as.factor(table$cluster_id_process))
      
      ari_raw = aricode::ARI(as.factor(table$cluster_id_tool),
                                     as.factor(table$cluster_id_process))
      
      # Now we need to only consider clusters != 'Subclonal' in tool
      table_no_tail = table %>% filter(cluster_is_tool!='Subclonal')
      
      nmi_interpreted_no_tail = randnet::NMI(as.factor(table_no_tail$cluster_id_tool),
                             as.factor(table_no_tail$cluster_id_process))
      
      ari_interpreted_no_tail = aricode::ARI(as.factor(table_no_tail$cluster_id_tool_interpreted),
                                     as.factor(table_no_tail$cluster_id_process))
      
      n_interpreted_clones_tool_no_tail = length(unique(table_no_tail$cluster_id_tool_interpreted))
      n_true_driver_clones_process_no_tail = length(unique(table_no_tail$cluster_id_process))
      
      df = data.frame(spn = spn, purity = purity, 
                      coverage = coverage, vcf_caller = vcf_caller, 
                      cna_caller=cna_caller,
                      n_raw_tool=n_raw_clones_tool,
                      nmi_raw=nmi_raw,
                      nmi_interpreted=nmi_interpreted,
                      ari_interpreted=ari_interpreted,
                      ari_raw=ari_raw,
                      nmi_interpreted_no_tail=nmi_interpreted_no_tail,
                      ari_interpreted_no_tail=ari_interpreted_no_tail,
                      # n_true_driver_tool = n_true_driver_clones_tool,
                      n_true_driver_process = n_true_driver_clones_process,
                      n_pipeline_interpreted_tool = n_interpreted_driver_clones_tool,
                      n_pipeline_interpreted_tool_no_tail=n_interpreted_clones_tool_no_tail,
                      n_true_driver_process_no_tail=n_true_driver_clones_process_no_tail,
                      tool=tool)#,
      # fraction_raw_clusters=fraction_raw_clusters,
      # fraction_true_driver_clusters=fraction_true_driver_clusters),
      # fraction_false_driver_clusters=fraction_false_driver_clusters,
      # fraction_pipeline_driver_clusters=fraction_all_driver_clusters)
      
    }
    
    metrics_table = rbind(metrics_table, df)
  }
}
saveRDS(metrics_table, file.path(save_path, "metrics_tables/table_clusters_metrics.rds"))


# New plots ####

metrics_table = readRDS(file.path(save_path, "metrics_tables/table_clusters_metrics.rds"))
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
    # aes(shape=purity, group = metric_label, fill = coverage),
    position = position_dodge(width = 0.8),
    # alpha = 1,
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
  # scale_fill_manual(values=palette_coverage, name='Coverage')+
  my_theme

nmi_plot
ggsave(file.path(save_path, "plots/metrics/nmi_plot_tools.png"), nmi_plot)

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
    # aes(shape=purity, group = metric_label, fill = coverage),
    position = position_dodge(width = 0.8),
    # alpha = 1,
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
  # scale_fill_manual(values=palette_coverage, name='Coverage')+
  my_theme


ari_plot

ggsave(file.path(save_path, "plots/metrics/ari_plot.png"), ari_plot)

## Absolute error number of clusters ####
metrics_table_long_numbers = metrics_table %>%
  select(spn, purity, coverage,tool,
         n_raw_tool,
         n_pipeline_interpreted_tool,
         n_true_driver_process) %>%
  mutate(raw_error =  n_raw_tool - n_true_driver_process,
         interpreted_error = n_pipeline_interpreted_tool - n_true_driver_process) %>% 
  pivot_longer(
    cols = c(raw_error,
             interpreted_error),
    names_to = "error",
    values_to = "value"
  ) %>% 
  mutate(error_label=case_when(
    error=='raw_error'~"Raw",
    error=='interpreted_error'~"Interpreted"
  ))

mean_abs_err = metrics_table_long_numbers %>%
  group_by(error_label) %>%
  summarize(mean_value = mean(value)) %>%
  mutate(
    mean_label = paste0(error_label, " (", round(mean_value, 3), ")")
  )


clusters_error_plot = metrics_table_long_numbers %>%
  filter() %>% 
  mutate(
    error_label = factor(error_label, levels = c("Raw", "Interpreted")),
    purity = factor(purity),
    coverage=factor(coverage)
  ) %>%
  ggplot(aes(x = n_true_driver_process, y = value, color=spn, shape=purity)) +
  geom_point(aes(alpha=coverage), size = 3) +
  ggh4x::facet_nested(~"Error" + error_label) +
  geom_hline(
    data = mean_abs_err,
    aes(yintercept = mean_value, linetype = "Mean Error"),
    color = 'black',
    linetype = 'dashed',
    linewidth = 0.5
  ) +
  xlab('Number of ProCESS clones') +
  ylab('Error (#tool - #ProCESS)') +
  scale_color_manual(values=palette_spn, name='SPN')+
  scale_shape_manual(values=c(16,17,3), name='Purity')+
  scale_alpha_manual(values=c(0.2,0.5,1), name='Coverage')+
  scale_x_continuous(breaks = scales::pretty_breaks())+
  scale_y_continuous(breaks = scales::pretty_breaks())+
  my_theme+
  theme(legend.position = 'right')

clusters_error_plot

ggsave(file.path(save_path, "plots/metrics/clusters_error.png"), clusters_error_plot)



# clusters_error_plot2 = metrics_table_long_numbers %>%
#   mutate(
#     error_label = factor(error_label, levels = c("Raw", "Interpreted")),
#     purity = factor(purity) 
#   ) %>%
#   ggplot(aes(x = n_raw_tool, y = value, color=spn, shape=purity)) +
#   geom_point(size = 3) +
#   ggh4x::facet_nested(~"Error" + error_label) +
#   xlab('Number of raw clusters tool') +
#   ylab('Error (#ProCESS - #tool)') +
#   scale_color_manual(values=palette_spn, name='SPN')+
#   scale_shape_manual(values=c(16,17,3), name='Purity')+
#   scale_x_continuous(breaks = scales::pretty_breaks())+
#   scale_y_continuous(breaks = scales::pretty_breaks())+
#   my_theme+
#   theme(legend.position = 'right')
# 
# clusters_error_plot2
# 
# ggsave(file.path(save_path, "plots/metrics/clusters_error_raw_base.png"), clusters_error_plot2)

# clusters_error_plot3 = metrics_table_long_numbers %>%
#   mutate(
#     error_label = factor(error_label, levels = c("Raw", "Interpreted")),
#     purity = factor(purity) 
#   ) %>%
#   ggplot(aes(x = n_pipeline_interpreted_tool, y = value, color=spn, shape=purity)) +
#   geom_point(size = 3) +
#   ggh4x::facet_nested(~"Error" + error_label) +
#   xlab('Number of interpreted clusters tool') +
#   ylab('Error (#ProCESS - #tool)') +
#   scale_color_manual(values=palette_spn, name='SPN')+
#   scale_shape_manual(values=c(16,17,3), name='Purity')+
#   scale_x_continuous(breaks = scales::pretty_breaks())+
#   scale_y_continuous(breaks = scales::pretty_breaks())+
#   my_theme+
#   theme(legend.position = 'right')
# 
# clusters_error_plot3
# 
# ggsave(file.path(save_path, "plots/metrics/clusters_error_inter_base.png"), clusters_error_plot3)


## Relative error number of clusters ####
metrics_table_long_relative_error = metrics_table %>%
  select(spn, purity, coverage,tool,
         n_raw_tool,
         n_pipeline_interpreted_tool,
         n_true_driver_process) %>%
  mutate(relative_raw_error = abs((n_raw_tool - n_true_driver_process)/n_true_driver_process),
         relative_interpreted_error = abs((n_pipeline_interpreted_tool - n_true_driver_process)/n_true_driver_process)) %>% 
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
  xlab('Number of ProCESS clones') +
  ylab('Relative error |(#tool - #ProCESS)/#ProCESS|') +
  scale_color_manual(values=palette_spn, name='SPN')+
  scale_shape_manual(values=c(16,17,3), name='Purity')+
  scale_alpha_manual(values=c(0.2,0.5,1), name='Coverage')+
  scale_x_continuous(breaks = scales::pretty_breaks())+
  # scale_y_continuous(breaks = scales::pretty_breaks())+
  my_theme+
  theme(legend.position = 'right')

clusters_relative_error_scatter

ggsave(file.path(save_path, "plots/metrics/clusters_relative_error_scatter.png"), clusters_relative_error_scatter)

### Boxplot ####
relative_error_plot = metrics_table_long_relative_error %>%
  mutate(
    error_label = factor(error_label, levels = c("Raw", "Interpreted")),
    purity = factor(purity),
    coverage=factor(coverage),
    spn=factor(spn)
  ) %>%
  ggplot(aes(x = spn, y = value)) +
  geom_boxplot(aes(color = error_label), position = position_dodge(width = 0.8)) +
  geom_point(
    aes(shape=purity, group = error_label, alpha = coverage),
    # aes(shape=purity, group = metric_label, fill = coverage),
    position = position_dodge(width = 0.8),
    # alpha = 1,
    size = 2.5
  ) +
  geom_hline(
    data = mean_relative_err,
    aes(yintercept = mean_value, linetype = "Mean Error", color=error_label),
    linetype = 'dashed',
    linewidth = 0.5
  ) +
  xlab('SPN') +
  ylab('Relative error |(#tool - #ProCESS)/#ProCESS|') +
  scale_color_manual(values=palette_raw_int, name='Error') +
  scale_shape_manual(values=c(16,17,3), name='Purity') +
  scale_alpha_manual(values=c(0.2,0.5,1), name='Coverage') +
  my_theme+
  theme(legend.position = 'right')
  
relative_error_plot

ggsave(file.path(save_path, "plots/metrics/clusters_relative_error_boxplot.png"), relative_error_plot)


p = metrics_table_long_relative_error %>%
  mutate(
    error_label = factor(error_label, levels = c("Raw", "Interpreted")),
    purity = factor(purity),
    coverage=factor(coverage),
    spn=factor(spn),
    n_true_driver_process=factor(n_true_driver_process)
  ) %>%
  ggplot(aes(x = n_true_driver_process, y = value)) +
  geom_boxplot(aes(color = spn), position = position_dodge(width = 0.8)) +
  geom_point(
    aes(shape=purity, group = error_label, alpha = coverage),
    # aes(shape=purity, group = metric_label, fill = coverage),
    position = position_dodge(width = 0.8),
    # alpha = 1,
    size = 2.5
  ) +
  ggh4x::facet_nested(~"Error" + error_label) +
  geom_hline(
    data = mean_relative_err,
    aes(yintercept = mean_value, linetype = "Mean Error", color=error_label),
    linetype = 'dashed',
    linewidth = 0.5
  ) +
  xlab('Number of ProCESS clones') +
  ylab('Relative error |(#tool - #ProCESS)/#ProCESS|') +
  scale_color_manual(values=palette_spn, name='Error') +
  scale_shape_manual(values=c(16,17,3), name='Purity') +
  scale_alpha_manual(values=c(0.2,0.5,1), name='Coverage') +
  my_theme+
  theme(legend.position = 'right')

p


## Fractions number of clusters ####
metrics_table_long = metrics_table %>%
  select(spn, purity, coverage, cna_caller, vcf_caller,
         fraction_raw_clusters,
         fraction_pipeline_driver_clusters,
         n_true_driver_process) %>%
  pivot_longer(
    cols = c(fraction_raw_clusters, fraction_pipeline_driver_clusters),
    names_to = "metric",
    values_to = "value"
  ) %>% 
  mutate(metric_label=case_when(
    metric=='fraction_raw_clusters'~"Raw",
    metric=='fraction_pipeline_driver_clusters'~"Interpreted"
  ))

clusters_fraction_plot = metrics_table_long %>%
  mutate(
    metric_label = factor(metric_label, levels = c("Raw", "Interpreted")),
    purity = factor(purity) 
  ) %>%
  ggplot(aes(x = n_true_driver_process, y = value, color=spn, shape=purity)) +
  # geom_line(aes(group = purity), color='black')+
  # geom_line(aes(group = interaction(spn, purity)))+
  geom_point(size = 3) +
  ggh4x::facet_nested(~"Error" + metric_label) +
  xlab('Number of ProCESS clones') +
  ylab('#clusters tool/#clusters ProCESS') +
  scale_color_manual(values=palette_spn, name='SPN')+
  scale_shape_manual(values=c(16,17,3), name='Purity')+
  scale_x_continuous(breaks = scales::pretty_breaks())+
  scale_y_continuous(breaks = scales::pretty_breaks())+
  my_theme+
  theme(legend.position = 'right')

clusters_fraction_plot

ggsave(file.path(save_path, "plots/metrics/clusters_fraction_error.png"), clusters_fraction_plot)



# Patchwork ####

p = patchwork::wrap_plots(
  # clusters_error_plot + labs(title="Number of clusters error"),
  clusters_relative_error_scatter+ labs(title="Number of clusters error"),
  ggplot() + labs(title="CCF"),
  nmi_plot + theme(legend.position = 'right') + labs(title="Clustering NMI"), 
  # ari_plot + labs(title="Clustering ARI"),
  ggplot() + labs(title="Signatures"),
  design="aabb\nccdd"
) &
  patchwork::plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(size=18, face="bold"))
p

ggsave(file.path(save_path, "plots/metrics/fig4.png"),p,
       width=35, height=25, units="cm")

# Old plots ####
metrics_table_long = metrics_table %>%
  select(spn, purity, coverage, cna_caller, vcf_caller,
         fraction_true_driver_clusters, 
         fraction_pipeline_driver_clusters, 
         fraction_false_driver_clusters) %>%
  pivot_longer(
    cols = c(fraction_true_driver_clusters, fraction_pipeline_driver_clusters, fraction_false_driver_clusters),
    names_to = "metric",
    values_to = "value"
  ) %>% 
  mutate(metric_label=case_when(
    metric=='fraction_true_driver_clusters'~"Interpreted clusters with T drivers",
    metric=='fraction_pipeline_driver_clusters'~"Interpreted clusters with T/F drivers",
    metric=='fraction_false_driver_clusters'~"Interpreted clusters with only F drivers"
  ))


p = metrics_table_long %>% filter(metric %in% c('fraction_true_driver_clusters', 'fraction_pipeline_driver_clusters')) %>%
  ggplot(aes(x = spn, y = value, fill = metric_label)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    x = "SPN",
    y = "#clusters tool/#clusters process",
    fill = "Metric"
  ) + 
  # theme_minimal() +
  theme(legend.position="bottom")+
  ggh4x::facet_nested(~"Purity" + purity)

p
ggsave(file.path(save_path, "plots/metrics/figure1.png"), p,
       width=25, height=18, units="cm")


p=metrics_table_long %>% filter(metric %in% c('fraction_false_driver_clusters')) %>% 
  ggplot(aes(x = spn, y = value, fill = metric_label)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    x = "SPN",
    y = "% of FP interpreted clusters", # numero diverso di cluster dovuti a falsi driver
    fill = "Metric"
  ) + 
  # facet_grid(~metric)+
  # scale_x_discrete(limits=spn_list) +
  # theme_minimal() +
  theme(legend.position="bottom")+
  # facet_grid(~purity)
  ggh4x::facet_nested(~"Purity" + purity)

p
ggsave(file.path(save_path, "plots/metrics/metrics_false_facet.png"), p)



metrics_table_long_numbers = metrics_table %>%
  select(spn, purity, coverage,n_true_driver_tool,
         n_true_driver_process,
         n_pipeline_interpreted_tool) %>%
  pivot_longer(
    cols = c(n_true_driver_tool, n_true_driver_process,
             n_pipeline_interpreted_tool),
    names_to = "n_clusters",
    values_to = "value"
  ) %>% 
  mutate(n_clusters_label=case_when(
    n_clusters=='n_true_driver_tool'~"Interpreted clusters with T drivers (tool)",
    n_clusters=='n_pipeline_interpreted_tool'~"Interpreted clusters with T/F drivers (tool)",
    n_clusters=='n_true_driver_process'~"Interpreted clusters with T drivers (process)"
  ))

p=metrics_table_long_numbers %>% 
  ggplot(aes(x = spn, y = value, fill = n_clusters_label)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    x = "SPN",
    y = "#clusters", # numero diverso di cluster dovuti a falsi driver
    fill = ""
  ) + 
  # facet_grid(~purity)+
  # scale_x_discrete(limits=spn_list) +
  # theme_minimal() +
  ggh4x::facet_nested(~"Purity" + purity)+
  theme(legend.position="bottom")

p
ggsave(file.path(save_path, "plots/metrics/number_clusters_w_purity.png"), p)




# Boxplots ####
metrics_table_long %>% filter(metric %in% c('fraction_false_driver_clusters')) %>% 
  ggplot(aes(x = spn, y = value, fill = metric_label)) +
  geom_boxplot() +
  ggh4x::facet_nested(~"Purity" + purity)+
  # ggh4x::facet_nested("Coverage" + coverage~"Purity" + purity)+
  ylim(0,0.8)




