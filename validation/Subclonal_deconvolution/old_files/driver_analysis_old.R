library(ggplot2)
library(tidyr)
library(dplyr)
library(ProCESS)
library(stringr)

spn = 'SPN03'
purity=0.9
vcf_caller = "mutect2"
cna_caller = "ascat"
coverage = 100

coverage_list = c(50,100,150)
purity_list = c(0.3, 0.6, 0.9)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN06', 'SPN07')

combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller=cna_caller_list,
                    spn=spn_list)


tool = 'viber'
tools = c('viber', 'pyclonevi')

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

metrics_drivers = data.frame()

for(tool in tools){
  print(tool)
  for(i in 1:nrow(combs)){
    
    coverage = combs[i, "coverage"]
    purity = combs[i, "purity"]
    vcf_caller = combs[i, "vcf_caller"]
    cna_caller = combs[i, "cna_caller"]
    spn = combs[i, "spn"]
    
    simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
    
    print(paste0(spn,'_', simulation_id))
    
    # Get interpreted table
    final_table = readRDS(file.path(main_path, "validation_subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds"))) # process table in folder tables/
    
    
    # Compute TP drivers (i.e. how many drivers found in viber are true drivers)
    # TP: cluster con almeno un driver giusto
    
    drivers_tool = final_table %>% 
      filter(is_driver_tool == TRUE) %>% 
      select(mutation_id, cluster_id_tool, cluster_id_process, 
             is_driver_tool, is_driver_process) %>% 
      distinct(mutation_id, cluster_id_tool, 
               is_driver_process, cluster_id_process)
    
    drivers_process = final_table %>% 
      filter(is_driver_process == TRUE) %>% 
      select(mutation_id, cluster_id_tool, cluster_id_process, 
             is_driver_tool, is_driver_process) %>% 
      distinct(mutation_id, cluster_id_tool, 
               is_driver_tool, cluster_id_process)
    
    ## Clusters analysis ####
    
    ### Interpreted ####
    #### TP #### 
    # clusters which are both in cluster_id_tool_interpreted and in cluster_id_tool_interpreted_drivers
    c_interpreted_drivers = setdiff(unique(final_table$cluster_id_tool_interpreted_driver), "Subclonal")
    
    which_TP_c = intersect(c_interpreted_drivers, unique(final_table$cluster_id_tool_interpreted))
    TP_c = length(which_TP_c)
    
    #### FP ####
    # clusters which are in cluster_id_tool_interpreted but not in cluster_id_tool_interpreted_drivers
    which_FP_c = setdiff(unique(final_table$cluster_id_tool_interpreted), c_interpreted_drivers)
    which_FP_c = setdiff(which_FP_c, "Subclonal")
    FP_c = length(which_FP_c)
    
    #### FN ####
    # setdiff(x, y) restituisce gli elementi presenti nel vettore x ma non nel vettore y
    which_FN_c = setdiff(c_interpreted_drivers, unique(drivers_tool$cluster_id_tool)) 
    FN_c = length(which_FN_c)
    
    ### Blind ####
    #### TP #### 
    # clusters which are both in cluster_id_tool_interpreted and in cluster_id_tool_interpreted_drivers
    c_interpreted_drivers = setdiff(unique(final_table$cluster_id_tool_interpreted_driver), "Subclonal")
    which_TP_c_blind = intersect(c_interpreted_drivers, unique(final_table$cluster_id_tool))
    TP_c_blind = length(which_TP_c_blind)
    
    #### FP ####
    # clusters which are in cluster_id_tool_interpreted but not in cluster_id_tool_interpreted_drivers
    which_FP_c_blind = setdiff(unique(final_table$cluster_id_tool), c_interpreted_drivers)
    FP_c_blind = length(which_FP_c_blind)
    
    #### FN ####
    FN_c_blind = 0
    
    
    ## Drivers analysis ####
      # Only consider TP clusters
    #### TP ####
    TP_table = final_table %>% 
      dplyr::filter(cluster_id_tool_interpreted %in% which_TP_c)
    
    drivers_TP_tool = TP_table %>% 
      filter(is_driver_tool == TRUE) %>% 
      select(mutation_id, cluster_id_tool, cluster_id_process, 
             is_driver_tool, is_driver_process) %>% 
      distinct(mutation_id, cluster_id_tool, 
               is_driver_process, cluster_id_process)
    
    drivers_TP_process = TP_table %>% 
      filter(is_driver_process == TRUE) %>% 
      select(mutation_id, cluster_id_tool, cluster_id_process, 
             is_driver_tool, is_driver_process) %>% 
      distinct(mutation_id, cluster_id_tool, 
               is_driver_tool, cluster_id_process)
    
    # TP: how many true drivers I have in the TP clusters
    
    TP_d = drivers_TP_tool %>%
      summarise(count = sum(is_driver_process, na.rm = TRUE))
    
    which_TP_d = drivers_TP_tool %>%
      filter(is_driver_process) %>%
      pull(cluster_id_tool)
    
    #### FP ####
    # FP: how many wrong drivers I have in the TP clusters (i.e. in the interpreted clusters)
      # tool clusters where all rows have is_driver_process == FALSE (i.e., zero TRUE in that cluster).
    FP_d = drivers_TP_tool %>%
      summarise(count = sum(!is_driver_process, na.rm = TRUE))
    
    which_FP_d = drivers_TP_tool %>%
      filter(!is_driver_process) %>%
      pull(cluster_id_tool)
    
    #### FN ####
    FN_d = drivers_TP_process %>%
      summarise(count = sum(!is_driver_tool, na.rm = TRUE))
    
    which_FN_d = drivers_TP_process %>%
      filter(!is_driver_tool) %>%
      pull(cluster_id_tool)
    
    
    df <- data.frame(
      spn = spn, purity = purity, coverage = coverage,
      vcf_caller = vcf_caller, cna_caller = cna_caller,
      tool = tool,
      
      TP_c_blind = TP_c_blind,
      FP_c_blind = FP_c_blind,
      FN_c_blind = FN_c_blind,
      
      TP_d = TP_d$count,
      TP_d_list = I(list(unname(which_TP_d))),
      FP_d = FP_d$count,
      FP_d_list = I(list(unname(which_FP_d))),
      FN_d = FN_d$count,
      FN_d_list = I(list(unname(which_FN_d))),
      
      TP_c = TP_c,
      TP_c_list = I(list(unname(which_TP_c))),
      FP_c = FP_c,
      FP_c_list = I(list(unname(which_FP_c))),
      FN_c = FN_c,
      FN_c_list = I(list(unname(which_FN_c))),
      stringsAsFactors = FALSE
    )
    
    df = df %>% 
      mutate(
        precision_c = ifelse(TP_c + FP_c > 0, TP_c / (TP_c + FP_c), NA_real_),
        recall_c    = ifelse(TP_c + FN_c > 0, TP_c / (TP_c + FN_c), NA_real_),
        precision_c_blind = ifelse(TP_c_blind + FP_c_blind > 0, TP_c_blind / (TP_c_blind + FP_c_blind), NA_real_),
        recall_c_blind  = 1,
        precision_d = ifelse(TP_d + FP_d > 0, TP_d / (TP_d + FP_d), NA_real_),
        recall_d    = ifelse(TP_d + FN_d > 0, TP_d / (TP_d + FN_d), NA_real_)
      )
      
    metrics_drivers = bind_rows(metrics_drivers, df)
  }
}
saveRDS(metrics_drivers, file.path(save_path, "metrics_tables/metrics_drivers.rds"))

metrics_drivers= readRDS(file.path(save_path, "metrics_tables/metrics_drivers.rds"))

# Import themes ####
source('/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/figures/figure3/utils_plot.R')

# palette = RColorBrewer::brewer.pal(n=3, name="Dark2")
# SPN_colors = RColorBrewer::brewer.pal(n=length(spn_list), name="Set1")
# SPN_colors[6] = 'gold'
SPN_colors <-c("01"='steelblue', "02"='seagreen', "03"='goldenrod', 
               "04"='coral', "05"="magenta4","06"='palevioletred', "07"='indianred3')

# Plots only interpreted ####
t = metrics_drivers %>%
  select(spn, purity, coverage, cna_caller, vcf_caller,
         tool,
         TP_d, FP_d, FN_d,
         TP_c, FP_c, FN_c,
         #TP_c_blind, FP_c_blind, FN_c_blind,
         precision_c, recall_c, 
         precision_d, recall_d) %>%
  mutate(spn = as.character(spn),
         spn = str_extract(spn, "\\d+"),
         spn = str_sub(spn, -2, -1),
         spn=factor(spn)) %>%
  pivot_longer(
    cols = matches("_[cd]$"),
    names_to = c("metric", "type"),
    names_pattern = "^(.*)_([cd])$",
    values_to = "value"
  ) %>% 
  mutate(type = ifelse(type=='d', 'Driver', 'Cluster'))

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

### Precision blind vs interpreted ####

p = t1 %>%
  filter(metric %in% c('precision')) %>%
  mutate(metric = factor(metric,
                         levels = c("precision")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("viber", "pyclonevi")),
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
    # ggh4x::facet_nested(~"Tool" + ~factor(tool, levels=c('viber', 'pyclonevi'))) +
    xlab("Cluster Filtering")+
    ylab("Mutation Clustering Precision")+
    scale_color_manual(values=SPN_colors, name='SPN')+
    my_ggplot_theme()

ggsave(filename = paste0(save_path, "plots/metrics/cluster_analysis/precision_blind_vs_interpreted.pdf"), 
       plot = p, device="pdf", width=4, height=3, units="in")

saveRDS(p, file = paste0(save_path, "plots/metrics/cluster_analysis/rds/precision_blind_vs_interpreted.rds"))

### Precision blind vs interpreted with tool ####
p = t1 %>%
  filter(metric %in% c('precision')) %>%
  mutate(metric = factor(metric,
                         levels = c("precision")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("viber", "pyclonevi")),
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
  ggh4x::facet_nested(~"Tool" + ~factor(tool, levels=c('viber', 'pyclonevi'))) +
  xlab("Cluster Filtering")+
  ylab("Mutation Clustering Precision")+
  scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()

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
         tool=factor(tool,
                     levels = c("viber", "pyclonevi")),
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

ggsave(filename = paste0(save_path, "plots/metrics/cluster_analysis/precision_recall_blind_vs_interpreted.pdf"), 
       plot = p, device="pdf", width=4, height=3, units="in")

saveRDS(p, file = paste0(save_path, "plots/metrics/cluster_analysis/rds/precision_recall_blind_vs_interpreted.rds"))

### Precision blind vs interpreted with tool ####
p = t1 %>%
  filter(metric %in% c('precision', 'recall')) %>%
  mutate(metric = factor(metric,
                         levels = c("precision", "recall")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("viber", "pyclonevi")),
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
    ~"Tool" + ~factor(tool, levels=c('viber', 'pyclonevi'))) +
  xlab("Cluster Filtering")+
  # ylab("Mutation Clustering Precision")+
  scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()
p
ggsave(filename = paste0(save_path, "plots/metrics/cluster_analysis/precision_recall_blind_vs_interpreted_w_tool.pdf"), 
       plot = p, device="pdf", width=4, height=3, units="in")
saveRDS(p, file = paste0(save_path, "plots/metrics/cluster_analysis/rds/precision_recall_blind_vs_interpreted_w_tool.rds"))


# OLD plots ####
### TP, FP, FN with tools, driver (only on TP clusters), clusters ####
t %>%
  filter(metric %in% c('TP', 'FP', 'FN')) %>% 
  mutate(metric = factor(metric,
                         levels = c("TP", "FP", "FN")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("viber", "pyclonevi")),
         type=factor(type,
                     levels = c("Driver", "Cluster"))) %>%
  ggplot(aes(x = metric, y = value)) +
  geom_boxplot(outliers = F) +
  stat_summary(
    aes(x = metric, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.4)
  )+
  ggh4x::facet_nested("Tool"+~factor(tool, levels=c('viber', 'pyclonevi'))
                      ~"Type"+~factor(type, levels=c("Driver", "Cluster"))
  ) +
  xlab("")+
  ylab("")+
  scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()


### TP, FP, FN without tool ####
t %>%
  filter(metric %in% c('TP', 'FP', 'FN')) %>% 
  mutate(metric = factor(metric,
                         levels = c("TP", "FP", "FN")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("viber", "pyclonevi")),
         type=factor(type,
                     levels = c("Driver", "Cluster"))) %>%
  ggplot(aes(x = metric, y = value)) +
  geom_boxplot(outliers = F) +
  stat_summary(
    aes(x = metric, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.4)
  )+
  ggh4x::facet_nested(~factor(type, levels=c("Driver", "Cluster"))) +
  xlab("")+
  ylab("")+
  scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()


### TP, FP, FN without tool only drivers ####
t %>%
  filter(metric %in% c('TP', 'FP', 'FN'), type == 'Driver') %>% 
  mutate(metric = factor(metric,
                         levels = c("TP", "FP", "FN")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("viber", "pyclonevi"))) %>%
  ggplot(aes(x = metric, y = value)) +
  geom_boxplot(outliers = F) +
  stat_summary(
    aes(x = metric, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.4)
  )+
  # ggh4x::facet_nested(~factor(type, levels=c("Driver", "Cluster"))) +
  xlab("Driver analysis")+
  ylab("Value")+
  scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()

### Precision vs recall ####
pr = t %>%
  filter(metric %in% c("TP", "FP", "FN")) %>%
  pivot_wider(
    names_from = metric,
    values_from = value
  ) %>%
  mutate(
    precision = TP / (TP + FP),
    recall    = TP / (TP + FN)
  ) %>% 
  group_by(spn, purity,type) %>% 
  mutate(mean_precision = mean(precision, na.rm = TRUE),
         mean_recall = mean(recall, na.rm = TRUE)) %>%
  ungroup()%>%
  distinct(
    spn, purity, tool, type, cna_caller, vcf_caller, # remove differences per coverage because I only want 1 point in the scatter
    mean_precision, mean_recall,
    .keep_all = TRUE
  )

pr %>% 
  mutate(purity = factor(purity)) %>% 
  ggplot(aes(x = mean_precision, y = mean_recall, color = spn, alpha = purity)) +
  geom_point(size=3,na.rm = TRUE) +
  ggh4x::facet_nested(
    # "Tool" ~ tool ~ "Error" + ~ type
    ~factor(type, levels=c("Driver", "Cluster"))
  ) +
  theme_minimal()+
  xlim(0,1)+
  ylim(0,1)+
  xlab('Precision')+
  ylab('Recall')+
  scale_color_manual(values=SPN_colors, name='SPN')+
  scale_alpha_manual(values = c(0.1, 0.5,1))+
  my_ggplot_theme()


### Precision and recall ####

t %>%
  filter(metric %in% c('precision', 'recall')) %>%
  mutate(metric=ifelse(metric=='precision', 'Precision', 'Recall'))  %>%
  mutate(metric = factor(metric,
                         levels = c("Precision", "Recall")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("viber", "pyclonevi"))) %>% 
  ggplot(aes(x = metric, y = value)) +
  geom_boxplot(outliers = F) +
  stat_summary(
    aes(x = metric, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3)
  ) +
  ggh4x::facet_nested(~factor(type, levels=c("Driver", "Cluster"))) +
  # ggh4x::facet_nested("Tool"+~factor(tool, levels=c('viber', 'pyclonevi'))
  #                     ~"Error"+~factor(type, levels=c("Driver", "Cluster"))
  #                     ) +
  xlab("")+
  ylab("")+
  ylim(0,1)+
  scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()

# Plots blind vs interpreted ####

### TP, FP, FN without tool ####
t1 %>%
  filter(metric %in% c('TP', 'FP', 'FN')) %>% 
  mutate(metric = factor(metric,
                         levels = c("TP", "FP", "FN")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("viber", "pyclonevi")),
         type=factor(type,
                     levels = c('Cluster interpreted', 'Cluster blind'))) %>%
  ggplot(aes(x = metric, y = value)) +
  geom_boxplot(outliers = F) +
  stat_summary(
    aes(x = metric, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.4)
  )+
  ggh4x::facet_nested(~factor(type, levels=c('Cluster interpreted', 'Cluster blind'))) +
  xlab("")+
  ylab("")+
  scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()

## Precision and recall ####
t1 %>%
  filter(metric %in% c('precision', 'recall')) %>%
  mutate(metric=ifelse(metric=='precision', 'Precision', 'Recall'))  %>%
  mutate(metric = factor(metric,
                         levels = c("Precision", "Recall")),
         purity = factor(purity),
         coverage = factor(coverage),
         tool=factor(tool,
                     levels = c("viber", "pyclonevi"))) %>% 
  ggplot(aes(x = metric, y = value)) +
  geom_boxplot(outliers = F) +
  stat_summary(
    aes(x = metric, color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3)
  ) +
  ggh4x::facet_nested("Error"~factor(type, levels=c('Cluster blind','Cluster interpreted'))) +
  xlab("")+
  ylab("")+
  ylim(0,1)+
  scale_color_manual(values=SPN_colors, name='SPN')+
  my_ggplot_theme()


