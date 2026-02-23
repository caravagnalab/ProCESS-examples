library(ggplot2)
library(tidyr)
library(dplyr)
library(ProCESS)

spn = 'SPN03'
purity=0.3
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
    
    # final_table %>% filter((is_driver_tool==T) | (is_driver_process==T)) %>%
    #   select(mutation_id, cluster_id_tool, 
    #          cluster_id_process, is_driver_tool, is_driver_process,
    #          driver_label_tool, driver_label_process) %>% View()
    # 
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
    
    ### Drivers analysis ####
    #### TP ####
    # TP: clusters with at least one true driver (i.e. at least one is_driver_process == TRUE)
    
    TP_d = drivers_tool %>% 
      group_by(cluster_id_tool) %>% 
      summarise(has_process_driver = any(is_driver_process == TRUE, na.rm = TRUE)) %>% 
      summarise(count = sum(has_process_driver))
    
    which_TP_d = drivers_tool %>% 
      group_by(cluster_id_tool) %>% 
      summarise(has_process_driver = any(is_driver_process == TRUE, na.rm = TRUE)) %>% 
      filter(has_process_driver)%>% pull(cluster_id_tool) 
    
    #### FP ####
    # FP: cluster with wrong drivers (and without true drivers)
      # tool clusters where all rows have is_driver_process == FALSE (i.e., zero TRUE in that cluster).
    FP_d = drivers_tool %>% 
      group_by(cluster_id_tool) %>% 
      summarise(all_false = all(is_driver_process == FALSE, na.rm = TRUE)) %>% 
      summarise(count = sum(all_false)) 
    
    which_FP_d = drivers_tool %>% 
      group_by(cluster_id_tool) %>% 
      summarise(all_false = all(is_driver_process == FALSE, na.rm = TRUE)) %>% 
      filter(all_false) %>% pull(cluster_id_tool) 
    
    #### FN ####
    d_interpreted_drivers = setdiff(unique(final_table$cluster_id_tool_interpreted_driver), "Subclonal")
    which_FN_d = setdiff(d_interpreted_drivers, unique(drivers_tool$cluster_id_tool)) 
    FN_d = length(which_FN_d)
    
    ### Clusters analysis ####
    
    #### TP #### 
      # clusters which are both in cluster_id_tool_interpreted and in cluster_id_tool_interpreted_drivers
    c_interpreted_drivers = setdiff(unique(final_table$cluster_id_tool_interpreted_driver), "Subclonal")
    which_TP_c = intersect(c_interpreted_drivers, unique(final_table$cluster_id_tool_interpreted))
    TP_c = length(which_TP_c)
    
    #### FP ####
      # clusters which are in cluster_id_tool_interpreted but not in cluster_id_tool_interpreted_drivers
    which_FP_c = setdiff(unique(final_table$cluster_id_tool_interpreted), c_interpreted_drivers)
    FP_c = length(which_FP_c)
    
    #### FN ####
    # setdiff(x, y) restituisce gli elementi presenti nel vettore x ma non nel vettore y
    which_FN_c = setdiff(c_interpreted_drivers, unique(drivers_tool$cluster_id_tool)) 
    FN_c = length(which_FN_c)
    
    # FN: process clusters with drivers but not found in tool (i.e. it's like FP but with opposite clusters)
      # so I need to take those clusters that do not have any true driver (i.e. are not in interpreted)
    
    # FN = drivers_process %>% 
    #   group_by(cluster_id_process) %>% 
    #   summarise(all_false = all(is_driver_tool == FALSE, na.rm = TRUE)) %>% 
    #   summarise(count = sum(all_false))
    # 
    # which_FN = drivers_process %>% 
    #   group_by(cluster_id_process) %>% 
    #   summarise(all_false = all(is_driver_tool == FALSE, na.rm = TRUE)) %>% 
    #   filter(all_false) %>% pull(cluster_id_process) 
    
    df <- data.frame(
      spn = spn, purity = purity, coverage = coverage,
      vcf_caller = vcf_caller, cna_caller = cna_caller,
      tool = tool,
      TP_d = TP_d$count,
      TP_d_list = I(list(unname(which_TP_d))),
      FP_d = FP_d$count,
      FP_d_list = I(list(unname(which_FP_d))),
      FN_d = FN_d,
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
        precision_d = ifelse(TP_d + FP_d > 0, TP_d / (TP_d + FP_d), NA_real_),
        recall_d    = ifelse(TP_d + FN_d > 0, TP_d / (TP_d + FN_d), NA_real_)
      )
      
    metrics_drivers = bind_rows(metrics_drivers, df)
  }
}
saveRDS(metrics_drivers, file.path(save_path, "metrics_tables/metrics_drivers.rds"))

metrics_drivers= readRDS(file.path(save_path, "metrics_tables/metrics_drivers.rds"))
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

palette = RColorBrewer::brewer.pal(n=3, name="Dark2")
palette_spn = RColorBrewer::brewer.pal(n=length(spn_list), name="Set1")
palette_spn[6] = 'gold'

# Plots ####
t = metrics_drivers %>%
  select(spn, purity, coverage, cna_caller, vcf_caller,
         tool,
         TP_d, FP_d, FN_d,
         TP_c, FP_c, FN_c,
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

### TP, FP, FN with tool ####
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
                      ~"Error"+~factor(type, levels=c("Driver", "Cluster"))
  ) +
  xlab("")+
  ylab("")+
  scale_color_manual(values=palette_spn, name='SPN')+
  my_theme

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
  ggh4x::facet_nested("Error"~factor(type, levels=c("Driver", "Cluster"))) +
  xlab("")+
  ylab("")+
  scale_color_manual(values=palette_spn, name='SPN')+
  my_theme

## Precision vs recall ####
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
  scale_color_manual(values=palette_spn, name='SPN')+
  scale_alpha_manual(values = c(0.1, 0.5,1))+
  my_theme


## Precision and recall ####

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
  # ggh4x::facet_nested("Error"~factor(type, levels=c("Driver", "Cluster"))) +
  ggh4x::facet_nested("Tool"+~factor(tool, levels=c('viber', 'pyclonevi'))
                      ~"Error"+~factor(type, levels=c("Driver", "Cluster"))
                      ) +
  xlab("")+
  ylab("")+
  ylim(0,1)+
  scale_color_manual(values=palette_spn, name='SPN')+
  my_theme
