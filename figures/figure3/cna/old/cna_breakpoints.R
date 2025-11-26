rm(list=ls())
library(ProCESS)
library(ggplot2)
library(dplyr)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
# source("compute_FGA.R")
scout_dir <-"/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
SPNS <- c("SPN01","SPN02","SPN03","SPN04", "SPN06", 'SPN07')
COVERAGES <- c("50","100", "150")
PURITIES <- c("0.3","0.6","0.9")
SPN_colors <-c("SPN01"='steelblue', "SPN02"='seagreen', "SPN03"='goldenrod', 
               "SPN04"='coral', "SPN06"='palevioletred', "SPN07"='indianred3')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")
WGD_samples <- c("SPN01_1.1","SPN01_1.3","SPN06_3.1","SPN06_3.2")

validation_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION"
params_grid = expand.grid(COVERAGES, PURITIES)
colnames(params_grid) = c("coverage", "purity")
df_all_SPN_cna <- list()
for (SPN in SPNS){
  spn <- SPN
  df = lapply(1:nrow(params_grid), function(i) {
    coverage = params_grid[i,]$coverage
    purity = params_grid[i,]$purity
    comb <- paste0(coverage,"x_",purity)
    samples <- get_sample_names(spn = spn)
    validation_dir_cna <- file.path(validation_dir,spn,"cna")
    all_metrics_comb <- list()
    for (sample in samples){
      metrics_filename <- paste0(validation_dir_cna,"/",comb,"/",sample,"/metrics.rds")
      metrics_bp_filename <- paste0(validation_dir_cna,"/",comb,"/",sample,"/metrics_bp.rds")
      if (!file.exists(metrics_filename)){
        message("File not found: ", metrics_filename)
        # all_metrics_comb[[sample]] <- NA
      } else{
        metrics_df <- readRDS(metrics_filename) %>% 
          mutate(delta_purity=as.numeric(true_purity)-as.numeric(purity)) %>% 
          mutate(delta_ploidy=as.numeric(true_ploidy)-as.numeric(ploidy)) %>% 
          filter(tool!="cnvkit")
        metrics_bp_df <- readRDS(metrics_bp_filename) %>% 
          filter(chr=="genome") %>% 
          filter(tool!="cnvkit")
        all_metrics_comb[[sample]] <- inner_join(metrics_df,metrics_bp_df,by=c("tool","sample","spn","coverage","fga","fgs","true_purity"))
      }
    }
    all_combinations_SPN <- do.call("rbind",all_metrics_comb)
  })
  df_all_SPN_cna[[spn]] <- do.call("rbind",df)  
}
df_all_combs_SPN_cna <- do.call("rbind",df_all_SPN_cna)
df_all_combs_SPN_cna <- df_all_combs_SPN_cna %>% 
  mutate(fga_class=case_when(fga>=30 ~ "High FGA",
                             TRUE ~ "Low FGA"))
df_all_combs_SPN_cna <- df_all_combs_SPN_cna %>% mutate(tool = case_when(
  tool == 'ascat' ~ 'ASCAT',
  tool == 'sequenza' ~ 'Sequenza',
  tool == 'battenberg' ~ 'Battenberg',
)) 



plt <- df_all_combs_SPN %>%
  ggplot() +
  geom_boxplot(aes(x = tool, y = correctness_clonal*100, fill = tool)) +
  scale_fill_manual('Caller',values = c('deepskyblue4', 'maroon', 'sienna')) +
  theme_minimal() +
  xlab('') +
  ylab('% correcly identified clonal CN') +
  
  df_all_combs_SPN %>%
  ggplot() +
  geom_boxplot(aes(x = tool, y = precision, fill = tool)) +
  scale_fill_manual('Caller', values = c('deepskyblue4', 'maroon', 'sienna')) +
  theme_minimal()  +
  #facet_grid(true_purity ~ coverage) +
  ylim(0,1) +
  ylab('Breakpoint precision') +
  xlab('') +
  
  df_all_combs_SPN %>%
  ggplot() +
  geom_boxplot(aes(x = tool, y = recall, fill = tool)) +
  scale_fill_manual('Caller',values = c('deepskyblue4', 'maroon', 'sienna')) +
  theme_minimal() +
  ylim(0,1) +
  #facet_grid(true_purity ~ coverage) +
  ylab('Breakpoint sensitivity') +
  xlab('') +
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom')
ggsave(filename = "cna_stat.pdf",plot = plt, width = 9,height = 3)


p1 <- df_all_combs_SPN %>% 
  group_by(tool, fga_class, coverage, true_purity) %>%
  summarise(
    mean_precision = mean(precision, na.rm = TRUE),
    sd_precision   = sd(precision, na.rm = TRUE),
    mean_recall    = mean(recall, na.rm = TRUE),
    sd_recall      = sd(recall, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(
    x = mean_recall,
    y = mean_precision
  )) +
  
  # First color scale: by tool for lines
  geom_line(linewidth = 1, aes(color = tool, group = interaction(tool, true_purity),alpha=true_purity)) +
  scale_color_manual(values = col_cna_tools, name = "Tool") +
  # scale_color_manual(values = col_cna_tools, name = "Tool") +
  ggnewscale::new_scale_color() +
  geom_point(
    aes(
      color = fga_class,
      # size = coverage
    ),
    alpha = 0.8,
    size = 2
  ) +  

  
  # geom_pointrange(
  #   aes(
  #     color = fga_class,
  #     size = coverage,
  #     ymin = mean_precision - sd_precision,
  #     ymax = mean_precision + sd_precision
  #   ),
  #   fatten = 2,    # makes the central point more visible
  #   alpha = 0.8
  # ) +
  
  scale_color_manual(
    name = "FGA class",
    values = c("High FGA" = "indianred3", "Low FGA" = "grey40")
  ) +
  scale_size_manual(
    name = "coverage",
    values = c("50" = 1, "100"=3,"150"=5)
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Mean recall",
    y = "Mean precision",
    caption = "Each point shows mean ± SD (dispersion) per combination"
  )
p1



p2 <- df_all_combs_SPN_cna %>% 
  group_by(tool, fga_class, coverage, true_purity) %>%
  summarise(
    mean_precision = mean(precision, na.rm = TRUE),
    sd_precision   = sd(precision, na.rm = TRUE),
    mean_recall    = mean(recall, na.rm = TRUE),
    sd_recall      = sd(recall, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(group_ellipse = interaction(fga_class, tool, drop = TRUE)) %>% 
  ggplot(aes(
    x = mean_recall,
    y = mean_precision
  )) +
  
  
  geom_point(
    aes(
      #shape = fga_class,
      fill = fga_class,
      alpha = true_purity,
      # size = coverage,
      color=fga_class
    ),size=3
  ) +
  # scale_fill_manual(values = col_cna_tools, name = "Tool") +
  scale_alpha_manual(values = c("0.3"=.5,"0.6"=.7,"0.9"=1), name = "Purity") +
  scale_fill_manual(values = c("High FGA" = "indianred2", "Low FGA"="dodgerblue3"), name = "FGA class") +
  scale_color_manual(values = c("High FGA" = "indianred2", "Low FGA"="dodgerblue3"), name = "FGA class") +
  theme_bw() +
  labs(
    x = "Mean recall",
    y = "Mean precision",
    caption = "Each point shows mean ± SD (dispersion) per combination"
  )+
  facet_wrap(~tool)
p2

# p2 <- df_all_combs_SPN %>% 
#   # group_by(tool, fga_class, coverage, true_purity) %>%
#   # summarise(
#   #   mean_precision = mean(precision, na.rm = TRUE),
#   #   sd_precision   = sd(precision, na.rm = TRUE),
#   #   mean_recall    = mean(recall, na.rm = TRUE),
#   #   sd_recall      = sd(recall, na.rm = TRUE),
#   #   .groups = "drop"
#   # ) %>%
#   mutate(group_ellipse = interaction(fga_class, tool, drop = TRUE)) %>% 
#   ggplot(aes(
#     x = recall,
#     y = precision
#   )) +
#   geom_density_2d(aes(colour = tool))+
#   scale_color_manual(values = col_cna_tools, name = "Tool") +
#   ggnewscale::new_scale_color() +
#   geom_point(
#     aes(
#       #shape = fga_class,
#       fill = fga_class,
#       alpha = true_purity,
#       #size = coverage,
#       color=fga_class
#     ),
#     # stroke = 3,
#     #    shape = 21,
#     size=2
#   ) +
#   
#   scale_alpha_manual(values = c("0.3"=.5,"0.6"=.7,"0.9"=1), name = "Purity") +
#   scale_fill_manual(values = c("High FGA" = "indianred2", "Low FGA" = "grey40"), name = "FGA class") +
#   scale_color_manual(values = c("High FGA" = "indianred2", "Low FGA" = "grey40"), name = "FGA class") +
#   theme_bw() +
#   labs(
#     x = "Mean recall",
#     y = "Mean precision",
#     caption = "Each point shows mean ± SD (dispersion) per combination"
#   )+
#   facet_grid(coverage~fga_class)
# p2
# 
# df_long <- df_all_combs_SPN %>%
#   pivot_longer(
#     cols = c(precision, recall),
#     names_to = "metric",
#     values_to = "metric_value"
#   ) %>%
#   group_by(tool, fga_class, true_purity, metric) %>%
#   summarise(
#     mean_value = mean(metric_value, na.rm = TRUE),
#     sd_value = sd(metric_value, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # Plot
# p_col <- ggplot(df_long, aes(x = fga_class, y = mean_value, fill = true_purity)) +
#   geom_col(position = position_dodge()) +
#   # geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
#   #               position = position_stack(vjust = 0.5),
#   #               width = 0.2, color = "black") +
#   facet_grid(metric~ tool) +
#   # scale_fill_manual(
#   #   name = "FGA class",
#   #   values = c("High FGA" = "indianred3", "Low FGA" = "grey40")
#   # ) +
#   scale_fill_manual(
#     name = "SImualted purity",
#     values = purity_colors
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "bottom",
#     strip.text = element_text(face = "bold")
#   ) +
#   labs(
#     x = "Tool",
#     y = "Mean value",
#     fill = "Metric",
#     title = "Precision & Recall stacked per tool, faceted by purity",
#     caption = "Stacked bars show mean precision and recall per tool within each purity level"
#   )
# 
# p_col
# 
# # plt_final/p1+plot_layout(guides="collect") & theme(legend.position = "bottom")
# 
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# 
# df_long <- df_all_combs_SPN %>%
#   pivot_longer(
#     cols = c(precision, recall),
#     names_to = "metric",
#     values_to = "metric_value"
#   ) %>%
#   group_by(tool, fga_class, true_purity, metric) %>%
#   summarise(
#     mean_value = mean(metric_value, na.rm = TRUE),
#     sd_value   = sd(metric_value, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # Plot with lines
# p_line <- ggplot(df_long, aes(
#   x = fga_class,
#   y = mean_value,
#   color = true_purity,
#   group = true_purity
# )) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 2) +
#   geom_errorbar(
#     aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
#     width = 0.2,
#     alpha = 0.7
#   ) +
#   facet_grid(metric ~ tool) +
#   scale_color_manual(
#     name = "Simulated purity",
#     values = purity_colors
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "bottom",
#     strip.text = element_text(face = "bold")
#   ) +
#   labs(
#     x = "FGA class",
#     y = "Mean value",
#     title = "Precision & Recall across FGA class, faceted by tool",
#     caption = "Lines show mean ± SD per purity level"
#   )
# 
# p_line
# 
