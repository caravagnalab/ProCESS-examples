# rm(list=ls())
library(ProCESS)
library(ggplot2)
library(dplyr)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
library(ggnewscale)
library(tidyr)

source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
# source("compute_FGA.R")
scout_dir <-"/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
SPNS <- c("SPN01","SPN02","SPN03","SPN04", "SPN06", 'SPN07')
COVERAGES <- c("50","100", "150")
PURITIES <- c("0.3","0.6","0.9")

source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
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
                             TRUE ~ "Low FGA")) %>% 
  mutate(tool = case_when(
    tool == 'ascat' ~ 'ASCAT',
    tool == 'sequenza' ~ 'Sequenza',
    tool == 'battenberg' ~ 'Battenberg',
  )) 



plt_breakpoint <- df_all_combs_SPN_cna %>% 
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
      fill = fga_class,
      alpha = true_purity,
      color=fga_class
    ),size=2
  ) +
  # scale_fill_manual(values = col_cna_tools, name = "Tool") +
  scale_alpha_manual(values = c("0.3"=.5,"0.6"=.7,"0.9"=1), name = "Purity") +
  scale_fill_manual(values = c("High FGA" = "indianred2", "Low FGA"="dodgerblue3"), name = "FGA class") +
  scale_color_manual(values = c("High FGA" = "indianred2", "Low FGA"="dodgerblue3"), name = "FGA class") +
  my_ggplot_theme()+
  labs(
    x = "Mean recall",
    y = "Mean precision",
    #caption = "Each point shows mean ± SD (dispersion) per combination"
  )+
  facet_wrap(~tool)
plt_breakpoint



df_all_combs_SPN_cna_long <- df_all_combs_SPN_cna %>%
  mutate(correct_purity_estimation_class = case_when(
    abs(delta_purity) < 0.2 ~ "correctly estimated",
    delta_purity > 0.2 ~ "underestimated",
    delta_purity < -0.2 ~ "overestimated"),
    correct_ploidy_estimation_class = case_when(
      abs(delta_ploidy) < 1 ~ "correctly estimated",
      delta_ploidy > 1 ~ "underestimated",
      delta_ploidy < -1 ~ "overestimated"),
    true_purity_label = paste0('Purity =', true_purity)
  ) %>%
  pivot_longer(
    cols = c(delta_purity, delta_ploidy),
    names_to = "delta_type",
    values_to = "delta_value"
  ) %>%
  mutate(
    correct_estimation_class = ifelse(delta_type == "delta_purity",
                                      correct_purity_estimation_class,
                                      correct_ploidy_estimation_class)
  ) %>% 
  mutate(delta_type=case_when(delta_type=="delta_ploidy" ~ "Ploidy",
                              delta_type=="delta_purity" ~ "Purity"))

plt_purity_ploidy <- ggplot(df_all_combs_SPN_cna_long, aes(
  x = reorder(tool, delta_value, \(x) sd(x, na.rm = TRUE)),
  y = delta_value,
  group = interaction(tool, true_purity_label)
)) +
  
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  # scale_color_manual(name = "CNA caller", values = col_cna_tools) +
  # scale_fill_manual(name = "CNA caller", values = col_cna_tools) +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  # Jitter points colored by FGA class
  geom_jitter(aes(color = fga_class, fill = fga_class, shape = correct_estimation_class),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              size = 2, alpha = 0.8) +
  scale_color_manual(name = "FGA class", values = c("High FGA"="indianred2","Low FGA"="dodgerblue3")) +
  scale_fill_manual(name = "FGA class", values = c("High FGA"="indianred2","Low FGA"="dodgerblue3")) +
  scale_shape_manual(name = "Estimation class",
                     values = c("correctly estimated"=16,
                                "underestimated"=25,
                                "overestimated"=24)) +
  
  # Horizontal lines
  geom_hline(aes(yintercept = 0), linetype = "solid", color = "grey40") +
  geom_hline(data = subset(df_all_combs_SPN_cna_long, delta_type=="Purity"),
             aes(yintercept = 0.2), linetype = "dashed", color = "grey") +
  geom_hline(data = subset(df_all_combs_SPN_cna_long, delta_type=="Purity"),
             aes(yintercept = -0.2), linetype = "dashed", color = "grey") +
  geom_hline(data = subset(df_all_combs_SPN_cna_long, delta_type=="Ploidy"),
             aes(yintercept = 1), linetype = "dashed", color = "grey") +
  geom_hline(data = subset(df_all_combs_SPN_cna_long, delta_type=="Ploidy"),
             aes(yintercept = -1), linetype = "dashed", color = "grey") +
  
  # Facet grid: delta_type as rows, true purity as columns
  facet_grid(delta_type ~ true_purity_label, scales = "free_y") +
  
  labs(y = "Deviation", x = "") + my_ggplot_theme()+
  theme(axis.text.x = element_text(angle = 30,hjust = 1),
        axis.title.x = element_blank())
plt_purity_ploidy
