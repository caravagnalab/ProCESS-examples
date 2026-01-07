library(ggplot2)
library(tidyverse)
library(patchwork)

base = "/orfeo/scratch/cdslab/shared/SCOUT/VALIDATION/"
SPNS <- paste0('SPN0', c(1,2,3,4,6,7))
CALLERS <- c('mutect2', 'strelka', 'freebayes')
COVERAGES <- c(50,100,150)
PURITY <- c(0.3, 0.6, 0.9)
MUTS <- c('SNV', 'INDEL')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
params_grid = expand.grid(COVERAGES, PURITY, MUTS, SPNS)
names(params_grid) <- c('cov', 'pur', 'muts', 'spn')

df_all_SPN_somatic = lapply(1:nrow(params_grid), function(i) {
  mut_type = params_grid[i,]$muts
  coverage = params_grid[i,]$cov
  purity = params_grid[i,]$pur
  spn = params_grid[i,]$spn

  file =  paste0('/orfeo/scratch/cdslab/shared/SCOUT/VALIDATION/', spn ,'/somatic/', coverage,'x_',purity,'p/report/',mut_type,'/metrics_new_binning.rds')

  if (file.exists(file)) {
    metrics = readRDS(file)

    parsed_metrics = lapply(names(metrics), function(caller_name) {
      caller_metrics = metrics[[caller_name]]
      sample_names = names(caller_metrics)
      lapply(sample_names, function(sample_id) {
        FP = caller_metrics[[sample_id]]$detection_summary["False Positive"]
        d = caller_metrics[[sample_id]]$performance_table %>%
          dplyr::mutate(sample_id = sample_id, SPN = spn) %>%
          dplyr::mutate(purity = as.numeric(purity), coverage = as.numeric(coverage))
        d$false_positive[1] = FP
        d
      }) %>% do.call("bind_rows", .) %>%
        dplyr::mutate(caller = caller_name)
    }) %>% do.call("bind_rows", .)

    dplyr::bind_cols(parsed_metrics, params_grid[i,])
  }
}) %>% do.call("bind_rows", .)

df_all_SPN_somatic <- df_all_SPN_somatic %>% 
  mutate(CCF_bin=case_when(CCF_bin=="Clonal"~"95-100%",
                           TRUE ~ CCF_bin)) %>% 
  mutate(CCF_bin=case_when(CCF_bin=="0-5%" ~ "0-0.05",
                           CCF_bin=="5-10%" ~ "0.05-0.10",
                           CCF_bin=="10-25%" ~ "0.10-0.25",
                           CCF_bin=="25-50%" ~ "0.25-0.50",
                           CCF_bin=="50-75%" ~ "0.50-0.75",
                           CCF_bin=="75-95%" ~ "0.75-0.95",
                           CCF_bin=="95-100%" ~ "0.95-1"))
median_df_tools <- df_all_SPN_somatic %>%
  # mutate(CCF_bin = factor(CCF_bin, levels=c("0-5%","5-10%","10-25%","25-50%","50-75%","75-95%","95-100%"))) %>%
  mutate(CCF_bin=factor(CCF_bin,levels=c("0-0.05","0.05-0.10","0.10-0.25","0.25-0.50","0.50-0.75","0.75-0.95","0.95-1"))) %>% 
  group_by(CCF_bin, caller, muts) %>%
  summarise(
    sensitivity = median(sensitivity, na.rm = TRUE),
    .groups = "drop"
  )

median_df_purity <- df_all_SPN_somatic %>%
  # mutate(CCF_bin = factor(CCF_bin, levels=c("0-5%","5-10%","10-25%","25-50%","50-75%","75-95%","95-100%"))) %>%
  mutate(CCF_bin=factor(CCF_bin,levels=c("0-0.05","0.05-0.10","0.10-0.25","0.25-0.50","0.50-0.75","0.75-0.95","0.95-1"))) %>% 
  group_by(CCF_bin, purity, muts) %>%
  summarise(
    sensitivity = median(sensitivity, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(purity=as.factor(purity))


median_df_coverage <- df_all_SPN_somatic %>%
  # mutate(CCF_bin = factor(CCF_bin, levels=c("0-5%","5-10%","10-25%","25-50%","50-75%","75-95%","95-100%"))) %>%
  mutate(CCF_bin=factor(CCF_bin,levels=c("0-0.05","0.05-0.10","0.10-0.25","0.25-0.50","0.50-0.75","0.75-0.95","0.95-1"))) %>% 
  group_by(CCF_bin, coverage, muts) %>%
  summarise(
    sensitivity = median(sensitivity, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(coverage=as.factor(coverage))

rect_df_ccf <- tibble::tibble(
  xmin = c("0-0.05","0.10-0.25","0.95-1"),
  xmax = c("0.05-0.10","0.75-0.95","0.95-1"),
  ymin = -0.15,
  ymax = -0.05,
  region = c(
    "Neutral mutations",
    "Subclonal mutations",
    "Clonal mutations")) %>%
  mutate(xmid = xmin)

rect_df_ccf_bin <- tibble::tibble(
  xmin = c("0-0.05","0.05-0.10","0.10-0.25","0.25-0.50","0.50-0.75","0.75-0.95","0.95-1"),
  xmax = c("0-0.05","0.05-0.10","0.10-0.25","0.25-0.50","0.50-0.75","0.75-0.95","0.95-1"),
  ymin = -0.15,
  ymax = -0.05,
  region = c("Neutral mutations","Neutral mutations",
             "Subclonal mutations","Subclonal mutations",
             "Subclonal mutations","Subclonal mutations",
             "Clonal mutations"),
  ccf_bin=c("0-0.05","0.05-0.10","0.10-0.25","0.25-0.50","0.50-0.75","0.75-0.95","0.95-1"))


#### plot tools
  
caller_line_boxes <- df_all_SPN_somatic %>% 
  mutate(
    CCF_bin = factor(
      CCF_bin,
      levels = c("0-0.05","0.05-0.10","0.10-0.25",
                 "0.25-0.50","0.50-0.75","0.75-0.95","0.95-1")
    )
  ) %>% 
  
  ggplot(aes(x = CCF_bin, y = sensitivity, col = as.factor(caller))) +
  geom_rect(
    data = rect_df_ccf_bin,
    aes(
      xmin = stage(xmin, after_scale = xmin-0.48),
      xmax = stage(xmax, after_scale = xmax+0.48),
      ymin = ymin, ymax = ymax
    ),
    inherit.aes = FALSE,
    alpha = 1,
    fill="white"
  ) +
  geom_rect(
    data = rect_df_ccf_bin,
    aes(
      xmin = stage(xmin, after_scale = xmin-0.48),
      xmax = stage(xmax, after_scale = xmax+0.48),
      ymin = ymin, ymax = ymax,
      fill = region
    ),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_text(
    data = rect_df_ccf_bin,
    aes(
      x = ccf_bin,
      y = (ymin + ymax) / 2,
      label = ccf_bin
    ),
    inherit.aes = FALSE,
    size = 3
  )+
  stat_summary(
    fun = "median",
    fun.min = function(x) boxplot.stats(x)$stats[2],
    fun.max = function(x) boxplot.stats(x)$stats[4],
    linewidth = 1,
    size = 0.5,
    position = position_dodge(width = 0.5)
  ) +
  geom_line(
    data = median_df_tools,
    aes(
      x = CCF_bin,
      y = sensitivity,
      group = caller,
      colour = caller
    ),
    linewidth = 0.3,
    position = position_dodge(width = 0.5)
  ) +
  scale_color_manual("Somatic caller", values = col_somatic_tools) +
  scale_fill_manual(
    name = "CCF class",
    values = c(
      "Neutral mutations"     = "grey40",
      "Subclonal mutations" = "#8c00ec",
      "Clonal mutations"         = "#ffae00",
      guide = guide_legend(nrow = 2, byrow = TRUE)
    )
  ) +
  ylab("Sensitivity") +
  xlab("CCF bin") +
  facet_wrap(~muts, ncol = 1,strip.position = "right") +
  my_ggplot_theme()+
  theme(
    axis.text.x = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0.05, "cm")
  )

caller_line_boxes
# 
# ### plot purity
# purity_line_boxes <- df_all_SPN_somatic %>% 
#   mutate(
#     CCF_bin = factor(
#       CCF_bin,
#       levels = c("0-0.05","0.05-0.10","0.10-0.25",
#                  "0.25-0.50","0.50-0.75","0.75-0.95","0.95-1")
#     )
#   ) %>% 
#   mutate(purity=as.factor(purity)) %>% 
#   ggplot(aes(x = CCF_bin, y = sensitivity, col = as.factor(purity))) +
#   geom_rect(
#     data = rect_df_ccf_bin,
#     aes(
#       xmin = stage(xmin, after_scale = xmin-0.48),
#       xmax = stage(xmax, after_scale = xmax+0.48),
#       ymin = ymin, ymax = ymax
#     ),
#     inherit.aes = FALSE,
#     alpha = 1,
#     fill="white"
#   ) +
#   geom_rect(
#     data = rect_df_ccf_bin,
#     aes(
#       xmin = stage(xmin, after_scale = xmin-0.48),
#       xmax = stage(xmax, after_scale = xmax+0.48),
#       ymin = ymin, ymax = ymax,
#       fill = region
#     ),
#     inherit.aes = FALSE,
#     alpha = 0.2
#   ) +
#   geom_text(
#     data = rect_df_ccf_bin,
#     aes(
#       x = ccf_bin,
#       y = (ymin + ymax) / 2,
#       label = ccf_bin
#     ),
#     inherit.aes = FALSE,
#     size = 3
#   )+
#   stat_summary(
#     fun = "median",
#     fun.min = function(x) boxplot.stats(x)$stats[2],
#     fun.max = function(x) boxplot.stats(x)$stats[4],
#     linewidth = 1,
#     size = 0.5,
#     position = position_dodge(width = 0.5)
#   ) +
#   geom_line(
#     data = median_df_purity,
#     aes(
#       x = CCF_bin,
#       y = sensitivity,
#       group = purity,
#       colour = purity
#     ),
#     linewidth = 0.3,
#     position = position_dodge(width = 0.5)
#   ) +
#   scale_color_manual("Somatic caller", values = purity_colors) +
#   scale_fill_manual(
#     name = "CCF class",
#     values = c(
#       "Neutral mutations"     = "grey40",
#       "Subclonal mutations" = "#8c00ec",
#       "Clonal mutations"         = "#ffae00"
#     )
#   ) +
#   ylab("Sensitivity") +
#   xlab("CCF bin") +
#   facet_wrap(~muts, ncol = 1,strip.position = "right") +
#   my_ggplot_theme()+
#   theme(
#     axis.text.x = element_blank(),
#     legend.position = "bottom"
#   )
# 
# purity_line_boxes
# 
# #### plot coverage 
# coverage_line_boxes <- df_all_SPN_somatic %>% 
#   mutate(
#     CCF_bin = factor(
#       CCF_bin,
#       levels = c("0-0.05","0.05-0.10","0.10-0.25",
#                  "0.25-0.50","0.50-0.75","0.75-0.95","0.95-1")
#     )
#   ) %>% 
#   mutate(coverage=as.factor(coverage)) %>% 
#   ggplot(aes(x = CCF_bin, y = sensitivity, col = as.factor(coverage))) +
#   geom_rect(
#     data = rect_df_ccf_bin,
#     aes(
#       xmin = stage(xmin, after_scale = xmin-0.48),
#       xmax = stage(xmax, after_scale = xmax+0.48),
#       ymin = ymin, ymax = ymax
#     ),
#     inherit.aes = FALSE,
#     alpha = 1,
#     fill="white"
#   ) +
#   geom_rect(
#     data = rect_df_ccf_bin,
#     aes(
#       xmin = stage(xmin, after_scale = xmin-0.48),
#       xmax = stage(xmax, after_scale = xmax+0.48),
#       ymin = ymin, ymax = ymax,
#       fill = region
#     ),
#     inherit.aes = FALSE,
#     alpha = 0.2
#   ) +
#   geom_text(
#     data = rect_df_ccf_bin,
#     aes(
#       x = ccf_bin,
#       y = (ymin + ymax) / 2,
#       label = ccf_bin
#     ),
#     inherit.aes = FALSE,
#     size = 3
#   )+
#   stat_summary(
#     fun = "median",
#     fun.min = function(x) boxplot.stats(x)$stats[2],
#     fun.max = function(x) boxplot.stats(x)$stats[4],
#     linewidth = 1,
#     size = 0.5,
#     position = position_dodge(width = 0.5)
#   ) +
#   geom_line(
#     data = median_df_coverage,
#     aes(
#       x = CCF_bin,
#       y = sensitivity,
#       group = coverage,
#       colour = coverage
#     ),
#     linewidth = 0.3,
#     position = position_dodge(width = 0.5)
#   ) +
#   scale_color_manual("Somatic caller", values = coverage_colors) +
#   scale_fill_manual(
#     name = "CCF class",
#     values = c(
#       "Neutral mutations"     = "grey40",
#       "Subclonal mutations" = "#8c00ec",
#       "Clonal mutations"         = "#ffae00"
#     )
#   ) +
#   ylab("Sensitivity") +
#   xlab("CCF bin") +
#   facet_wrap(~muts, ncol = 1,strip.position = "right") +
#   my_ggplot_theme()+
#   theme(
#     axis.text.x = element_blank(),
#     legend.position = "bottom"
#   )
# 
# coverage_line_boxes
# 
# # somatic_panel <- caller_line_boxes + purity_line_boxes + coverage_line_boxes + plot_layout(nrow = 3,heights = 1, guides = 'collect') &
# #   theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) 
# 
# ggsave(plot = caller_line_boxes, filename = '/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/somatic/caller_line_boxes.pdf',
#        width = 8, height = 5, units = 'in')
# #ggsave(plot = line, filename = '/orfeo/cephfs/scratch/area/lvaleriani/tesi/sarek/all_sarek.pdf',width = 6, height = 8, units = 'in')
# 
# 
# # tmp <- df_all_SPN_somatic %>% 
# #   dplyr::mutate(CCF_bin = case_when(
# #     CCF_bin == "Clonal" ~ "Clonal",
# #     CCF_bin %in% c("10-25%", "25-50%", "50-99%") ~ "10-99%",
# #     TRUE ~ "0-10%"
# #   )) %>% 
# #   select(sample_id, CCF_bin, pur, spn, cov, muts, sensitivity, caller) %>% 
# #   group_by(sample_id, CCF_bin, pur, spn, cov, muts, caller) %>%
# #   summarise(sensitivity = mean(sensitivity))
# #   
# # caller <- tmp %>% 
# #   dplyr::group_by(sample_id, CCF_bin, pur, spn, cov, muts) %>% 
# #   dplyr::select(sensitivity, caller) %>% 
# #   dplyr::mutate(normalized_sens = sensitivity / max(sensitivity)) %>% 
# #   ggplot(mapping = aes(x = CCF_bin, y = normalized_sens, fill = caller)) +
# #   geom_boxplot(outliers = F) +
# #   ylab('Normalized sensitivity') + 
# #   scale_fill_manual('Caller', values = method_colors) + 
# #   facet_wrap(~muts) +
# #   theme_minimal() +
# #   theme(axis.text.x = element_text(angle = 30, vjust = 1.2, hjust=1),axis.title.x = element_blank())
# # 
# # 
# # purity <- tmp %>% 
# #   dplyr::group_by(sample_id, CCF_bin, caller, spn, cov, muts) %>% 
# #   dplyr::select(sensitivity, pur) %>% 
# #   dplyr::mutate(normalized_sens = sensitivity / max(sensitivity)) %>% 
# #   ggplot(mapping = aes(x = CCF_bin, y = normalized_sens, fill = as.factor(pur))) +
# #   geom_boxplot(outliers = F) +
# #   scale_fill_manual('Purity', values = purity_colors) + 
# #   ylab('Normalized sensitivity') + 
# #   facet_wrap(~muts) +
# #   theme_minimal()+
# #   theme(axis.text.x = element_text(angle = 30, vjust = 1.2, hjust=1),axis.title.x = element_blank())
# # 
# # coverage <- tmp %>% 
# #   dplyr::group_by(sample_id, CCF_bin, caller, spn, pur, muts) %>% 
# #   dplyr::select(sensitivity, cov) %>% 
# #   dplyr::mutate(normalized_sens = sensitivity / max(sensitivity)) %>% 
# #   ggplot(mapping = aes(x = CCF_bin, y = normalized_sens, fill = as.factor(cov))) +
# #   geom_boxplot(outliers = F) +
# #   scale_fill_manual('Coverage', values = coverage_colors) + 
# #   ylab('Normalized sensitivity') + 
# #   facet_wrap(~muts) +
# #   theme_minimal()+
# #   theme(axis.text.x = element_text(angle = 30, vjust = 1.2, hjust=1),axis.title.x="CCF bin")
# # 
# # boxplot <- caller + purity + coverage + plot_layout(nrow = 3)
# # 
# # tmp_fp <- df %>% 
# #   dplyr::mutate(CCF_bin = case_when(
# #     CCF_bin == "Clonal" ~ "Clonal",
# #     CCF_bin %in% c("10-25%", "25-50%", "50-99%") ~ "10-99%",
# #     TRUE ~ "0-10%"
# #   )) %>% 
# #   select(sample_id, CCF_bin, pur, spn, cov, muts, false_positive, true_positives, caller) %>% 
# #   group_by(sample_id, CCF_bin, pur, spn, cov, muts, caller) %>%
# #   dplyr::summarise(FP = sum(false_positive, na.rm = TRUE), TP = sum(true_positives, na.rm = TRUE)) %>% 
# #   dplyr::mutate(FDR = FP / (TP + FP))
# # 
# # fp <- tmp_fp %>% 
# #   dplyr::group_by(caller, muts, cov, pur, spn) %>% 
# #   dplyr::select(FDR, caller) %>% 
# #   dplyr::mutate(normalized_fdr = FDR / max(FDR)) %>% 
# #   ggplot(mapping = aes(x = CCF_bin, y = normalized_fdr, fill = caller)) +
# #   geom_boxplot() + 
# #   facet_wrap(~muts) +
# #   theme_minimal() +
# #   scale_fill_manual(values = method_colors) +
# #   labs(x = "Caller", y = "False Discovery Rate", fill = "Caller")
# 
# 
# 
