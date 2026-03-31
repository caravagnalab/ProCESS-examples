library(ProCESS)
library(ggplot2)
library(dplyr)
library(ggalluvial)
library(tidyr)
library(tibble)
library(purrr)
library(ggpubr)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
library(ggnewscale)
library(tidyr)

source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
# source("compute_FGA.R")
scout_dir <-"/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
SPNS <- c("SPN01","SPN02","SPN03","SPN04", "SPN05","SPN06", 'SPN07')
COVERAGES <- c("50","100", "150")
PURITIES <- c("0.3","0.6","0.9")

source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
WGD_samples <- c("SPN01_1.1","SPN01_1.3","SPN06_3.2")

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
fga_fgs_df <- readRDS(file = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/fga_df.rds") %>% 
  dplyr::rename(sample=spn)
df_all_combs_SPN_cna <- df_all_combs_SPN_cna %>% 
  dplyr::select(!c("fga","fgs")) %>% 
  left_join(fga_fgs_df) %>% 
  mutate(tool = case_when(
    tool == 'ascat' ~ 'ASCAT',
    tool == 'sequenza' ~ 'Sequenza',
    tool == 'battenberg' ~ 'Battenberg',
  ))  %>% 
  mutate(abs_delta_purity=abs(delta_purity)) %>% 
  mutate(abs_delta_ploidy=abs(delta_ploidy)) %>% 
  mutate(tool = factor(tool, levels = c("ASCAT","Sequenza","Battenberg"))) %>% 
  mutate(class=case_when(abs(delta_purity) > 0.20 & abs(delta_ploidy) < 0.7 ~ 'uncorrect purity',
                         abs(delta_purity) < 0.20 & abs(delta_ploidy) > 0.7 ~ 'uncorrect ploidy',
                         abs(delta_purity) <= 0.20 & abs(delta_ploidy) <= 0.7 ~ 'correct',
                         TRUE ~ "uncorrect"))
  
summary_stats <-df_all_combs_SPN_cna %>% 
  # mutate(class = ifelse(abs(delta_purity) >= 0.20 | abs(delta_ploidy) > 0.7, 'not correct', 'correct')) %>% 
  group_by(tool, class) %>% 
  summarise(n = n(), .groups = "drop_last") %>% 
  mutate(percentage = round(100 * n / sum(n),0)) %>% 
  mutate(label_percentage=paste0(as.character(percentage),"%")) %>% 
  ungroup()

df_all_combs_SPN_cna_all_metrics <- df_all_combs_SPN_cna %>% 
  left_join(summary_stats,relationship = "many-to-many")

# rect_df <- tibble::tibble(
#   xmin = c(0,   0.2, 0,   0.2),
#   xmax = c(0.2, Inf, 0.2, Inf),
#   ymin = c(0,   0,   0.7, 0.7),
#   ymax = c(0.7, 0.7, Inf, Inf),
#   class = c(
#     "correct",
#     "not correct",
#     "not correct",
#     "not correct"
#   ),
#   subregion = c("Correct",
#                 "Correct Ploidy",
#                 "Correct Purity",
#                 "Incorrect Both")
# )

rect_df <- tibble::tibble(
  xmin = c(0,   0.2, 0,   0.2),
  xmax = c(0.2, max(df_all_combs_SPN_cna_all_metrics$abs_delta_purity), 0.2, max(df_all_combs_SPN_cna_all_metrics$abs_delta_purity)),
  ymin = c(0,   0,   0.7, 0.7),
  ymax = c(0.7, 0.7, max(df_all_combs_SPN_cna_all_metrics$abs_delta_ploidy), max(df_all_combs_SPN_cna_all_metrics$abs_delta_ploidy)),
  class = c(
    "correct",
    "uncorrect purity",
    "uncorrect ploidy",
    "not correct"
  ),
  subregion = c("Correct",
                "Incorrect Purity",
                "Incorrect Ploidy",
                "Incorrect Both")
)

# rect_df <- rect_df %>% 
#   mutate(xmid=case_when(subregion=="Correct"~((xmax-xmin)/2),
#                         # subregion=="Incorrect Both" ~ (max(df_all_combs_SPN_cna$abs_delta_purity)/2),
#                         TRUE ~ NA)) %>% 
#   mutate(ymid=case_when(subregion=="Correct"~((ymax-ymin)/2),
#                         # subregion=="Incorrect Both" ~ (max(df_all_combs_SPN_cna$abs_delta_ploidy)/2),
#                         TRUE ~ NA))

rect_df <- rect_df %>% 
    mutate(xmid=case_when(subregion=="Correct"~(xmax-xmin)/2,
                          subregion=="Incorrect Purity"~0.5,
                          subregion=="Incorrect Ploidy"~0.1,
                          # subregion=="Incorrect Both" ~ (max(df_all_combs_SPN_cna$abs_delta_purity)/2),
                          TRUE ~ NA)) %>%
  mutate(ymid=case_when(subregion=="Correct"~(ymax-ymin)/2,
                        subregion=="Incorrect Purity"~0.3,
                        subregion=="Incorrect Ploidy"~1.5,
                        # subregion=="Incorrect Both" ~ (max(df_all_combs_SPN_cna$abs_delta_purity)/2),
                        TRUE ~ NA))
                          


df_all_combs_SPN_cna_all_metrics <- df_all_combs_SPN_cna_all_metrics %>% 
  left_join(rect_df,by = "class",relationship = "many-to-many")


##### Plot breakpoint
plt_breakpoint <- df_all_combs_SPN_cna %>% 
  group_by(tool, fga_class, coverage, true_purity) %>%
  summarise(
    mean_precision = mean(precision, na.rm = TRUE),
    sd_precision   = sd(precision, na.rm = TRUE),
    mean_recall    = mean(recall, na.rm = TRUE),
    sd_recall      = sd(recall, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # mutate(group_ellipse = interaction(fga_class, tool, drop = TRUE)) %>% 
  mutate(tool = factor(tool, levels = c("ASCAT","Sequenza","Battenberg"))) %>% 
  ggplot(aes(
    x = mean_recall,
    y = mean_precision
  )) +
  geom_smooth(colour = "grey60",fill="grey70",method = "lm")+
  geom_point(
    aes(
      fill = fga_class,
      color=fga_class
    ),size=2
  ) +
  
  scale_fill_manual(values = c("High FGA" = "indianred2", "Low FGA"="dodgerblue3"), name = "FGA class") +
  scale_color_manual(values = c("High FGA" = "indianred2", "Low FGA"="dodgerblue3"), name = "FGA class") +
  my_ggplot_theme()+
  labs(
    x = "Mean breakpoint recall",
    y = "Mean breakpoint precision"
  )+
  facet_wrap(~tool) +
  theme(legend.position = 'bottom',
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.box.background = element_blank(), 
        legend.spacing.x = unit(0.1, "cm"))

plt_breakpoint
#ggsave(filename = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/cna/panel_D.pdf",plot = plt_breakpoint,width = 10,height = 4)

##### Plot purity/ploidy 

scatter_purity_ploidy <- df_all_combs_SPN_cna_all_metrics %>%
  mutate(tool = factor(tool, levels = c("ASCAT","Sequenza","Battenberg"))) %>% 
  ggplot(aes(x = abs(delta_purity), y = abs(delta_ploidy))) +
  
  geom_rect(
    data = rect_df,
    aes(
      xmin = xmin, xmax = xmax,
      ymin = ymin, ymax = ymax,
      fill = class
    ),
    inherit.aes = FALSE,
    alpha=0.2
  ) +
  
  geom_point(aes(col = fga_class)) +
  # geom_text(aes(
  #     x = xmid,
  #     y = ymid,
  #     label = label_percentage
  #   ),
  #   inherit.aes = FALSE,
  #   size = 3,
  #   color="gray30"
  # )+
  facet_wrap(~tool) +
  
  geom_vline(xintercept = 0.2, col = "gray60", linetype = "dashed") + 
  geom_hline(yintercept = 0.7, col = "gray60", linetype = "dashed") +
  
  scale_fill_manual(
    name = "Caller estimate",
    values = c(
      "uncorrect ploidy"   = 'coral4',
      "uncorrect purity" = "#EEAD0E",
      "correct" = "#548B54",
      "not correct" = "transparent"
      # "not correct" = "#EEAD0E"
    )
  ) +
  
  scale_color_manual(
    name = "FGA class",
    values = c(
      "High FGA" = "indianred2",
      "Low FGA"  = "dodgerblue3"
    )
  ) +
  
  my_ggplot_theme() +
  scale_y_continuous(breaks = c(0,0.7,1.5,2))+
  xlab("Absoulte error (purity)") +
  ylab("Absoulte error (ploidy)") +
  theme(legend.position = 'bottom',
        legend.box = "horizontal",
        legend.spacing.x = unit(0.1, "cm"),
        panel.grid.minor = element_blank(),
        legend.spacing.y = unit(0.00001, "cm"))+
  guides(
    col  = guide_legend(nrow = 2),
    fill = guide_legend(nrow = 2, override.aes = list(alpha = 1))
  )
scatter_purity_ploidy
#ggsave(filename = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/cna/panel_B.pdf",plot = scatter_purity_ploidy,width = 10,height = 4)

##### Alluvial plot
df <- df_all_combs_SPN_cna %>% 
  # mutate(class = ifelse(abs(delta_purity) >= 0.20 | abs(delta_ploidy) > 0.7,
  #                       'not correct', 'correct')) %>% 
  mutate(class=case_when(abs(delta_purity) > 0.20 & abs(delta_ploidy) < 0.7 ~ 'uncorrect purity',
                         abs(delta_purity) < 0.20 & abs(delta_ploidy) > 0.7 ~ 'uncorrect ploidy',
                         abs(delta_purity) <= 0.20 & abs(delta_ploidy) <= 0.7 ~ 'correct',
                         TRUE ~ "uncorrect")) %>% 
  mutate(id=paste0(sample,"_",coverage,"_",true_purity)) %>% 
  mutate(comb=paste0(coverage,"_",true_purity)) %>% 
  dplyr::select(tool,class,id,fga_class) %>% 
  mutate(tool = factor(tool, levels = c("ASCAT","Sequenza","Battenberg")))

alluvial_plot <- ggplot(
  df,
  aes(
    x = tool,
    stratum = class,
    alluvium = id,
    fill = class
  )
) +
  stat_stratum(
    decreasing = TRUE,
    alpha = .4,
    width = 0.25,
    color = "grey40"
  ) +
  stat_flow(
    decreasing = TRUE,
    alpha = .2,
    width = 0.25,
    color = "grey40"
  ) +
  
  ## 🔑 percentage labels per stratum (per tool)
  stat_stratum(
    geom = "text",
    aes(
      label = scales::percent(after_stat(prop), accuracy = 2)
    ),
    size = 3.2,
    color = "black",
    decreasing = TRUE
  ) +
  
  scale_fill_manual(
    name = "Caller estimate",
    values = c(
      "uncorrect ploidy" = "coral4",
      "uncorrect purity" = "#EEAD0E",
      "correct" = "#548B54",
      "not correct" = "transparent"
    )
  ) +
  my_ggplot_theme() +
  theme(
    legend.position = "bottom",
    legend.box.spacing = unit(0, "cm"),
    panel.grid = element_blank(),
    legend.spacing.x = unit(0.1, "cm"),
    panel.background = element_blank()
  ) +
  xlab("") +
  ylab("Frequency") +
  scale_x_discrete(expand = c(0.1, 0.1))

alluvial_plot

#ggsave(filename = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/cna/panel_C.pdf",plot = alluvial_plot,width = 10,height = 4)


#plt_breakpoint/scatter_purity_ploidy/alluvial_plot
#write.table(x = df_all_combs_SPN_cna_all_metrics,file = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/cna/cna_metrics.tsv",append = F,quote = F,sep = "\t",row.names = F)

# ggsave(plot = scatter_purity_ploidy, filename = '/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/cna/scatter_ploidy_purity.pdf',
#        width = 6, height = 3, units = 'cm')
# ggsave(plot = alluvial_plot, filename = '/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/cna/alluvial_plot.pdf',
#        width = 15, height = 5, units = 'cm')

############# OLD STAFF ################
# df <- df_all_combs_SPN_cna %>% 
#   mutate(class = ifelse(abs(delta_purity) >= 0.20 | abs(delta_ploidy) > 0.7,
#                         'not correct', 'correct')) %>% 
#   select(tool,class,sample) %>% 
#   ggplot(aes(x = tool, stratum = class, alluvium = sample,
#              fill = class, label = class)) +
#   scale_fill_brewer(type = "qual", palette = "Set2") +
#   geom_flow() +
#   geom_stratum() +
#   theme(legend.position = "bottom") 
# 
# 
# pie_chart <- df_all_combs_SPN_cna %>% 
#   mutate(tool = factor(tool, levels = c("ASCAT","Sequenza","Battenberg"))) %>% 
#   mutate(class = ifelse(abs(delta_purity) >= 0.20 | abs(delta_ploidy) > 0.7, 'not correct', 'correct')) %>% 
#   group_by(tool, class) %>% 
#   summarize(n = n()) %>% 
#   ggplot() + 
#   geom_bar(aes(x = '', y = n, fill = class), stat ='identity',alpha=0.8) +
#   scale_fill_manual('Prediction', values = c('correct' = "#548B54", 'not correct' = '#EEAD0E')) +
#   coord_polar(theta = "y", start=0) +
#   facet_wrap(.~tool) + 
#   xlab('') +
#   ylab('') +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         legend.key.size=unit(0.3, "cm"),
#         panel.background=element_rect(fill="white"),
#         axis.title=element_text(size=10),
#         legend.text=element_text(size=8),
#         legend.title=element_text(size=10),
#         strip.text = element_text(color = "white"),
#         text=element_text(size=10),
#         plot.title=element_text(size=12))
#   
# 
# 
# 
# 
# df_all_combs_SPN_cna_long <- df_all_combs_SPN_cna %>%
#   mutate(estimation_class = case_when(
#     abs_delta_purity < 0.15 & abs_delta_ploidy < 0.7 ~ "correctly estimated",
#     abs_delta_purity > 0.15 & abs_delta_ploidy < 0.7 ~ "uncorrect purity",
#     abs_delta_purity < 0.15 & abs_delta_ploidy > 0.7 ~ "uncorrect ploidy",
#     TRUE ~ "all uncorrect")
#     ) %>% 
#   #   correct_ploidy_estimation_class = case_when(
#   #     abs(delta_ploidy) < 0.7 ~ "correctly estimated",
#   #     delta_ploidy > 0.7 ~ "underestimated",
#   #     delta_ploidy < -0.7 ~ "overestimated"),
#   #   true_purity_label = paste0('Purity =', true_purity)
#   # ) %>%
#   pivot_longer(
#     cols = c(delta_purity, delta_ploidy),
#     names_to = "delta_type",
#     values_to = "delta_value"
#   ) %>%
#   mutate(delta_type=case_when(delta_type=="delta_ploidy" ~ "Ploidy",
#                               delta_type=="delta_purity" ~ "Purity"))
# 
# df_all_combs_SPN_cna_long <- df_all_combs_SPN_cna %>%
#   mutate(correct_purity_estimation_class = case_when(
#     abs(delta_purity) < 0.15 ~ "correctly estimated",
#     delta_purity > 0.15 ~ "underestimated",
#     delta_purity < -0.15 ~ "overestimated"),
#     correct_ploidy_estimation_class = case_when(
#       abs(delta_ploidy) < 0.7 ~ "correctly estimated",
#       delta_ploidy > 0.7 ~ "underestimated",
#       delta_ploidy < -0.7 ~ "overestimated"),
#     true_purity_label = paste0('Purity =', true_purity)
#   ) %>%
#   pivot_longer(
#     cols = c(delta_purity, delta_ploidy),
#     names_to = "delta_type",
#     values_to = "delta_value"
#   ) %>%
#   mutate(
#     correct_estimation_class = ifelse(delta_type == "delta_purity",
#                                       correct_purity_estimation_class,
#                                       correct_ploidy_estimation_class)
#   ) %>% 
#   mutate(delta_type=case_when(delta_type=="delta_ploidy" ~ "Ploidy",
#                               delta_type=="delta_purity" ~ "Purity"))
# 
# 
# 
# df_plot <- df_all_combs_SPN_cna_long %>%
#   count(tool, estimation_class) %>%
#   group_by(tool) %>%
#   mutate(prop = n / sum(n))
# 
# ggplot(df_plot, aes(x = "", y = prop, fill = estimation_class)) +
#   geom_col(width = 1, color = "white") +
#   coord_polar(theta = "y") +
#   facet_wrap(~ tool) +
#   labs(y = NULL, x = NULL, fill = "Estimation class") +
#   scale_y_continuous(labels = scales::percent) +
#   theme_void() +
#   theme(strip.text = element_text(size = 12))
# 
# ggplot(df_plot, aes(x = tool, y = prop, fill = estimation_class)) +
#   geom_bar(stat = "identity", position = "fill") +
#   coord_polar()+
#   # facet_wrap(~tool) +
#   labs(y = "Percentage", x = "Tool") +
#   scale_y_continuous(labels = scales::percent) +
#   scale_fill_manual(
#     values = c(
#       "overestimated" = "indianred2",
#       "underestimated" = "dodgerblue3",
#       "correctly estimated" = "springgreen4"
#     )
#   ) +
#   my_ggplot_theme() +
#   stat_pvalue_manual(
#     pval_df,
#     label = "p.adj.signif",
#     tip.length = 0.02,
#     bracket.shorten = 0.05
#   ) +
#   coord_cartesian(ylim = c(0, 1.1), clip = "off")
# 
# 
# pval_df <- df_plot %>%
#   group_by(delta_type) %>%
#   summarise(
#     p = list({
#       mat <- cur_data() %>%
#         select(fga_class, correct_estimation_class, n) %>%
#         pivot_wider(names_from = correct_estimation_class, values_from = n) %>%
#         column_to_rownames("fga_class") %>%
#         as.matrix()
#       fisher.test(mat)$p.value
#     })
#   ) %>%
#   mutate(
#     p = unlist(p),
#     group1 = "High FGA",
#     group2 = "Low FGA",
#     p.adj = p,
#     p.adj.signif = case_when(
#       p < 0.0001 ~ "****",
#       p < 0.001  ~ "***",
#       p < 0.01   ~ "**",
#       p < 0.05   ~ "*",
#       TRUE       ~ "ns"
#     ),
#     y.position = 1.05
#   )
# p_purity_ploidy <-ggplot(df_plot, aes(x = fga_class, y = prop, fill = correct_estimation_class)) +
#   geom_bar(stat = "identity", position = "fill") +
#   facet_wrap(~delta_type) +
#   labs(y = "Percentage", x = "FGA class") +
#   scale_y_continuous(labels = scales::percent) +
#   scale_fill_manual(
#     values = c(
#       "overestimated" = "indianred2",
#       "underestimated" = "dodgerblue3",
#       "correctly estimated" = "springgreen4"
#     )
#   ) +
#   my_ggplot_theme() +
#   stat_pvalue_manual(
#     pval_df,
#     label = "p.adj.signif",
#     tip.length = 0.02,
#     bracket.shorten = 0.05
#   ) +
#   coord_cartesian(ylim = c(0, 1.1), clip = "off")
# 
# saveRDS(object = plt_breakpoint,file = "plt_breakpoint.rds")
# plt_purity_ploidy_old <- ggplot(df_all_combs_SPN_cna_long, aes(
#   x = reorder(tool, delta_value, \(x) sd(x, na.rm = TRUE)),
#   # x = true_purity_label,
#   y = delta_value,
#   group = interaction(tool, true_purity_label)
# )) +
#   
#   geom_boxplot(alpha = 0.3, outlier.shape = NA) +
#   # scale_color_manual(name = "CNA caller", values = col_cna_tools) +
#   # scale_fill_manual(name = "CNA caller", values = col_cna_tools) +
#   
#   ggnewscale::new_scale_color() +
#   ggnewscale::new_scale_fill() +
#   
#   # Jitter points colored by FGA class
#   geom_jitter(aes(color = fga_class, fill = fga_class, shape = correct_estimation_class),
#               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
#               size = 2, alpha = 0.8) +
#   scale_color_manual(name = "FGA class", values = c("High FGA"="indianred2","Low FGA"="dodgerblue3")) +
#   scale_fill_manual(name = "FGA class", values = c("High FGA"="indianred2","Low FGA"="dodgerblue3")) +
#   scale_shape_manual(name = "Estimation class",
#                      values = c("correctly estimated"=16,
#                                 "underestimated"=25,
#                                 "overestimated"=24)) +
#   
#   # Horizontal lines
#   geom_hline(aes(yintercept = 0), linetype = "solid", color = "grey40") +
#   geom_hline(data = subset(df_all_combs_SPN_cna_long, delta_type=="Purity"),
#              aes(yintercept = 0.2), linetype = "dashed", color = "grey") +
#   geom_hline(data = subset(df_all_combs_SPN_cna_long, delta_type=="Purity"),
#              aes(yintercept = -0.2), linetype = "dashed", color = "grey") +
#   geom_hline(data = subset(df_all_combs_SPN_cna_long, delta_type=="Ploidy"),
#              aes(yintercept = 1), linetype = "dashed", color = "grey") +
#   geom_hline(data = subset(df_all_combs_SPN_cna_long, delta_type=="Ploidy"),
#              aes(yintercept = -1), linetype = "dashed", color = "grey") +
#   
#   # Facet grid: delta_type as rows, true purity as columns
#   facet_grid(delta_type ~ true_purity_label, scales = "free_y") +
#   labs(y = "Deviation", x = "") + my_ggplot_theme()+
#   theme(axis.text.x = element_text(angle = 30,hjust = 1),
#         axis.title.x = element_blank())
# 
# #saveRDS(object = plt_purity_ploidy_old,file = "plt_purity_ploidy_old.rds")
# 
# 
# 
# plt_purity_ploidy <- df_all_combs_SPN_cna_long %>% 
#   ggplot(aes(
#     x = factor(true_purity),
#     y = delta_value,
#     group = interaction(factor(true_purity), fga_class),
#   )) +
#   geom_boxplot(
#     aes(color = fga_class),
#     alpha = 0.3, 
#     outlier.shape = NA) +
#   geom_jitter(
#     aes(color = fga_class, fill = fga_class, shape = correct_estimation_class),
#     position = position_jitterdodge(jitter.width = 0.2), #jitter.width = 0.2,
#     size = 2, alpha = 0.5
#   ) +
#   scale_color_manual(name = "FGA class", values = c("High FGA"="indianred2","Low FGA"="dodgerblue3")) +
#   scale_fill_manual(name = "FGA class", values = c("High FGA"="indianred2","Low FGA"="dodgerblue3")) +
#   scale_shape_manual(name = "Estimation class",
#                      values = c("correctly estimated"=16,
#                                 "underestimated"=25,
#                                 "overestimated"=24)) +
#   geom_hline(aes(yintercept = 0), linetype = "solid", color = "grey40") +
#   geom_hline(data = subset(df_all_combs_SPN_cna_long, delta_type=="Purity"),
#              aes(yintercept =  0.2), linetype = "dashed", color = "grey") +
#   geom_hline(data = subset(df_all_combs_SPN_cna_long, delta_type=="Purity"),
#              aes(yintercept = -0.2), linetype = "dashed", color = "grey") +
#   geom_hline(data = subset(df_all_combs_SPN_cna_long, delta_type=="Ploidy"),
#              aes(yintercept =  1), linetype = "dashed", color = "grey") +
#   geom_hline(data = subset(df_all_combs_SPN_cna_long, delta_type=="Ploidy"),
#              aes(yintercept = -1), linetype = "dashed", color = "grey") +
#   facet_wrap(~ delta_type, scales = "free_y",
#              labeller = as_labeller(c(
#                Ploidy = "Estimated Ploidy",
#                Purity = "Estimated Purity"
#              ))) +
#   labs(y = "Deviation", x = "Purity") +
#   my_ggplot_theme()
#   
# # df_all_combs_SPN_cna_long %>% 
# #   ggplot(aes(
# #   # x = reorder(tool, delta_value, \(x) sd(x, na.rm = TRUE)),
# #   x = true_purity,
# #   y = delta_value
# #   # group = interaction(tool, true_purity_label)
# # )) +
# #   geom_boxplot(aes(color = fga_class),alpha = 0.3, outlier.shape = NA) +
# #   # scale_color_manual(name = "CNA caller", values = col_cna_tools) +
# #   # scale_fill_manual(name = "CNA caller", values = col_cna_tools) +
# #   
# #   # ggnewscale::new_scale_color() +
# #   # ggnewscale::new_scale_fill() +
# #   
# #   # Jitter points colored by FGA class
# #   geom_jitter(aes(color = fga_class, fill = fga_class, shape = correct_estimation_class),
# #               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
# #               size = 2, alpha = 0.5) +
# #   scale_color_manual(name = "FGA class", values = c("High FGA"="indianred2","Low FGA"="dodgerblue3")) +
# #   scale_fill_manual(name = "FGA class", values = c("High FGA"="indianred2","Low FGA"="dodgerblue3")) +
# #   scale_shape_manual(name = "Estimation class",
# #                      values = c("correctly estimated"=16,
# #                                 "underestimated"=25,
# #                                 "overestimated"=24)) +
# #   
# #   # Horizontal lines
# #   geom_hline(aes(yintercept = 0), linetype = "solid", color = "grey40") +
# #   geom_hline(data = subset(df_all_combs_SPN_cna_long, delta_type=="Purity"),
# #              aes(yintercept = 0.2), linetype = "dashed", color = "grey") +
# #   geom_hline(data = subset(df_all_combs_SPN_cna_long, delta_type=="Purity"),
# #              aes(yintercept = -0.2), linetype = "dashed", color = "grey") +
# #   geom_hline(data = subset(df_all_combs_SPN_cna_long, delta_type=="Ploidy"),
# #              aes(yintercept = 1), linetype = "dashed", color = "grey") +
# #   geom_hline(data = subset(df_all_combs_SPN_cna_long, delta_type=="Ploidy"),
# #              aes(yintercept = -1), linetype = "dashed", color = "grey") +
# #   
# #   # Facet grid: delta_type as rows, true purity as columns
# #   # facet_grid(delta_type ~ true_purity_label, scales = "free_y") +
# #   facet_wrap(~delta_type, scales = "free_y", labeller = as_labeller(
# #     c(Ploidy = "Estimated Ploidy", Purity = "Estimated Purity"))) +
# #   labs(y = "Deviation", x = "Purity") + 
# #   my_ggplot_theme() #+
# #   #theme(#axis.text.x = element_text(angle = 30,hjust = 1),
# #   #      axis.title.x = element_blank()) 
# 
# plt_purity_ploidy
# #saveRDS(object = plt_purity_ploidy,file = "plt_purity_ploidy.rds")
