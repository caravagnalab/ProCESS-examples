source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
somatic_performance <- function(df_all_SPN_somatic,type,median_df){
  if (type=="caller"){
    palette <- col_somatic_tools
    legend_title <- "Somatic caller"
  }
  if (type=="purity"){
    palette <- purity_colors
    legend_title <- "Purity"
  } 
  if (type=="coverage"){
    palette <- coverage_colors
    legend_title <- "Coverage"
  }
  rect_df_ccf <- tibble::tibble(
    xmin = c("0-0.05","0.10-0.25","0.95-1"),
    xmax = c("0.05-0.10","0.75-0.95","0.95-1"),
    ymin = -0.15,
    ymax = -0.05,
    region = c(
      "Subclonal low frequency",
      "Subclonal high frequency",
      "Clonal")) %>%
    mutate(xmid = xmin)
  
  rect_df_ccf_bin <- tibble::tibble(
    xmin = c("0-0.05","0.05-0.10","0.10-0.25","0.25-0.50","0.50-0.75","0.75-0.95","0.95-1"),
    xmax = c("0-0.05","0.05-0.10","0.10-0.25","0.25-0.50","0.50-0.75","0.75-0.95","0.95-1"),
    ymin = -0.15,
    ymax = -0.05,
    region = c("Subclonal low frequency","Subclonal low frequency",
               "Subclonal high frequency","Subclonal high frequency",
               "Subclonal high frequency","Subclonal high frequency",
               "Clonal"),
    ccf_bin=c("0-0.05","0.05-0.10","0.10-0.25","0.25-0.50","0.50-0.75","0.75-0.95","0.95-1"))
  plt <- df_all_SPN_somatic %>% 
    mutate(
      CCF_bin = factor(
        CCF_bin,
        levels = c("0-0.05","0.05-0.10","0.10-0.25",
                   "0.25-0.50","0.50-0.75","0.75-0.95","0.95-1")
      )
    ) %>% 
    ggplot(aes(x = CCF_bin, y = sensitivity, col = as.factor(.data[[type]]))) +
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
      alpha = 0.3
    ) +
    geom_text(
      data = rect_df_ccf_bin,
      aes(
        x = ccf_bin,
        y = (ymin + ymax) / 2,
        label = ccf_bin
      ),
      inherit.aes = FALSE,
      size = 2
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
      data = median_df,
      aes(
        x = CCF_bin,
        y = sensitivity,
        group = as.factor(.data[[type]]),
        colour = as.factor(.data[[type]])
      ),
      linewidth = 0.3,
      position = position_dodge(width = 0.5)
    ) +
    scale_color_manual(legend_title, values = palette) +
    scale_fill_manual(
      name = "Mutation class",
      values = c(
        "Subclonal low frequency"     = "grey40",
        "Subclonal high frequency" = "mediumpurple",
        "Clonal"         = "#ffae00"    ),
      breaks = c(
        "Subclonal low frequency",
        "Subclonal high frequency",
        "Clonal")
    ) +
    ylab("Sensitivity") +
    xlab("CCF bin") +
    facet_wrap(~muts, ncol = 1,strip.position = "right") +
    my_ggplot_theme()+
    theme(
      axis.text.x = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.minor = element_blank(),
      legend.spacing.y = unit(0.05, "cm"),
      legend.spacing.x = unit(0.1, "cm")
    )+
    guides(
      col  = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3)
    )
  return(plt)
  
}