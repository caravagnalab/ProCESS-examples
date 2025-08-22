library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)

purities <- c(0.3,0.6,0.9)
coverages <- c(50,100)
contexts <- c("SBS96","ID83")
all_combs <- list()
all_metrics <- list()
all_cosine  <- list()
spns <- c("SPN01","SPN03")

plot_all_spns <- lapply(spns,function(spn){
  setwd(paste0("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION/",spn,"/signature/"))
  all_combs <- list()
  all_metrics <- list()
  all_cosine  <- list()
  for (purity in purities){
    for (coverage in coverages){
      comb <- paste0(coverage,"x_",purity)
      metrics_spn <- list()
      cosine_spn  <- list()
      
      for (ctx in contexts){
        rds_metrics <- file.path(comb,paste0("metrics_",ctx,"_spn.rds"))
        rds_cosine  <- file.path(comb,paste0("cosine_mse_",ctx,".rds"))
        
        metrics_spn[[ctx]] <- readRDS(rds_metrics) %>% 
          dplyr::mutate(context=ctx,
                        purity=purity,
                        coverage=coverage)
        
        cosine_spn[[ctx]] <- readRDS(rds_cosine) %>% 
          dplyr::mutate(context=ctx,
                        purity=purity,
                        coverage=coverage)
      }
      
      # bind SBS96 + ID83 for this purity/coverage
      metrics_df <- dplyr::bind_rows(metrics_spn)
      cosine_df  <- dplyr::bind_rows(cosine_spn)
      
      # save into global list
      all_metrics[[comb]] <- metrics_df
      all_cosine[[comb]]  <- cosine_df
    }
  }
  
  # final dataframes
  final_metrics <- dplyr::bind_rows(all_metrics)
  final_cosine  <- dplyr::bind_rows(all_cosine)
  color_caller <- list("ID83"=c('BASCULE' = '#E1BFF8','SigProfiler' = '#F6AF92'),
                       "SBS96"=c('BASCULE' = 'purple1','SigProfiler' = 'orange2'))
  
  plots <- list()
  for (ctx in contexts){
    
    cosine_plot <-final_cosine %>%
      filter(context==ctx) %>% 
      ggplot(aes(x = sample, y = cosine, fill = caller)) +
      geom_col(position = position_dodge()) +
      geom_errorbar(
        aes(ymin = cosine - mse, ymax = cosine + mse),
        position = position_dodge(width = 0.9),
        width = 0.3
      ) +
      scale_fill_manual(values = color_caller[[ctx]]) +
      facet_grid(as.numeric(coverage) ~ as.numeric(purity)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      ylab(label = paste0("cosine ",ctx))
    
    
    metrics_plot <- final_metrics %>%
      dplyr::filter(context==ctx) %>% 
      ggplot(aes(x = metric, y = mean, fill = caller)) +
      geom_col(position = position_dodge()) +
      geom_errorbar(
        aes(x=metric,ymin = mean - sd, ymax = mean + sd),
        position = position_dodge(width = 0.9),
        width = 0.3
      ) +
      scale_fill_manual(values = color_caller[[ctx]]) +
      facet_grid(as.numeric(coverage) ~ as.numeric(purity)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      ylab(label = paste0("mean ",ctx))
    
    plots[[ctx]] <- wrap_plots(list(cosine_plot,metrics_plot))+plot_layout(guides = "collect")
  }

   wrap_plots(plots,nrow = 2)+plot_annotation(title = spn)

})

