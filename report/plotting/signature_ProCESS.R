library(ProCESS)
library(ggalluvial)
library(dplyr)
library(ggplot2)
library(patchwork)
source("plotting/utils.R")

get_colors_for <- function(values, pal_name = "Set3") {
  colors <- NULL
  if (length(values) > 0) {
    num_of_values <- values %>% length()
    if (num_of_values < 3) {
      colors <- RColorBrewer::brewer.pal(3, pal_name)
      colors <- colors[seq_len(num_of_values)]
    } else {
      colors <- RColorBrewer::brewer.pal(num_of_values, pal_name)
    }
    names(colors) <- values
  }
  return(colors)
}

plot_exposure_evolution <- function(sample_forest,phylo_forest,snapshot){
  simulation <- ProCESS::recover_simulation(snapshot)
  color_map_clones <- get_clone_map(sample_forest)
  muller <- plot_muller(simulation,color_map = color_map_clones)
  
  signature <- phylo_forest$get_exposures()
  palette_signature <-  get_colors_for(signature %>% pull(signature) %>% unique())
  
  df_sign <- get_exposure_ends(phylo_forest) %>% 
    mutate(time = round(time),
           end_time = round(end_time)-1)
  start_time <- 0
  end_time <- max(df_sign$end_time)
  
  df_final <- tibble()
  for (i in seq(start_time, end_time, 2)){
    for (irow in 1:nrow(df_sign)){
      row <- df_sign[irow,]
      if (i %in% seq(row$time, row$end_time)){
        tmp <- tibble(time = i, signature = row$signature, exposure = row$exposure, type = row$type)
        df_final <- bind_rows(df_final, tmp)
      }
    }
  }

  sign_1 <- ggplot(df_final, aes(x = factor(time), stratum = signature, alluvium = signature, y = exposure)) +
    geom_alluvium(aes(fill = signature), width = 0.001, curve_type = 'linear', decreasing = NA, knot.pos = 0, alpha = 1) +
    scale_x_discrete(breaks = seq(0, 300, 20)) +
    scale_fill_manual('Signature', values = palette_signature)+
    labs(legend = 'Signature',
         x = "",
         y = "Exposure") +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.margin = margin(0, 0, 0, 0, "pt")) + facet_grid(type~.)
  
  
  sign_muller <- sign_1 + muller + theme(legend.position = 'right') + plot_layout(nrow = 2)
  
  sign_tree <- plot_forest(sample_forest,color_map = color_map_clones) %>% 
    annotate_forest(forest = phylo_forest, samples = T, MRCAs = T, drivers = F, add_driver_label = F) + 
    theme(legend.position = 'left') + 
    ggplot2::guides(size = "none",
                    shape = "none",
                    color = ggplot2::guide_legend("Species")) +
    df_final %>%
    group_by(time, type) %>%
    reframe(exposure = exposure/sum(exposure), across(everything())) %>%
    ggplot()+
    geom_bar(aes(x = time, y = exposure, fill = signature), stat = "identity", width = 10) +
    ylab('Exposure') +
    xlab('') +
    scale_x_reverse() +
    scale_fill_manual('Signature', values = palette_signature)+
    theme_minimal() +
    coord_flip() + theme(axis.title.y=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks.y=element_blank(),
                         plot.margin = margin(0, 0, 0, 0, "pt")) + facet_grid(.~type)  + 
    plot_layout(nrow = 1, guides = 'collect')
  return(list("sign_tree"=sign_tree,"sign_muller"=sign_muller))
}
