library(ggplot2)
library(patchwork)
library(ProCESS)
library(dplyr)
library(tidyverse)

setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/report')
source("../getters/process_getters.R")
source("plotting/signature_ProCESS.R")
source("plotting/plot_genome_wide.R")
source("plotting/dynamics_ProCESS.R",local =T)
source("plotting/tables.R", local = knitr::knit_global())
metadata <- read.table(file = "SCOUT_metadata.csv",header = T,sep = "\t")


sample_forest <- load_sample_forest(get_sample_forest('SPN03'))
phylo_forest <- load_phylogenetic_forest(get_phylo_forest('SPN03'))
sample_names <- sapply(get_process_cna('SPN03'), function(path) {
  base_name <- basename(path)
  sub("_cna.rds$", "", base_name)
})
number_of_samples <- length(sample_names)
cna_data <- lapply(get_process_cna('SPN03'),readRDS)

info_spn <- metadata %>% filter(SPN_ID=='SPN03')

# sample dynamics
setwd('/orfeo/scratch/cdslab/shared/SCOUT/SPN03/process/')
simulation <- ProCESS::recover_simulation('SPN03')
simulation_info <- simulation$get_samples_info()
color_map_clones <- c('Clone 1' = '#4f5d75', 'Clone 2' = '#a0c1b9', 
                      'Clone 3' = '#f4d35e', 'Wild-type' = 'white', 
                      `NA` = 'white')

timing <- simulation_info %>% mutate(what = '')

muller <- plot_muller(simulation) + 
  geom_vline(xintercept = timing$time, linetype = 'dashed') + 
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = color_map_clones) 

muller

samples_name <- simulation_info$name
plot_sample <- lapply(samples_name, FUN = function(s){
  plot_tissue(simulation = simulation, before_sample = s,color_map = color_map_clones) +
    ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank())
})
plot_sampling <- patchwork::wrap_plots(plot_sample, nrow = 1, ncol = 4, guides = 'collect') & theme(legend.position = 'none')


forest <- plot_forest(sample_forest) %>% 
  annotate_forest(forest = sample_forest) +
  scale_color_manual(values = color_map_clones) +
  labs(y = "Generation") +
  theme_bw() +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none"
  )



plot_net_growth_rate <- function(df, last_t) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # Ensure proper types
  df <- df %>%
    mutate(
      time = as.numeric(time),
      rate = as.numeric(rate),
      event = tolower(event)
    )
  
  # Create complete grid and fill missing values (LOCF)
  full_df <- df %>%
    select(time, mutant, event, rate) %>%
    complete(mutant, event = c("growth", "death"), time = unique(df$time)) %>%
    arrange(mutant, event, time) %>%
    group_by(mutant, event) %>%
    tidyr::fill(rate, .direction = "down") %>%
    ungroup() %>%
    mutate(rate = ifelse(is.na(rate), 0, rate))
  
  # Pivot to wide format with net growth rate
  net_df <- full_df %>%
    pivot_wider(names_from = event, values_from = rate, values_fill = 0) %>%
    mutate(net_growth = growth - death) %>%
    arrange(mutant, time)
  
  # Add a row at last_t for each mutant if missing
  padded_df <- net_df %>%
    group_by(mutant) %>%
    summarize(last_time = max(time[time <= last_t]), .groups = "drop") %>%
    left_join(net_df, by = c("mutant", "last_time" = "time")) %>%
    filter(!is.na(net_growth)) %>%
    mutate(time = last_t) %>%
    select(names(net_df))  # ensure column order is same
  
  # Append to main data and filter only up to last_t
  net_df <- bind_rows(net_df, padded_df) %>%
    arrange(mutant, time) %>%
    filter(time <= last_t)
  
  # Plot
  ggplot(net_df, aes(x = time, y = net_growth, color = mutant)) +
    geom_step(size = 1) +
    labs(
      x = "Generation",
      y = "Growth Rate",
      color = "Mutant"
    ) +
    theme_minimal(base_size = 14)
}

rates_plot <- plot_net_growth_rate(df = simulation$get_rates_update_history(), last_t = simulation$get_clock()) +
  scale_color_manual(values = color_map_clones) +
  theme_bw() +
  theme(legend.position = "none")


p = plot_sampling / ((muller / rates_plot) | forest) +
  plot_layout(heights = c(.9, 2)) & 
  theme(text = element_text(size = 16))
#p

base <- '/orfeo/cephfs/scratch/area/lvaleriani/tesi/spn03/'
ggsave(file.path(base, "sim_2.pdf"), plot = p, width = 16, height = 12, units = "in", device = cairo_pdf)
