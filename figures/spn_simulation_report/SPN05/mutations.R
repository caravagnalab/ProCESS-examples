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


# sample_forest <- load_sample_forest(get_sample_forest('SPN05'))
# #phylo_forest <- load_phylogenetic_forest(get_phylo_forest('SPN05'))
# sample_names <- get_sample_names('SPN05')
# number_of_samples <- length(sample_names)
# 
# data <- readRDS(get_mutations(spn = 'SPN05', coverage = 50, purity = 0.9, type = 'tumour'))
# # normal <- readRDS(get_mutations(spn = 'SPN05', coverage = 30, purity = 1, type = 'normal')) %>%
# #   filter(classes == 'germinal') %>% 
# #   seq_to_long()
# 
# # germline
# # germline <- data %>% 
# #   filter(classes == 'germinal') 
# 
# # cna
# # cna <- phylo_forest$get_bulk_allelic_fragmentation() 
# 
# # n_cov = 30
# # t_cov = 50
# # ratio = n_cov/t_cov
# # g_sample <- germline %>%   
# #   ungroup() %>% 
# #   left_join(normal %>% ungroup(), join_by(chr,from,to), suffix = c('_T', '_N')) %>% 
# #   filter(VAF_N > 0.3) %>% 
# #   filter(VAF_N < 0.7) %>%
# #   filter(DP_T > 15) %>% 
# #   sample_n(size = 1e4) %>% 
# #   mutate(DR = (DP_T/DP_N)*ratio) %>% 
# #   left_join(cna) %>% 
# #   filter(from >= begin, to <= end) %>% 
# #   mutate(CN = paste(major, minor, sep =':'))
# # 
# # col <- c('1:1'='slategray4', '1:0'='goldenrod3')
# # germline <- g_sample %>% 
# #   ggplot() +
# #   geom_point(aes(x = from, y = VAF_T, col = CN), size = .3) +
# #   facet_grid(.~sample_name_T) +
# #   scale_color_manual(values = col)+
# #   ylab('BAF') +
# #   xlab('chr13') +
# #   ylim(0,1) + 
# #   guides(color = guide_legend(override.aes = list(size = 3))) +
# #   
# #   g_sample %>% 
# #   ggplot() +
# #   geom_point(aes(x = from, y = DR, , col = CN), size = .3) +
# #   facet_grid(.~sample_name_T) +
# #   scale_color_manual(values = col)+
# #   ylab('DR') +
# #   xlab('chr13') +
# #   ylim(0,3) +
# #   guides(color = guide_legend(override.aes = list(size = 3))) +
# #   plot_layout(nrow = 2, guides = 'collect') & theme_minimal()
# # 
# # germline
# 
# 
# somatic <- data %>%
#   dplyr::filter(classes!="germinal") %>%
#   dplyr::filter(!stringr::str_detect(causes, 'errors')) %>%
#   ProCESS::seq_to_long()
# 
muts <- readRDS('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/validation_subclonal/tables/table_process_SPN05_50x_0.9p_mutect2_ascat.rds') %>% 
  group_by(cluster_id_process, sample_id) %>%
  mutate(is_clonal_process=replace(FALSE, ccf_process > 0.9, TRUE)) %>% 
  ungroup() %>%
  mutate(cluster_id_process_full = cluster_id_process) %>%
  mutate(cluster_id_process = replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal')) %>% 
  select(-ccf_process, -n_clones_process, -n_mutations_process,-coverage, -purity, -patient_id)
  # Parse mutation_id: <sample>:<chr>:<pos>:<ref>
  
muts <- muts %>%
  separate(
    mutation_id,
    into = c("sample_prefix", "chr", "pos", "ref"),
    sep = ":",
    remove = FALSE
  ) %>%
  mutate(
    from = as.integer(pos),
    to   = as.integer(pos)
  ) %>%
  select(
    chr, from, to,
    causes, is_driver_process, driver_label_process,
    sample_id, vaf_process, cluster_id_process
  ) 


dups <- muts |>
  dplyr::summarise(n = dplyr::n(), 
                   .by = c(chr, from, to, causes, is_driver_process, driver_label_process,
                                           sample_id)) |>
  dplyr::filter(n > 1L) %>% 
  mutate(id = paste(chr, from, to,sep = ':'))


wider_muts <- muts %>%
  mutate(id = paste(chr, from, to,sep = ':')) %>% 
  filter(!(id %in% dups$id)) %>% 
  pivot_wider(
  names_from  = sample_id,
  values_from = vaf_process,
  values_fill = 0
)

wider_muts <- wider_muts %>% 
  dplyr::rename(VAF.SPN05_1.2 = SPN05_SPN05_1.2,
                VAF.SPN05_1.1 = SPN05_SPN05_1.1,
                VAF.SPN05_1.3 = SPN05_SPN05_1.3) %>% 
  mutate(driver_label_process = case_when(
    driver_label_process == "SPN05:13:32340437:T" ~ "BRCA2",
    driver_label_process == "SPN05:17:7674221:A" ~ "TP53"
  ))


# define colors
labs <- muts$cluster_id_process %>% unique()
color <- c(
  "plum4", # purple
  "chocolate1", # orange
  "palegreen4", # green
  "steelblue3", # blue
  "palevioletred", # pink
  "goldenrod", # grey
  "sienna4", # brown
  "indianred" # red
)
palette <- setNames(color[seq_along(labs)], labs)

plt_multi <- list()
couples <- combn(c('VAF.SPN05_1.1', 'VAF.SPN05_1.2', 'VAF.SPN05_1.3'), 2, simplify = FALSE)
c = couples[[1]]
for (c in couples){
  s1 = c[[1]]
  s2 = c[[2]]
  
  
  plt_multi[[paste(s1, s2, sep = '-')]] <-  ggplot() + 
    geom_point(data = wider_muts, aes(x = .data[[s1]], y = .data[[s2]], col = cluster_id_process), size = .1, alpha = .5) +
    geom_point(
      data = subset(wider_muts, is_driver_process == T),
      aes(x = .data[[s1]], y = .data[[s2]]),
      color = "black", size = 2, shape = 15
    ) +
    scale_color_manual('', values = palette) +
    ggrepel::geom_label_repel(
      data = subset(wider_muts, is_driver_process == T),
      aes(x = .data[[s1]], y = .data[[s2]], label = driver_label_process, color = cluster_id_process),
      size = 3,
      nudge_y = 0,
      nudge_x = 0,
      show.legend = FALSE
    ) +
    ggplot2::coord_cartesian(clip = 'off') +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    theme_minimal() +
    theme(legend.position = 'bottom')
  
}
p_mult <- wrap_plots(plt_multi,guides = 'collect') &theme(legend.position = 'bottom')

ggsave(plot = p_mult, 
       filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/spn_plot/SPN05/multivariate.png',
       height = 3.5, 
       width = 9, units = 'in', limitsize = F)



plt_marg <- list()
for (s in c('VAF.SPN05_1.1', 'VAF.SPN05_1.2', 'VAF.SPN05_1.3')){
  plt_marg[[s]] <- ggplot(data = wider_muts) +
    geom_histogram(aes(x = .data[[s]], fill = cluster_id_process), binwidth = 0.01, position = 'identity', alpha = .5) +
    scale_fill_manual(values = palette) +
    geom_label_repel(
      data = subset(wider_muts, is_driver_process == TRUE),
      aes(
        x = .data[[s]],
        y = 0,
        label = driver_label_process,
        color = cluster_id_process
      ),
      nudge_y = 3000,        # Push the label up to the 2000 mark
      direction = "y",       # Force movement vertically
      min.segment.length = 0, # Keeps segments even for short distances
      show.legend = FALSE,
      segment.size = 0.25,
      segment.curvature = 0
    ) + 
    scale_color_manual(values = palette) +
    xlim(0.07, 1.01)+
    ylab('') + 
    theme_minimal() +
    theme(legend.position = 'none')
}


p_marg <- wrap_plots(plt_marg, guides = 'collect', nrow = 1) 
ggsave(plot = p_marg, 
       filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/spn_plot/SPN05/marginal.pdf',       
       height = 3.2, 
       width = 9, units = 'in', limitsize = F)

