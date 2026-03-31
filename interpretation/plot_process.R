setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')
source("../validation/SCOUT/colors.R")
source("../figures/figure3/utils_plot.R")

out = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/"
base = '/orfeo/cephfs/scratch/cdslab/shared/SCOUT/validation_subclonal_new/tables/'


spn='SPN04'
cov=100
pur=.9

color_palette_process = RColorBrewer::brewer.pal(n = 8, name = "Dark2") 
names(color_palette_process) <- c('Clonal', paste0('Clone ', 1:7))
color_palette_process['Subclonal'] = 'gray70'

color_palette_process_uni = RColorBrewer::brewer.pal(n = 8, name = "Set1") 
names(color_palette_process_uni) <- c('Clonal', paste0('Clone ', 1:7))
color_palette_process_uni['Subclonal'] = 'gray70'

get_exposure <- function(table){
  
  nmuts <- table %>% 
    group_by(patient_id, coverage, purity, cluster_id_process) %>% 
    summarise(
      tot_nmuts = n_distinct(mutation_id),
      .groups = "drop"
    )
  
  exposure_raw <- table %>% 
    group_by(patient_id, coverage, purity, causes, cluster_id_process) %>% 
    summarise(
      nmuts_cause = n_distinct(mutation_id),
      .groups = "drop"
    )
  
  exposure <- exposure_raw %>% 
    left_join(
      nmuts,
      by = c("patient_id", "coverage", "purity", "cluster_id_process")
    ) %>% 
    mutate(
      exposure = nmuts_cause / tot_nmuts
    )
  return(exposure)
}

plot_exposure <- function(df, sig_type, tool, col, type = '', signature = ''){
  plot <- df %>% 
    ggplot(aes(fill=causes, y=exposure, x=as.factor(cluster_id_process))) + 
    geom_bar(position="fill", stat="identity")+
    coord_flip()+
    scale_fill_manual(values=c(sbs_colors, id_colors))+
    theme_bw()+
    theme(legend.position = "bottom")+
    xlab("Cluster")+
    ylab("Exposures")+
    ggtitle(label = paste(sig_type, type, signature))
  return(plot)
}

true_drivers_table = readRDS(file.path("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/drivers", spn, "process_drivers.rds")) %>% as_tibble()

true_drivers_table = true_drivers_table  %>% rowwise() %>% 
  mutate(code=replace(code, is.na(code), paste0(type, '_', CNA_type, '_', chr, '_', start, '_', end)),
         mutation_id=paste0(spn, ":", chr, ":", start,  ":", alt)) %>% 
  ungroup() %>%
  select(mutation_id, code)

table_multi <- readRDS(paste0(base, 'table_process_', spn, '_', cov, 'x_', pur, 'p_mutect2_ascat.rds')) %>% 
  filter(!str_detect(causes, "errors")) 

table_multi <- table_multi %>% 
  left_join(true_drivers_table, by = 'mutation_id') %>% 
  select(-driver_label_process) %>% 
  mutate(driver_label_process=code)


if(spn=='SPN03'){
  table_multi = table_multi %>% mutate(mutation_id = if_else(
    mutation_id=="SPN03:9:136496197:",
    'SPN03:9:136496196:C',
    mutation_id
  ))
}


n_samples = length(table_multi$sample_id %>% unique())

c = table_multi %>%
  distinct(cluster_id_process, sample_id) %>%   # keep unique pairs
  dplyr::count(cluster_id_process) %>%
  filter(n==n_samples) %>% 
  select(cluster_id_process)

table_multi = table_multi %>% 
  group_by(cluster_id_process) %>% 
  mutate(
    is_clonal_process = if_else(
      dplyr::first(cluster_id_process) %in% c$cluster_id_process, # first returns the first value in the current group
      all(ccf_process > 0.95),
      FALSE
    )
  ) %>% 
  # mutate(is_clonal_process=replace(FALSE, all(ccf_process > 0.95), TRUE)) %>% ungroup() %>% 
  dplyr::mutate(cluster_id_process_full = cluster_id_process) %>% 
  dplyr::mutate(cluster_id_process=replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal')) %>% 
  dplyr::mutate(cluster_id_process=ifelse(cluster_id_process=='Truncal', 'Clonal', cluster_id_process))


table_uni <- readRDS(paste0(base, 'table_process_univariate_w_private_', spn, '_', cov, 'x_', pur, 'p_mutect2_ascat.rds')) %>% 
  filter(!str_detect(causes, "errors"))  %>%   
  left_join(true_drivers_table, by = 'mutation_id') %>% 
  select(-driver_label_process) %>% 
  mutate(driver_label_process=code)

if(spn=='SPN03'){
  table_uni = table_uni %>% mutate(mutation_id = if_else(
    mutation_id=="SPN03:9:136496197:",
    'SPN03:9:136496196:C',
    mutation_id
  ))
}


table_multi_plot <- table_multi %>% select(sample_id, 
                                           mutation_id, 
                                           is_driver_process, 
                                           cluster_id_process, 
                                           vaf_process, 
                                           driver_label_process)

if (spn == 'SPN07'){
  dups <- table_multi_plot %>% 
    ungroup() %>% 
    dplyr::summarise(n = dplyr::n(), .by = c(mutation_id, is_driver_process, cluster_id_process,
                                             driver_label_process, sample_id)) %>% 
    dplyr::filter(n > 1L) 
  table_multi_plot <- table_multi_plot %>% filter(!mutation_id%in%dups$mutation_id)
}

table_wide <- table_multi_plot %>%
  pivot_wider(
    names_from  = sample_id,
    values_from = vaf_process,
    values_fill = 0
  )


plots_pair <- list()
pairs <- combn(unique(table_multi_plot$sample_id), 2, simplify = FALSE)
for (i in 1:length(pairs)){
  s1 <- pairs[[i]][1]
  s2 <- pairs[[i]][2]

  p_process = table_wide %>%
    ggplot(aes(x =.data[[s1]], y = .data[[s2]], color=cluster_id_process))+
    geom_point( alpha=0.4, size = .7)+
    xlim(0,1)+
    ylim(0,1)+
    xlab(pairs[[i]][1]) +
    ylab(pairs[[i]][2])+
    theme_bw() +
    scale_color_manual('ProCESS clusters', values = color_palette_process) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1) ) )

  p_process = p_process + ggrepel::geom_label_repel(
    data = table_wide %>% filter(is_driver_process == TRUE),
    aes(
      x = .data[[s1]],
      y = .data[[s2]],
      label = driver_label_process,
      colour = cluster_id_process,
    ),
    show.legend = F,
    inherit.aes = FALSE,
    size = 3,
    min.segment.length = 0,
    box.padding = 1)

  plots_pair[[i]] <- p_process
}

nc=2

plt_mutlivariate <- wrap_plots(plots_pair,
           ncol = nc,
           guides = 'collect')


table_uni <- table_uni %>% 
  group_by(cluster_id_process, sample_id) %>%
  mutate(is_clonal_process=replace(FALSE, ccf_process > 0.9, TRUE)) %>% 
  ungroup() %>%
  mutate(cluster_id_process_full = cluster_id_process) %>%
  mutate(cluster_id_process = replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal')) %>% 
  separate(cluster_id_process, sep = '_', into = c('cluster_id_process', 'tmp')) %>% 
  select(-tmp) 

samples <- unique(table_uni$sample_id)
plots_sample <- list()
for (s in samples){
  plt <- table_uni %>% 
    filter(sample_id == s) %>% 
    ggplot() +
    geom_histogram(aes(x = vaf_process, fill = cluster_id_process), alpha=0.6, binwidth = 0.01, position = "identity") +
    xlim(0.05,1.01)+
    xlab(paste0('VAF ', s)) + 
    my_ggplot_theme() +
    scale_fill_manual('ProCESS clusters', values = color_palette_process_uni) 
  
  plots_sample[[s]] <- plt + ggrepel::geom_label_repel(
    data = table_uni %>% filter(is_driver_process == TRUE) %>% filter(sample_id == s) %>% filter(vaf_process>0),
    aes(
      x = vaf_process,
      y = 10,
      label = driver_label_process,
      colour = cluster_id_process, 
    ),
    show.legend = F,
    inherit.aes = FALSE,
    size = 3,
    min.segment.length = 0,
    box.padding = 1) +
    scale_color_manual('ProCESS clusters', values = color_palette_process_uni)
}
plt_univariate <- wrap_plots(plots_sample, ncol = 2, guides = 'collect')

table <- table_multi #%>% filter(cluster_id_process %in% n_muts$cluster_id_process)

table_sbs <- table %>% 
  filter(str_detect(causes, "SBS"))

table_id <- table %>% 
  filter(str_detect(causes, "ID"))

exposure_sbs <- get_exposure(table_sbs)
exposure_id <- get_exposure(table_id)

plt_sbs <- plot_exposure(df = exposure_sbs, sig_type = '')
plt_id <- plot_exposure(df = exposure_id, sig_type = '')

plt <- wrap_plots(plt_univariate,
       plt_mutlivariate, 
       wrap_plots(plt_sbs, plt_id, nrow=1),
       nrow = 3, design = 'A\nA\nA\nB\nB\nB\nB\nB\nC')

ggsave(plot = plt,
       filename = paste0(out, 'plot_process/', spn, '.png'),
       dpi = 200,
       width = 8,
       height = 12)#15
