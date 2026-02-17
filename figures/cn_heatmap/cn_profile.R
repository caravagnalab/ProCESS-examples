library(dplyr)
library(tidyverse)
library(GenomicRanges)
library(circlize)
library(ComplexHeatmap)
library(ProCESS)
library(EnrichedHeatmap)

source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/report/plotting/utils.R")
setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/cn_heatmap/')
source('../../getters/process_getters.R')
source('utils.R')

upper <- 0.89
lower <- 0.11

cna_data_all <- tibble()
for (spn in c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN05', 'SPN06', 'SPN07')){
  samples <- get_sample_names(spn)
  cna_data <- lapply(samples, FUN = function(s){
    readRDS(get_process_cna(spn = spn, sample = s)) %>% mutate(sample = s)
  }) %>% bind_rows()
  
  x <- cna_data %>%
    dplyr::rename("from"="begin") %>% 
    dplyr::rename("to"="end") %>% 
    dplyr::rename(Major=major) %>% 
    dplyr::rename(sample_id=sample) %>% 
    mutate(segment_id=paste(chr,from,to,sep=":")) %>% 
    mutate(CN_type = ifelse(ratio < upper & ratio > lower, 'sub-clonal', 'clonal')) %>% 
    filter(!(CN_type == 'sub-clonal' & ratio  > upper)) %>% 
    filter(!(CN_type == 'clonal' & ratio <= lower))  %>% 
    mutate(ratio = ifelse(ratio < upper & ratio > lower, ratio, 1)) %>% 
    mutate(ratio = round(ratio, digits=1)) %>% 
    filter(ratio != 0) 
  
  cna_data_all <- bind_rows(cna_data_all, x)
}


out <- segment_fixed_windows_all(cna_data_all, window_size = 5e5)# %>% 
  #mutate(sample_id = factor(sample_id, levels = samples_levels)) %>% 
  #arrange(sample_id, desc(ratio))

segments <- out %>% 
  filter(ratio ==1) %>% 
  mutate(CNA = dplyr::case_when(
    minor == 0 & Major == 1 ~ "Deletion",
    minor == 0 & Major == 2 ~ "Deletion",
    minor == 0 & (Major > 2 | Major == 0) ~ "Deletion",
    minor > 0 & ((Major + minor) >= 3) ~ "Gain",
    TRUE ~ "None"
  )
)

segments_binned = segments %>%
  dplyr::group_by(chr, from, to, CNA) %>%
  distinct() %>% 
  dplyr::summarise(n = n()) %>%
  ungroup() %>% 
  dplyr::mutate(n = ifelse(CNA == 'Deletion',-1 * n, n)) %>% 
  filter(CNA != 'None') %>% 
  filter(chr %in% 1:22) %>% 
  filter(n <= 24)

relative_to_absolute_coordinates <- function(data){
  reference_genome = CNAqc:::get_reference('GRCh38')
  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr
  data %>% mutate(from = from + vfrom[chr], to = to + vfrom[chr])
}
segments_binned_abs = relative_to_absolute_coordinates(segments_binned %>% 
                                                         mutate(chr = paste0('chr',chr)))




reference_coordinates = CNAqc:::get_reference('GRCh38') %>% filter(chr %in% paste0('chr',1:22))
low = min(reference_coordinates$from)
upp = max(reference_coordinates$to)

driver <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/cn_heatmap/driver.rds')
driver <- tibble(gene = unique(driver))
coords_gene <- CNAqc::gene_coordinates_GRCh38
driver <- left_join(driver, coords_gene) %>% 
  filter(chr %in% paste0('chr',1:22)) %>% 
  relative_to_absolute_coordinates

cn_plot <- ggplot2::ggplot(reference_coordinates) +
  my_ggplot_theme() +
  ggplot2::geom_rect(
    ggplot2::aes(
      xmin = centromerStart,
      xmax = centromerEnd,
      ymin = -Inf,
      ymax = Inf
    ),
    alpha = .3,
    colour = 'gainsboro'
  ) +
  ggplot2::geom_segment(
    data = reference_coordinates,
    ggplot2::aes(
      x = from,
      xend = from,
      y = -Inf,
      yend = Inf
    ),
    size = .1,
    color = 'black',
    linetype = 8
  ) +
  ggplot2::geom_hline(yintercept = 0,
                      size = 1,
                      colour = 'gainsboro') +
  ggplot2::labs(x = "Chromosome", y = "Cases (out of 24 )") +
  ggplot2::scale_x_continuous(
    breaks = c(0, reference_coordinates$from, upp),
    labels = c(
      "",
      gsub(
        pattern = 'chr',
        replacement = '',
        reference_coordinates$chr
      ),
      ""
    )
  ) + 
ggplot2::geom_rect(
  data = segments_binned_abs,
  ggplot2::aes(
    xmin = from,
    xmax = to,
    ymin = 0,
    ymax = n,
    fill = CNA
  ),
  alpha = .8
) +
  ggrepel::geom_label_repel(
    data = driver,
    aes(
      x = from,
      y = 0,
      label = gene
    ),
    inherit.aes = FALSE,
    show.legend = FALSE,
    
    # text
    size = 2.5,
    fill = "white",
    label.size = 0,
    
    # repel behavior
    #direction = "y",
    ylim = c(-15, 20),
    box.padding = 1,
    #point.padding = 0.2,
    max.overlaps = Inf,
    
    # segments
    segment.curvature = 1e-50,
    min.segment.length = 0,
    segment.size = 0.3,
    segment.color = "grey30",
    segment.linetype = 1,
    segment.ncp = 2,
    alpha = .8,
    segment.angle = 90
  ) + 
  ggplot2::scale_y_continuous(labels = c(10,0,10,20,30),breaks = c(-10,0,10,20,30)) + 
  ggplot2::scale_fill_manual(values = c('Gain' = 'indianred', 'Deletion' = 'steelblue')) +
  ggplot2::guides(fill = ggplot2::guide_legend('', override.aes = list(alpha = 1))) +
  my_ggplot_theme()
cn_plot
