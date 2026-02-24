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
sample_names <- get_sample_names('SPN03')
number_of_samples <- length(sample_names)

data <- readRDS(get_mutations(spn = 'SPN03', coverage = 50, purity = 0.9, type = 'tumour'))
normal <- readRDS(get_mutations(spn = 'SPN03', coverage = 30, purity = 1, type = 'normal')) %>%
  filter(classes == 'germinal') %>% 
  filter(chr == 13) %>% 
  seq_to_long()

# germline
germline <- data %>% 
  filter(classes == 'germinal') %>% 
  filter(chr == 13) %>% 
  seq_to_long()

cna <- phylo_forest$get_bulk_allelic_fragmentation() %>% filter(chr == 13)

n_cov = 30
t_cov = 50
ratio = n_cov/t_cov
g_sample <- germline %>%   
  ungroup() %>% 
  left_join(normal %>% ungroup(), join_by(chr,from,to), suffix = c('_T', '_N')) %>% 
  filter(VAF_N > 0.3) %>% 
  filter(VAF_N < 0.7) %>%
  filter(DP_T > 15) %>% 
  sample_n(size = 1e4) %>% 
  mutate(DR = (DP_T/DP_N)*ratio) %>% 
  left_join(cna) %>% 
  filter(from >= begin, to <= end) %>% 
  mutate(CN = paste(major, minor, sep =':'))

col <- c('1:1'='slategray4', '1:0'='goldenrod3')
germline <- g_sample %>% 
  ggplot() +
  geom_point(aes(x = from, y = VAF_T, col = CN), size = .3) +
  facet_grid(.~sample_name_T) +
  scale_color_manual(values = col)+
  ylab('BAF') +
  xlab('chr13') +
  ylim(0,1) + 
  guides(color = guide_legend(override.aes = list(size = 3))) +
  
g_sample %>% 
  ggplot() +
  geom_point(aes(x = from, y = DR, , col = CN), size = .3) +
  facet_grid(.~sample_name_T) +
  scale_color_manual(values = col)+
  ylab('DR') +
  xlab('chr13') +
  ylim(0,3) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  plot_layout(nrow = 2, guides = 'collect') & theme_minimal()

germline


somatic <- data %>%
  dplyr::filter(classes!="germinal") %>%
  dplyr::filter(!stringr::str_detect(causes, 'errors')) %>%
  ProCESS::seq_to_long()

somatic <- somatic %>%
  tidyr::pivot_wider(
    names_from = sample_name,
    values_from = c(NV, DP, VAF),
    names_glue = "{sample_name}.{.value}"
  )

mutations <- somatic %>% 
  tidyr::as_tibble() %>%
  dplyr::mutate(mutation_id=paste(chr, from, to, ref, alt, sep="_"), chr_pos=from)%>%
  dplyr::select(mutation_id, everything())

driver <- phylo_forest$get_driver_mutations() %>% filter(type == 'SID') %>% select(chr, start, end, code) %>% separate(code, sep = ' ', into = c('gene', 'code'))
mutations <- mutations %>% left_join(driver %>% rename(from = start))

mutations <- mutations %>%
  dplyr::filter(rowSums(dplyr::select(., ends_with(".NV")) != 0) > 0)

model_df <- mutations %>%
  dplyr::filter(!stringr::str_detect(causes, 'errors')) %>%
  dplyr::select(mutation_id, contains("VAF"), classes, gene) %>% 
  mutate(classes = ifelse(classes == 'driver', 'driver', NA))

mutation_ids <- unique(model_df$mutation_id)


labels <- c()
for (mut in mutation_ids){
  muts_data <- model_df %>% 
    filter(mutation_id==mut) %>% 
    select(contains("VAF"))
  
  no_zeros <- muts_data[which(muts_data != 0)]
  if (ncol(no_zeros)>1){
    samples <- gsub(pattern = ".VAF",replacement = "",x = colnames(no_zeros))
    label <- paste0(samples, collapse = "-")
    label <- label #paste0("SHARED_",label)
    labels <- c(labels,label)
  } else if (ncol(no_zeros)==1){
    samples <- gsub(pattern = ".VAF",replacement = "",x = colnames(no_zeros))
    label <- paste0(samples, collapse = "-")
    label <- label #paste0("PRIVATE_",label)
    labels <- c(labels,label)
  }
}
model_df$label <- labels

# define colors
labs <- model_df$label %>% unique()
color <- c(
  "sienna4", # brown
  "indianred", # red
  "plum4", # purple
  "chocolate1", # orange
  "palegreen4", # green
  "steelblue3", # blue
  "palevioletred", # pink
  "goldenrod"  # grey
)
palette <- setNames(color[seq_along(labs)], labs)

plt_multi <- list()
couples <- combn(paste0(sample_names, '.VAF'), 2, simplify = FALSE)
for (c in couples){
  s1 = c[[1]]
  s2 = c[[2]]
  

  plt_multi[[paste(s1, s2, sep = '-')]] <-  ggplot() + 
    geom_point(data = model_df, aes(x = .data[[s1]], y = .data[[s2]], col = label), size = .3) +
    geom_point(
      data = subset(model_df, classes == "driver"),
      aes(x = .data[[s1]], y = .data[[s2]]),
      color = "black", size = 2, shape = 15
    ) +
    scale_color_manual('', values = palette) +
    ggrepel::geom_label_repel(
      data = subset(model_df, classes == "driver"),
      aes(x = .data[[s1]], y = .data[[s2]], label = gene, color = label),
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
       filename = '/orfeo/cephfs/scratch/area/lvaleriani/tesi/spn03/multivariate.png',height = 5.5, width = 12, units = 'in', limitsize = F)

plt_marg <- list()
for (s in paste0(sample_names, '.VAF')){
  plt_marg[[s]] <- ggplot(data = model_df) +
    geom_histogram(aes(x = .data[[s]], fill = label), binwidth = 0.02, position = 'identity', alpha = .5) +
    scale_fill_manual(values = palette) +
    ggrepel::geom_label_repel(
      data = subset(model_df, classes == "driver"),
      aes(x = .data[[s]], y = 1, label = gene, color = label),
      #color = 'black',
      ylim = c(300 * .7, 300 * .9),
      size = 3,
      nudge_x = 0,
      show.legend = FALSE,
      segment.size = 0.25,
      segment.linetype = 2,
      segment.curvature = 1,
      segment.ncp = 1,
      segment.square = TRUE,
      segment.inflect = TRUE
    ) + 
    scale_color_manual(values = palette) +
    xlim(0.02, 1.01)+
    ylab('') + 
    theme_minimal() +
    theme(legend.position = 'none')
}


p_marg <- wrap_plots(plt_marg, guides = 'collect', nrow = 2) 


p_final <- wrap_plots(p_marg, p_mult, design = 'B\nB\nA', guides = 'collect')
ggsave(plot = p_final, 
       filename = '/orfeo/cephfs/scratch/area/lvaleriani/tesi/spn03/all.png',height = 8, width = 12, units = 'in', limitsize = F)
