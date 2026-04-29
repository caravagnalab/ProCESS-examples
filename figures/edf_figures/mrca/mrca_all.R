setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils.R")
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/SCOUT/colors.R')


spn04 <- readRDS("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN04/process/treatment_info.rds")
spn06 <- readRDS("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN06/process/chemo_timing.rds")
treatment = tibble(spn = c('SPN04', 'SPN06', 'SPN06', 'SPN07'),
                   start = c(302.4552, 289.873, 367.6605, 158.128),
                   end = c(344.6101, 330.5606, 401.4993, 173))

all_abs = lapply(paste0('SPN0', 1:7), FUN = function(sp){
  forest = load_sample_forest(get_sample_forest(spn = sp))
  sample_time = forest$get_samples_info() %>% select(name, time) %>% dplyr::rename(sample = name, sample_time = time)
  sample_names = forest$get_samples_info() %>% dplyr::pull(.data$name) 
  MRCAs <- lapply(sample_names, function(s) {
    forest$get_coalescent_cells(forest$get_nodes() %>% 
                                  dplyr::filter(sample %in% s) %>% dplyr::pull(.data$cell_id)) %>% 
      dplyr::mutate(sample = s)
  }) %>% Reduce(f = dplyr::bind_rows) %>% dplyr::group_by(.data$cell_id) %>% 
    dplyr::mutate(cell_id = paste(.data$cell_id)) %>% 
    ungroup() %>% 
    select(sample, birth_time)
  df = sample_time %>% left_join(MRCAs) %>% mutate(spn = sp)
  return(df)
  }
) %>% bind_rows()

all <- all %>% left_join(treatment)

treatment_norm <- all %>%
  filter(!is.na(start)) %>%
  distinct(spn, start, end) %>%
  group_by(spn) %>%
  mutate(
    spn_max    = max(all$sample_time[all$spn == cur_group()$spn]),
    norm_start = start / spn_max,
    norm_end   = end   / spn_max
  ) %>%
  ungroup()

# Normalise per SPN
all <- all_abs %>%
  mutate(start_time = 0) %>% 
  group_by(spn) %>%
  mutate(
    norm_sample = sample_time / max(sample_time),
    norm_birth  = birth_time  / max(sample_time)
  ) %>%
  ungroup()

all <- all %>%
  mutate(sub_sample = sub("^SPN\\d+_", "", sample))

# Pivot to long format so both point types share the same colour aesthetic
long <- all %>%
  pivot_longer(
    cols      = c(norm_sample, norm_birth),
    names_to  = "event_type",
    values_to = "norm_time"
  ) %>%
  mutate(
    event_type = recode(event_type,
                        norm_sample = "Sampling",
                        norm_birth  = "MRCA")
  ) %>% 
  mutate(sub_sample = factor(sub_sample, levels = c('1.1', '1.2', '1.3', '2.1', '2.2', '3.1', '3.2', '4.1')))

# Jitter position — same seed keeps circle/square pairs aligned vertically
pos <- position_jitter(height = 0, width = 0, seed = 42)

pp <- ggplot(long, aes(x = norm_time, y = sub_sample, fill = sub_sample, col = sub_sample)) +
  geom_rect(
    data        = treatment_norm,
    aes(xmin = norm_start, xmax = norm_end, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill        = "coral", alpha = 0.2
  ) +
  geom_point(aes(shape = event_type),
    size = 3, stroke = 0.25,
    position = pos
  ) +
  scale_x_continuous(
    expand = c(0.02, 0)
  ) +
  xlim(-0.001, 1.001) +
  scale_y_discrete(limits = rev) + 
  scale_color_brewer(palette = "Dark2", name = "Sample") +
  scale_fill_brewer(palette = "Dark2", name = "Sample") +
  scale_shape_manual(name = "Time-point", values = c(21,22)) +
  labs(
    x        = "Normalised time",
    y        = NULL
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 22, size = 3)
    )
  ) +
  my_ggplot_theme() + 
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  facet_grid(spn~., scales = 'free_y') 

pp
ggsave(filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure5/mrca_timeline.pdf',
       plot = pp, 
       width = 3.5, 
       height = 4.5, 
       units = 'in' )
