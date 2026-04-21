setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
library(ggalluvial)
library(purrr)
library(stringr)
library(ggtext)  

source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')
source("../validation/SCOUT/colors.R")
source("../figures/figure3/utils_plot.R")
source("../figures/figure3/utils.R")


base_subclonal = '/orfeo/cephfs/scratch/cdslab/shared/SCOUT/validation_subclonal/tables/'

option_list <- list(make_option(c("--cna_caller"), type = "character", default = 'ascat'),
                    make_option(c("--vcf_caller"), type = "character", default = 'mutect2')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

vcf_caller = opt$vcf_caller #"mutect2"
cna_caller = opt$cna_caller #"sequenza"
base_int = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/interpretation_', vcf_caller, "_", cna_caller)
out = paste0(base_int, '/plot_performance/')
dir.create(out, showWarnings = F, recursive = T)

color_palette_process = RColorBrewer::brewer.pal(n = 8, name = "Dark2") 
names(color_palette_process) <- c('Clonal', paste0('Clone ', 1:7))
color_palette_process['Subclonal'] = 'gray70'

color_class <- c("palegreen4","darkseagreen", "goldenrod", "indianred", "gainsboro")
names(color_class) <- c('Clonal', 'Strong\nExpansion', 'Medium\nExpansion', 'Low\nExpansion', 'Neutral\nExpansion')

color_score <- c('plum4', 'darkseagreen4', 'cadetblue4', 'salmon2')
score_types <- c("all", 'no_tail', 'no_driver', 'no_sign')
names(color_score) <- score_types

colors_cluster = c('indianred', 
                   'steelblue', 
                   'forestgreen', 
                   'goldenrod', 
                   'darkorange3', 
                   'palevioletred', 
                   'mediumpurple', 
                   'cornsilk4', 
                   'olivedrab3', 
                   'steelblue4', 
                   'indianred4',
                   'aquamarine3',
                   'saddlebrown',
                   'deeppink2',
                   'cornflowerblue',
                   'black')
names(colors_cluster) = paste0('C',0:15)

coverage_list = c(50, 100, 150)
purity_list = c(0.9, 0.6, 0.3)
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04','SPN05', 'SPN06', 'SPN07')
tool_list = c( 'viber', 'pyclonevi')
sig_tool_list = c('BASCULE', 'SigProfiler')

combs = expand.grid(spn = spn_list,
                    coverage=coverage_list,
                    purity=purity_list,
                    tool=tool_list,
                    sig_tool = sig_tool_list)

df_all <- lapply(1:nrow(combs), FUN = function(i){
  cov = combs[i, "coverage"]
  pur = combs[i, "purity"]
  s_tool = combs[i, "sig_tool"]
  tool = combs[i, "tool"]
  sp = combs[i, "spn"]
  
  df_tool = readRDS(paste0(base_int, '/res_int/', cov, 'x_', pur, 'p_', s_tool,'_',tool, '_df.rds')) %>% 
    filter(spn == sp) 
  
  return(df_tool)
}) %>% bind_rows()


sp = 'SPN04'
cov = 50
pur = 0.9
tool = 'viber'
s_tool = 'SigProfiler'

df_long <- df_all %>% 
  filter(spn == sp, 
         dec_tool == tool,
         sig_tool == s_tool | is_clonal_tool == T) %>% 
  pivot_longer(
    cols = c(score_driver, score_all, score_tail,
            # score_no_driver, score_no_tail, score_no_sign, 
            score_sign),
    names_to  = "score_type",
    values_to = "score"
  ) %>% 
  mutate(
    score_type = factor(score_type, levels = c(
      "score_driver", "score_tail", "score_sign",
      #"score_no_driver", "score_no_tail", "score_no_sign", 
      "score_all"
    )),
    purity   = factor(purity),
    coverage = factor(coverage)
  ) %>% 
  mutate(
    driver_status = case_when(
      contains_driver_process &  contains_driver_tool ~ "both",
      contains_driver_process & !contains_driver_tool ~ "process only",
      !contains_driver_process &  contains_driver_tool ~ "tool only",
      TRUE                                             ~ "neither"
    ) %>% factor(levels = c("neither", "process only", "tool only", "both"))
  ) %>% 
  mutate(
    cluster = factor(
      cluster,
      levels = stringr::str_sort(unique(cluster), numeric = TRUE)
    )
  )


score_labels <- c(
  score_driver    = "Driver",
  score_tail      = "Tail",
  score_sign      = "Signature",
  # score_no_driver = "Signature + Tail",
  # score_no_tail   = "Driver + Signature",
  # score_no_sign   = "Driver + Tail",
  score_all       = "All"
)

ggplot(df_long, aes(x = purity, y = coverage)) +
  
  # filled tiles: score value + border colour from driver_status
  geom_tile(
    aes(fill = score, colour = driver_status),
    linewidth = 1
  ) +
  
  # score value text
  geom_text(
    aes(label = round(score, 2)),
    size   = 2, colour = "white"
  ) +
  
  facet_grid(
    score_type ~ cluster,
    labeller = labeller(score_type = score_labels)
  ) +
  
  # colour scale: driver status border
  scale_colour_manual(
    name   = "Driver",
    values = c(
      "neither"      = "white",
      "process only" = "steelblue",
      "tool only"    = "mediumpurple4",
      "both"         = "goldenrod"
    ),
    # show solid-border swatches in the legend (not score-coloured)
    guide = guide_legend(
      override.aes = list(fill = "grey80", linewidth = 1.2)
    )
  ) +
  
  # clonal label inside tile
  # geom_text(data = ~ filter(.x, is_clonal_tool),
  #           aes(label = "C"),
  #           size = 2.8, colour = "white", fontface = "bold") +
  
  # score value text (optional — remove if tiles are small)
  geom_text(aes(label = round(score, 2)),
            size = 2, colour = "black") +
  
  facet_grid(score_type ~ cluster,
             labeller = labeller(score_type = score_labels)) +
  
  scale_fill_distiller(
    palette  = "RdYlGn",
    direction = 1,
    limits   = c(0, 1),
    name     = "Score"
  ) +
  
  labs(x = "Purity", y = "Coverage") +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    #strip.text.x = element_text(face = "bold"),   # SPN on top
    strip.text.y = element_text(angle = 0, hjust = 0), # score on right, readable
    legend.position = "bottom",
    axis.text    = element_text(size = 8)
  ) +
  ggtitle(paste(sp, tool, s_tool))
