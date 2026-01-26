setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

indir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/assing_signature/"

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN02'),
                    make_option(c("--purity"), type = "double", default = 0.9),
                    make_option(c("--coverage"), type = "integer", default = 50),
                    make_option(c("--cna_caller"), type = "character", default = 'ascat'),
                    make_option(c("--vcf_caller"), type = "character", default = 'mutect2'),
                    make_option(c("--signature"), type = "character", default = 'BASCULE')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

spn = opt$spn_id
cov = opt$coverage
pur = opt$purity
signature_tool =  opt$signature
cna_caller = opt$cna_caller
mut_caller = opt$vcf_caller

tool = 'viber'

# new table
s1 <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/50x_0.9p_mutect2_ascat_process_driver/driver_annotation/annotate_driver/SCOUT/SPN02/SPN02_SPN02_1.1/SCOUT_SPN02_SPN02_SPN02_1.1_driver.rds')$SPN02_SPN02_1.1$mutations
s2 <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/50x_0.9p_mutect2_ascat_process_driver/driver_annotation/annotate_driver/SCOUT/SPN02/SPN02_SPN02_1.2/SCOUT_SPN02_SPN02_SPN02_1.2_driver.rds')$SPN02_SPN02_1.2$mutations
s1 <- s1 %>% 
  select(chr, from, to, ref, alt, SYMBOL, VAF, is_driver) %>% 
  filter(VAF > 0) %>% 
  mutate(mut_id = paste(chr, from, to, ref, alt, sep = ':'))

s2 <- s2 %>% 
  select(chr, from, to, ref, alt, SYMBOL, VAF, is_driver) %>% 
  filter(VAF > 0) %>% 
  mutate(mut_id = paste(chr, from, to, ref, alt, sep = ':'))

join <- full_join(s1, s2, by = join_by(chr, from, to, ref, alt, mut_id, SYMBOL), suffix = c('_SPN02_1.1', '_SPN02_1.2')) %>% 
  mutate(is_driver = ifelse(is_driver_SPN02_1.1 == T | is_driver_SPN02_1.2 == T, T, F)) %>% 
  mutate(driver_label = ifelse(is_driver == T, SYMBOL, '')) %>% 
  mutate(across(starts_with("VAF_"), ~ tidyr::replace_na(.x, 0))) %>% 
  mutate(across(starts_with("is_driver_"), ~ tidyr::replace_na(.x, FALSE))) %>% 
  mutate(across(starts_with("driver_label"), ~ tidyr::replace_na(.x, ''))) 
  

driver <- unique(join$driver_label)
color_palette_process = c( 'gray50', RColorBrewer::brewer.pal(max(3,  length(driver)-1), name = "Dark2"))%>%
  setNames(driver)

spn02_process <- join %>% 
  ggplot() + 
    geom_point(aes(x = VAF_SPN02_1.1, y = VAF_SPN02_1.2, col = driver_label), size = .5, alpha = .5) + 
    ggrepel::geom_label_repel(
          data = join %>% filter(is_driver == TRUE),
          aes(
            x = VAF_SPN02_1.1,
            y = VAF_SPN02_1.2,
            label = driver_label,
            colour = driver_label,
          ),
          max.overlaps = 50,
          show.legend = F,
          inherit.aes = FALSE,
          size = 3,
          min.segment.length = 0,
          box.padding = 1) +
  scale_color_manual('Driver ProCESS', values = color_palette_process) +
  theme_minimal()

ggsave(plot = spn02_process, filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/spn02_driver_process.png', dpi = 200, width = 8, height = 7, units = 'in')


# old table
s1 <- readRDS('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN02/tumourevo/50x_0.9p_mutect2_ascat/driver_annotation/annotate_driver/SCOUT/SPN02/SPN02_SPN02_1.1/SCOUT_SPN02_SPN02_SPN02_1.1_driver.rds')$SPN02_SPN02_1.1$mutations
s2 <- readRDS('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN02/tumourevo/50x_0.9p_mutect2_ascat/driver_annotation/annotate_driver/SCOUT/SPN02/SPN02_SPN02_1.2/SCOUT_SPN02_SPN02_SPN02_1.2_driver.rds')$SPN02_SPN02_1.2$mutations
s1 <- s1 %>% 
  select(chr, from, to, ref, alt, SYMBOL, VAF, is_driver) %>% 
  filter(VAF > 0) %>% 
  mutate(mut_id = paste(chr, from, to, ref, alt, sep = ':'))

s2 <- s2 %>% 
  select(chr, from, to, ref, alt, SYMBOL, VAF, is_driver) %>% 
  filter(VAF > 0) %>% 
  mutate(mut_id = paste(chr, from, to, ref, alt, sep = ':'))

join <- full_join(s1, s2, by = join_by(chr, from, to, ref, alt, mut_id, SYMBOL), suffix = c('_SPN02_1.1', '_SPN02_1.2')) %>% 
  mutate(is_driver = ifelse(is_driver_SPN02_1.1 == T | is_driver_SPN02_1.2 == T, T, F)) %>% 
  mutate(driver_label = ifelse(is_driver == T, SYMBOL, '')) %>% 
  mutate(across(starts_with("VAF_"), ~ tidyr::replace_na(.x, 0))) %>% 
  mutate(across(starts_with("is_driver_"), ~ tidyr::replace_na(.x, FALSE))) %>% 
  mutate(across(starts_with("driver_label"), ~ tidyr::replace_na(.x, ''))) 


driver <- unique(join$driver_label)
color_palette_process = c( 'gray50', RColorBrewer::brewer.pal(9, name = "Set1"), RColorBrewer::brewer.pal(8, name = "Dark2"))%>%
  setNames(driver)

spn02_tumourevo <- join %>% 
  ggplot() + 
  geom_point(aes(x = VAF_SPN02_1.1, y = VAF_SPN02_1.2, col = driver_label), size = .5, alpha = .5) + 
  ggrepel::geom_label_repel(
    data = join %>% filter(is_driver == TRUE),
    aes(
      x = VAF_SPN02_1.1,
      y = VAF_SPN02_1.2,
      label = driver_label,
      colour = driver_label,
    ),
    max.overlaps = 50,
    show.legend = F,
    inherit.aes = FALSE,
    size = 3,
    min.segment.length = 0,
    box.padding = 1) +
  scale_color_manual('Driver ProCESS', values = color_palette_process) +
  theme_minimal()


ggsave(plot = spn02_tumourevo, filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/spn02_driver_tumourevo.png', dpi = 200, width = 8, height = 7, units = 'in')


process <- readRDS('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN02/sequencing/tumour/purity_0.9/data/mutations/seq_results_muts_merged_coverage_50x.rds') %>% 
  filter(classes != 'germinal')

process <- process %>% 
  #dplyr::mutate(label = ifelse(causes %in% c('Clone 3', 'Clone 2', 'Clone 1'), causes, '')) %>% 
  dplyr::mutate(label = ifelse(label == 'Clone 1', 'BRAF', label)) %>% 
  dplyr::mutate(label = ifelse(label == 'Clone 2', 'PIK3CA', label)) %>% 
  dplyr::mutate(label = ifelse(label == 'Clone 3', 'POLE', label)) 

driver <- unique(process$label) 
color_palette_process = c( 'gray50', RColorBrewer::brewer.pal(9, name = "Set1"), RColorBrewer::brewer.pal(8, name = "Dark2"))%>%
  setNames(driver)

spn02_true_process <- process %>% 
  ggplot() + 
  geom_point(aes(x = SPN02_1.1.VAF, y = SPN02_1.2.VAF, col = label), size = .5, alpha = .5) + 
  ggrepel::geom_label_repel(
    data = process %>% filter(classes == 'driver'),
    aes(
      x = SPN02_1.1.VAF,
      y = SPN02_1.2.VAF,
      label = label,
      colour = label,
    ),
    max.overlaps = 50,
    show.legend = F,
    inherit.aes = FALSE,
    size = 3,
    min.segment.length = 0,
    box.padding = 1) +
  scale_color_manual('ProCESS', values = color_palette_process) +
  theme_minimal()
  
ggsave(plot = spn02_true_process, filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/spn02_process.png', dpi = 200, width = 8, height = 7, units = 'in')
