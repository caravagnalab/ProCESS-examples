library(ggplot2)
library(tidyverse)
library(patchwork)

reticulate::use_condaenv("/orfeo/scratch/cdslab/ggandolfi/miniconda/envs/bascule-env")
py = reticulate::import("pybascule")
library(bascule)

setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/')
source('getters/tumourevo_getters.R')
source('getters/process_getters.R')
source('figures/figure3/utils_plot.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/SCOUT/colors.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/utils_plot.R')

base = "/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/"
y_breaks <- function(x) {
  x <- x[x != 0]                # exclude 0
  if(length(x) == 0) return(NULL)
  max_val <- max(x, na.rm = TRUE) - 0.5
  half_val <- max_val / 2
  max_val <- round(max_val, 0)
  half_val <- round(half_val, 0)
  return(c(half_val, max_val))
}


spn = 'SPN03'
data_join_all <- lapply(paste0('SPN0',1:7),FUN = function(spn){
  
  data_sbs = read.table(paste0('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/sign/', 
                               spn, 
                               '/output/SBS/',spn, '.SBS96.all'), header = T, sep = '\t')
  
  rownames(data_sbs) <- data_sbs$MutationType
  data_sbs <- data_sbs %>% select(!MutationType)
  data_sbs_counts <- t(data_sbs) %>% as.data.frame()
  rownames(data_sbs_counts) <- sub("^[^_]+_(.*)", "\\1", rownames(data_sbs_counts))
  data_sbs_counts = data_sbs_counts[rowSums(data_sbs_counts) > 0,]
  
  data <- data_sbs_counts %>%
    as.data.frame() %>%
    rownames_to_column("samples") %>%
    pivot_longer(
      cols = -samples,
      names_to = "features",
      values_to = "value"
    ) %>% 
    group_by(samples) %>%
    mutate(
      percentage = 100 * value / sum(value, na.rm = TRUE)
    ) %>%
    ungroup() %>% 
    mutate(spn = spn)
  
  process_data <- readRDS(paste0('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/sign/', 
                                 spn, 
                                 '/mut_class.rds')) %>% 
    separate(mut_id, into = c('chr', 'pos', 'ref', 'alt'), sep = ':', remove = F) %>% 
    filter(chr %in% as.character(1:22))
  
  
  muts_data <- lapply(1:22, FUN = function(chr){
    tmp_file <- read.table(paste0('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/sign/', 
                                  spn, '/output/vcf_files/SNV/',chr, '_seqinfo.txt'), header = F, sep = '\t')
  }) %>% bind_rows()
  colnames(muts_data) <- c('sample', 'chr', 'pos', 'context', '_')
  
  muts_data$features <- sub(
    ".*([ACGT])\\[([ACGT]>[ACGT])\\]([ACGT]).*",
    "\\1[\\2]\\3",
    muts_data$context
  )
  
  muts_data <- muts_data %>% select(sample, chr, pos, features) %>% distinct()
  muts_data %>% 
    left_join(process_data %>% mutate(pos = as.integer(pos), chr = as.integer(chr))) %>% 
    select(sample, chr, pos, features, causes) %>% 
    group_by(sample, features, causes) %>% 
    summarise(value = n()) %>% 
    group_by(sample) %>%
    mutate(
      percentage = 100 * value / sum(value, na.rm = TRUE)
    ) %>%
    ungroup() %>% 
    mutate(spn = spn) %>% 
    bascule:::reformat_contexts(what = "SBS") %>%
    separate(variant, into = c("v1", "v2"), sep = ">", remove = FALSE) %>% 
    separate(context, into = c("c1", "c2"), sep = "_", remove = FALSE) %>% 
    mutate(context = paste0(c1, v1, c2)) %>%  
    separate(sample, c("spn", "sample_id"), sep = "_")
}) %>% bind_rows()


plt_sbs <- ggplot(data_join_all, 
                  aes(x = context, y = percentage, fill = causes)) +
  geom_col() +
  ggh4x::facet_nested(spn + sample_id~variant, 
                      scales="free", 
                      space="free_x") + 
  ylab("Number of SBS") +
  #theme_bw() + 
  my_ggplot_theme() +
  scale_y_continuous(breaks = y_breaks) +  # control y ticks
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size=6),
    axis.text.y=element_text(size=5),
    panel.spacing.x = unit(0.05, "in"), 
    panel.spacing.y = unit(0.01, "in"),
    strip.text.x = element_text(size = 8, margin = margin(t = 0, b = 0), colour = 'gray20'),
    strip.text.y = element_text(size = 8, margin = margin(l = 0, r = 0), colour = 'gray20'),
    strip.background = element_rect(colour = "gray70", fill = "gainsboro"),
    #panel.background = element_blank(),      # removes grey background
    panel.grid.major = element_blank(),      # removes major grid lines
    panel.grid.minor = element_blank(),      # removes minor grid lines
    #panel.border = element_rect(color = "black", fill = NA)
  ) + scale_fill_manual(values = sbs_colors)
plt_sbs
ggsave(plot = plt_sbs, filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/sbs_new.pdf', 
       width = 10, height = 7, units = 'in')

# ID ####
data_join_id <- lapply(paste0('SPN0',1:7),FUN = function(spn){
  
  data_sbs = read.table(paste0('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/sign/', 
                               spn, 
                               '/output/ID/',spn, '.ID83.all'), header = T, sep = '\t')
  
  rownames(data_sbs) <- data_sbs$MutationType
  data_sbs <- data_sbs %>% select(!MutationType)
  data_sbs_counts <- t(data_sbs) %>% as.data.frame()
  rownames(data_sbs_counts) <- sub("^[^_]+_(.*)", "\\1", rownames(data_sbs_counts))
  data_sbs_counts = data_sbs_counts[rowSums(data_sbs_counts) > 0,]
  
  data <- data_sbs_counts %>%
    as.data.frame() %>%
    rownames_to_column("samples") %>%
    pivot_longer(
      cols = -samples,
      names_to = "features",
      values_to = "value"
    ) %>% 
    group_by(samples) %>%
    mutate(
      percentage = 100 * value / sum(value, na.rm = TRUE)
    ) %>%
    ungroup() %>% 
    mutate(spn = spn)
  
  process_data <- readRDS(paste0('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/sign/', 
                                 spn, 
                                 '/mut_class.rds')) %>% 
    separate(mut_id, into = c('chr', 'pos', 'ref', 'alt'), sep = ':', remove = F) %>% 
    filter(chr %in% as.character(1:22))
  
  
  muts_data <- lapply(1:22, FUN = function(chr){
    tmp_file <- read.table(paste0('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/sign/', 
                                  spn, '/output/vcf_files/ID/',chr, '_seqinfo.txt'), header = F, sep = '\t')
  }) %>% bind_rows()
  colnames(muts_data) <- c('sample', 'chr', 'pos', 'context', '_')
  
  muts_data$features <- sub("^[^:]*:", "", muts_data$context)
  
  muts_data <- muts_data %>% select(sample, chr, pos, features) %>% distinct()
  muts_data %>% 
    left_join(process_data %>% mutate(pos = as.integer(pos), chr = as.integer(chr))) %>% 
    select(sample, chr, pos, features, causes) %>% 
    group_by(sample, features, causes) %>% 
    summarise(value = n()) %>% 
    group_by(sample) %>%
    mutate(
      percentage = 100 * value / sum(value, na.rm = TRUE)
    ) %>%
    ungroup() %>% 
    mutate(spn = spn) %>% 
    bascule:::reformat_contexts(what = "ID") %>%
    separate(col = variant, into = c('N', 'Type', 'Value')) %>% 
    mutate(name = paste(N, Type, sep = ' ')) %>% 
    separate(sample, c("spn", "sample_id"), sep = "_") %>% 
    mutate(name = ifelse(Value == 'R' & Type == 'Del', '> 1bp Deletion', name)) %>%
    mutate(name = ifelse(Value == 'R' & Type == 'Ins', '> 1bp Insertion', name)) %>% 
    mutate(Value = ifelse(Value == 'R', N, Value)) %>% 
    mutate(name = ifelse(Value == 'M', 'Microhomology', name)) %>% 
    mutate(Value = ifelse(Value == 'M', N, Value)) %>% 
    mutate(name = case_when(
      name == '1 Del' ~ '1bp Deletion',
      name == '1 Ins' ~ '1bp Insertion',.default = name
    ))
}) %>% bind_rows()


data_join_id$name <- factor(data_join_id$name, levels = c('1bp Deletion', '1bp Insertion', '> 1bp Deletion', '> 1bp Insertion', 'Microhomology'))


plt_id <- ggplot(data_join_id, 
                 aes(x = context, y = percentage, fill = causes)) +
  geom_col() +
  ggh4x::facet_nested(spn + sample_id~name+Value, 
                      scales="free", 
                      space="free_x") + 
  ylab("Percentage of ID") +
  my_ggplot_theme() +
  scale_y_continuous(breaks = y_breaks) +  # control y ticks
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1,size=6),
    axis.text.y=element_text(size=5),
    panel.spacing.x = unit(0.05, "in"), 
    panel.spacing.y = unit(0.01, "in"),
    strip.text.x = element_text(size = 8, margin = margin(t = 0, b = 0), colour = 'gray20'),
    strip.text.y = element_text(size = 8, margin = margin(l = 0, r = 0), colour = 'gray20'),
    strip.background = element_rect(colour = "gray70", fill = "gainsboro"),
    #panel.background = element_blank(),      # removes grey background
    panel.grid.major = element_blank(),      # removes major grid lines
    panel.grid.minor = element_blank(),      # removes minor grid lines
    #panel.border = element_rect(color = "black", fill = NA)
  )  + scale_fill_manual(values = id_colors)
plt_id

ggsave(plot = plt_id, filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/id_new.pdf', 
       width = 10, height = 7, units = 'in')

