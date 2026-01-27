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

base = "/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/"


SPN_colors <-c("SPN01"='steelblue', "SPN02"='seagreen', "SPN03"='goldenrod', 
               "SPN04"='coral', "SPN05"="magenta4","SPN06"='palevioletred', "SPN07"='indianred3')


cov = 100
pur = 0.9

cn_caller = 'ascat'
mut_caller = 'mutect2'

df_all_sbs = lapply(names(SPN_colors), FUN = function(spn){
  data_sbs = read.table(get_tumourevo_signatures(spn = spn, coverage = cov, purity = pur, 
                                         tool = 'SigProfiler', context = 'SBS96', 
                                         vcf_caller = mut_caller, 
                                         cna_caller = cn_caller)$context_matrix, header = T, sep = '\t')
  
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
}) %>% bind_rows()

df_all_id = lapply(names(SPN_colors), FUN = function(spn){
  data_id = read.table(get_tumourevo_signatures(spn = spn, coverage = cov, purity = pur, 
                                                 tool = 'SigProfiler', context = 'ID83', 
                                                 vcf_caller = mut_caller, 
                                                 cna_caller = cn_caller)$context_matrix, header = T, sep = '\t')
  
  rownames(data_id) <- data_id$MutationType
  data_id <- data_id %>% select(!MutationType)
  data_id_counts <- t(data_id) %>% as.data.frame()
  rownames(data_id_counts) <- sub("^[^_]+_(.*)", "\\1", rownames(data_id_counts))
  data_id_counts = data_id_counts[rowSums(data_id_counts) > 0,]
  
  data <- data_id_counts %>%
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
}) %>% bind_rows()



df_plot_sbs <- df_all_sbs %>% 
  bascule:::reformat_contexts(what = "SBS") %>%
  separate(variant, into = c("v1", "v2"), sep = ">", remove = FALSE) %>% 
  separate(context, into = c("c1", "c2"), sep = "_", remove = FALSE) %>% 
  mutate(context = paste0(c1, v1, c2)) %>%  
  separate(samples, c("spn", "sample_id"), sep = "_")

y_breaks <- function(x) {
  x <- x[x != 0]                # exclude 0
  if(length(x) == 0) return(NULL)
  max_val <- max(x, na.rm = TRUE) - 0.5
  half_val <- max_val / 2
  max_val <- round(max_val, 0)
  half_val <- round(half_val, 0)
  return(c(half_val, max_val))
}

plt_sbs <- ggplot(df_plot_sbs, 
                 aes(x = context, y = percentage)) +
  geom_col() +
  ggh4x::facet_nested(spn + sample_id~variant, 
                      scales="free", 
                      space="free_x") + 
  ylab("Percentage of SBS") +
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
  )

ggsave(plot = plt_sbs, filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/sbs.pdf', 
       width = 10, height = 5, units = 'in')

df_plot_id <- df_all_id %>% 
  bascule:::reformat_contexts(what = "ID") %>%
  separate(col = variant, into = c('N', 'Type', 'Value')) %>% 
  mutate(name = paste(N, Type, sep = ' ')) %>% 
  separate(samples, c("spn", "sample_id"), sep = "_") %>% 
  mutate(name = ifelse(Value == 'R' & Type == 'Del', '> 1bp Deletion', name)) %>%
  mutate(name = ifelse(Value == 'R' & Type == 'Ins', '> 1bp Insertion', name)) %>% 
  mutate(Value = ifelse(Value == 'R', N, Value)) %>% 
  mutate(name = ifelse(Value == 'M', 'Microhomology', name)) %>% 
  mutate(Value = ifelse(Value == 'M', N, Value)) %>% 
  mutate(name = case_when(
    name == '1 Del' ~ '1bp Deletion',
    name == '1 Ins' ~ '1bp Insertion',.default = name
  ))

df_plot_id$name <- factor(df_plot_id$name, levels = c('1bp Deletion', '1bp Insertion', '> 1bp Deletion', '> 1bp Insertion', 'Microhomology'))


plt_id <- ggplot(df_plot_id, 
                  aes(x = context, y = percentage)) +
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
  )
  
ggsave(plot = plt_id, filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/id.pdf', 
       width = 10, height = 5, units = 'in')
