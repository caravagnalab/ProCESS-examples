library(dplyr)
library(tidyverse)
library(GenomicRanges)
library(circlize)
library(ComplexHeatmap)
library(ProCESS)
library(ggrepel)
library(EnrichedHeatmap)


source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/report/plotting/utils.R")

upper <- 0.89
lower <- 0.11

setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/cn_heatmap/')
source('../../getters/process_getters.R')
source('utils.R')

samples_levels <- lapply(c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN05', 'SPN06', 'SPN07'), function(spn){
  samples <- get_sample_names(spn)
}) %>% unlist()
 
cna_data_all <- tibble()
for (spn in c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN05', 'SPN06', 'SPN07')){
  samples <- get_sample_names(spn)
  cna_data <- lapply(samples, FUN = function(s){
    readRDS(get_process_cna(spn = spn, sample = s)) %>% mutate(sample = s)
  }) %>% bind_rows()
  x <- cna_data %>%
    dplyr::mutate(chr=paste0("chr",chr)) %>%
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

out <- segment_fixed_windows_all(cna_data_all, window_size = 1e6) %>%
  mutate(sample_id = factor(sample_id, levels = samples_levels)) %>%
  arrange(sample_id, desc(ratio)) %>%
  mutate(sample_id = paste(sample_id, ratio, sep = ':')) %>%
  select(-segment_id) %>%
  distinct()
# saveRDS(out, file = 'data_heatmap.rds')


#out <- readRDS('data_heatmap.rds') #%>% 
  #separate(sample_id,into = c('spn','tmp'), '_', remove = F) %>% 
  #filter(spn %in% c('SPN02'))

levels_type <- c('clonal', 'sub-clonal')
######## for karotype
wide_df_kar <- out %>%
  mutate(CN=Major+minor) %>% 
  filter(!is.na(ratio)) %>% 
  mutate(type = ifelse(ratio == 1, 'clonal', 'sub-clonal')) %>% 
  tidyr::separate(sample_id, into = c('sample', 'r'), sep = ':') %>% 
  mutate(sample_id = paste(sample, type, sep = ':')) %>% 
  mutate(sample_id = as.character(sample_id)) %>% 
  mutate(type = factor(type, levels = levels_type)) %>% 
  arrange(sample_id, desc(type)) %>%
  select(chr, from, to, sample_id, CN) %>% # Select only relevant columns
  distinct() %>% 
  pivot_wider(names_from = sample_id, values_from = CN,values_fn = max) %>%  #, values_fn = max
  as.data.frame() %>% 
  replace(is.na(.), 'NA')%>% 
  filter(chr %in% paste0("chr", 1:22)) 
  
chr_level = paste0("chr", 1:22)
wide_df_kar <- wide_df_kar %>% mutate(chr = factor(chr, levels = chr_level)) 

chr <- wide_df_kar$chr
subgroup = tibble(sample_id = colnames(wide_df_kar %>% select(-chr, -from, -to)))
subgroup = subgroup %>%
  tidyr::separate(col = sample_id, into = c('sample_id', 'type'), sep = ':', convert = T) %>% 
  tidyr::separate(col = sample_id, into = c('spn', 'tmp'), sep = '_', remove = F) %>% 
  mutate(type = factor(type, levels = levels_type)) %>% 
  select(-tmp) %>% 
  arrange(sample_id, type)


# col_spn <- c("SPN01"='steelblue', "SPN02"='seagreen', "SPN03"='goldenrod2', 
#              "SPN04"='darkorange', "SPN05"="mediumpurple3","SPN06"='palevioletred', "SPN07"='indianred3')
# #col_proportion <- colorRamp2(c(0, 1), c("white", 'dodgerblue4'))
# col_proportion <- c('clonal' = 'navy', 'sub-clonal' = 'lightsteelblue')
# col_sample <-  c('lightblue', 
#                  'steelblue1', 
#                  'steelblue4', 
#                  'darkseagreen', 
#                  'seagreen4', 
#                  'lightgoldenrod1', 
#                  'lightgoldenrod3', 
#                  'goldenrod3', 
#                  'goldenrod4',
#                  'tan1', 
#                  'darkorange3', 
#                  'mediumpurple1',
#                  'mediumpurple3',
#                  'mediumpurple4',
#                  'pink', 
#                  'palevioletred1', 
#                  'palevioletred', 
#                  'violetred3', 
#                  'palevioletred4', 
#                  'indianred1', 
#                  'firebrick2', 
#                  'firebrick',
#                  'indianred3', 
#                  'indianred4')
# names(col_sample) <- samples_levels
# col_sample <- col_sample[!is.na(names(col_sample))]

#assemble final plot
row_ha = rowAnnotation(proportion = subgroup$type,
                       #sample = subgroup$sample_id,
                       #spn = subgroup$spn,
                       col=list(proportion = col_proportion
                                #sample = col_sample,
                                #spn = col_spn
                                )
)

# Create a named vector to assign colors to each chromosome level
colors <- rep(c("gray", "gray30"), length.out = length(levels(chr)))
color_map <- setNames(colors, levels(chr))
chr_colors <- color_map[chr]
column_ha = HeatmapAnnotation(chr = chr, 
                              col=list(chr=chr_colors),
                              # label = anno_mark(at = at, labels = labels),
                              show_legend = F)

my_palette <- c(
  "0"="#3590d2",
  "1"="#7fbbe3",
  "2" = "#d4cfcd",
  "3" = "#fdc38dff",
  "4" = "#fca16cff",
  "5" = "#f67b51ff",
  "6" = "#e7533aff",
  "7" = "#cf2518ff",
  "8" = "#ad0000ff",
  "9" = "#c4198aff",
  "NA" = 'white'
)


tmp <- wide_df_kar %>% select(-chr, -from, -to)
ht <- Heatmap(t(tmp),
        name = "CN",
        col = my_palette,
        column_split = chr,
        top_annotation = column_ha,
        row_split = subgroup$spn,
        right_annotation = row_ha,
        row_gap = unit(1.5, "mm"),  
        cluster_columns = T,
        column_title_gp = gpar(fontsize = 10),
        border = TRUE,
        column_gap = unit(0, "points"), 
        row_title = unique(subgroup$spn), 
        row_labels = sub(".*_", "", subgroup$sample_id), 
        #column_title = ifelse(1:22 %% 2 == 0, paste0("\n", chr_level), paste0(chr_level, "\n")),
        column_title = 1:22, 
        show_row_names = FALSE,
        heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop"), 
        use_raster = F)

pdf('cn_heatmap_spn.pdf', width = 13, height = 7)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = 'right') #annotation_legend_list = sign_legends)
dev.off()

