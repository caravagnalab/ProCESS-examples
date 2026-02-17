library(dplyr)
library(tidyverse)
library(GenomicRanges)
library(circlize)
library(ComplexHeatmap)
library(ProCESS)
library(EnrichedHeatmap)

source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/report/plotting/utils.R")

upper <- 0.89
lower <- 0.11

col_fga = c('High FGA' = 'indianred2', "Low FGA"="dodgerblue3")
col_wgd = c('Yes' =  '#555555','No' = 'gray95')
col_proportion <- colorRamp2(c(0, 1), c("white", 'goldenrod'))

my_palette <- c(
  "1:0"="#3590d2",
  "1:1" = "#d4cfcd",
  "2:0"=  "#7fbbe3",
  "2:1" = "#FDC38D",
  "2:2" = "#F7B583",
  "3:0" = "#F2A779",
  "3:1" = "#EC996F",
  "3:2" = "#E68B65",
  "3:3" = "#E07D5B",
  "4:0" = "#DB6F51",
  "4:1" = '#D56247',
  "4:2" = '#CF543C',
  "4:4" = '#CA4632',
  "5:0" = '#C43828',
  "5:1" = '#BE2A1E',
  "6:0" = '#B81C14',
  "6:2" = '#B30E0A',
  "8:0" = '#AD0000'
)


setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/cn_heatmap/')
source('../../getters/process_getters.R')
source('utils.R')

samples_levels <- lapply(c('SPN01', 'SPN02', 'SPN03', 'SPN04','SPN05', 'SPN06', 'SPN07'), function(spn){
  samples <- get_sample_names(spn)
}) %>% unlist()

fga_fgs <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/fga_df.rds')
ann_data = read.delim('../oncoprint/top_ann_info.csv', sep = ',') %>% 
  bind_rows(c(Sample = 'SPN05_1.1', SPN = 'SPN05', Tumour_type = 'BRCA', Treatment = 'No_treatment', Gender = 'XX', WGD = 'No', Hypermutant = 'No', Clonal_status = 'Monoclonal')) %>% 
  bind_rows(c(Sample = 'SPN05_1.2', SPN = 'SPN05', Tumour_type = 'BRCA', Treatment = 'No_treatment', Gender = 'XX', WGD = 'No', Hypermutant = 'No', Clonal_status = 'Polyclonal')) %>% 
  bind_rows(c(Sample = 'SPN05_1.3', SPN = 'SPN05', Tumour_type = 'BRCA', Treatment = 'No_treatment', Gender = 'XX', WGD = 'No', Hypermutant = 'No', Clonal_status = 'Monoclonal')) %>% 
  arrange(SPN) %>% 
  dplyr::rename(Samples = Clonal_status)

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


out <- segment_fixed_windows_all(cna_data_all, window_size = 1e6) %>% 
  mutate(sample_id = factor(sample_id, levels = samples_levels)) %>% 
  arrange(sample_id, desc(ratio)) %>% 
  mutate(sample_id = paste(sample_id, ratio, sep = ':')) %>% 
  tidyr::separate(sample_id, into = c('spn','tmp'), '_', remove = F)

df_kar <- out %>%
  mutate(
    CN    = paste0(Major, ":", minor),
    totCN = Major + minor,
    sample_id = as.character(sample_id)
  ) %>%
  na.omit() %>%
  group_by(chr, from, to, sample_id) %>%
  slice_max(
    order_by = totCN,
    n = 1,
    with_ties = FALSE   # guarantees one row
  ) %>%
  ungroup() %>%
  select(chr, from, to, sample_id, CN) %>% 
  unique() %>% 
  mutate(sid = sub(":.*", "", sample_id)) %>% 
  mutate(spn = sub("_.*", "", sid)) 

wide_df_kar <- df_kar %>% 
  mutate(ratio = as.numeric(sub(".*:", "", sample_id))) %>% 
  arrange(spn, sid, desc(ratio)) %>% 
  select(-sid, -spn, -ratio) %>% 
  pivot_wider(names_from = sample_id, values_from = CN) %>% 
  as.data.frame() %>% 
  replace(is.na(.), 'NA') %>% 
  filter(chr %in%  1:22) %>% 
  mutate(chr = paste0('chr', chr))

chr_level = paste0("chr", 1:22)
wide_df_kar <- wide_df_kar %>% mutate(chr = factor(chr, levels = chr_level)) #%>% arrange(chr)
chr <- wide_df_kar$chr


subgroup = tibble(sample_id = colnames(wide_df_kar %>% select(-chr, -from, -to))) %>% 
  mutate(sid = sub(":.*", "", sample_id)) %>% 
  mutate(spn = sub("_.*", "", sid)) %>% 
  left_join(fga_fgs %>% dplyr::rename(sid = spn)) %>% 
  left_join(ann_data %>% dplyr::rename(sid = Sample)) %>% 
  tidyr::separate(col = sample_id, into = c('tmp', 'ratio'), sep = ':', convert = T) %>% 
  select(-tmp) #%>% 
  #arrange(spn, sid, desc(ratio))

row_ha = rowAnnotation(
  proportion = subgroup$ratio,
  WGD = subgroup$WGD,
  `FGA class` = subgroup$fga_class,
  col = list(proportion = col_proportion, WGD = col_wgd, `FGA class` = col_fga),
  show_annotation_name = FALSE,
  simple_anno_size = unit(1.5, "mm"),
  
  # This is where you shrink the text
  annotation_legend_param = list(
    proportion = list(
      title_gp = gpar(fontsize = 3), # Title size
      labels_gp = gpar(fontsize = 3),
      grid_height = unit(2, "mm"),
      grid_width  = unit(2.5, "mm")
    ),
    WGD = list(
      title_gp = gpar(fontsize = 3),
      labels_gp = gpar(fontsize = 3),
      grid_height = unit(2, "mm"),
      grid_width  = unit(2.5, "mm")
    ),
    `FGA class` = list(
      title_gp = gpar(fontsize = 3),
      labels_gp = gpar(fontsize = 3),
      grid_height = unit(2, "mm"),
      grid_width  = unit(2.5, "mm")
    )
  )
)


cn_names <- unique(out %>%
                     mutate(CN=paste0(Major, ':',minor)) %>% pull(CN))

my_palette_cn = my_palette[which(names(my_palette) %in% cn_names)]
my_palette_cn = c(my_palette_cn, 'NA' = 'white')

tmp <- wide_df_kar %>% select(-chr, -from, -to)
ht <- Heatmap(t(tmp),
              name = "CN",
              col = my_palette_cn,
              column_split = chr,
              row_split = subgroup$spn,
              right_annotation = row_ha,
              cluster_columns = T,
              column_title_side = "bottom",
              #row_title = sub(".*_", "", sub(":.*", "", names(tmp))), 
              row_labels = rep('', nrow(t(tmp))),
              column_title = 1:22,
              heatmap_legend_param = list(direction = "horizontal", 
                                          title_position = "lefttop",
                                          title_gp  = gpar(fontsize = 3),
                                          labels_gp = gpar(fontsize = 3), 
                                          grid_height = unit(2, "mm"),
                                          grid_width  = unit(2.5, "mm")), 
              use_raster = F,
              show_column_names = F,
              
              row_gap = unit(0.9, "mm"),  # Slightly increased gap for a cleaner look
              column_gap = unit(0.3, "mm"),
              
              # --- Aesthetics for a 'Lighter' Theme ---
              border = gpar(col = "grey30", lwd = 1),                
              
              column_title_gp = gpar(fontsize = 3, col = "grey30"), # Soften title color
              row_title_gp = gpar(fontsize = 3, col = "grey30"),
              row_title_rot = 0,
              show_heatmap_legend = T)

pdf('cn_heatmap_last.pdf', width = 5, height = 1.5)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
)
dev.off()
