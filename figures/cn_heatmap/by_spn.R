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


true_ploidy<- tibble(
  sample = c(
    "SPN01_1.1", "SPN01_1.2", "SPN01_1.3",
    "SPN02_1.1", "SPN02_1.2",
    "SPN03_1.1", "SPN03_2.1", "SPN03_3.1", "SPN03_4.1",
    "SPN04_1.1", "SPN04_2.1",
    "SPN05_1.1", "SPN05_1.2", "SPN05_1.3",
    "SPN06_1.1", "SPN06_1.2", "SPN06_2.1", "SPN06_3.1","SPN06_3.2",
    "SPN07_1.1", "SPN07_1.2", "SPN07_1.3", "SPN07_2.1", "SPN07_2.2"
  ),
  true_ploidy = c(
    2.97, 2.06, 4.12,
    1.99, 2.04,
    1.93, 1.92, 1.91, 1.92,
    1.96, 1.96,
    2.02, 2.04, 2.01, 
    1.96, 1.95, 1.94, 1.95, 3,
    2.53, 2.53, 2.53, 2.53, 2.53
  )
)

fga_fgs <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/fga_df.rds')

names(true_ploidy$true_ploidy) = true_ploidy$sample
fga_fgs_ploidy <- fga_fgs %>% 
  left_join(true_ploidy %>% dplyr::rename(spn = sample)) %>% 
  arrange(spn)


ann_data = read.delim('../oncoprint/top_ann_info.csv', sep = ',') %>% 
  bind_rows(c(Sample = 'SPN05_1.1', SPN = 'SPN05', Tumour_type = 'BRCA', Treatment = 'No_treatment', Gender = 'XX', WGD = 'No', Hypermutant = 'No', Clonal_status = 'Monoclonal')) %>% 
  bind_rows(c(Sample = 'SPN05_1.2', SPN = 'SPN05', Tumour_type = 'BRCA', Treatment = 'No_treatment', Gender = 'XX', WGD = 'No', Hypermutant = 'No', Clonal_status = 'Polyclonal')) %>% 
  bind_rows(c(Sample = 'SPN05_1.3', SPN = 'SPN05', Tumour_type = 'BRCA', Treatment = 'No_treatment', Gender = 'XX', WGD = 'No', Hypermutant = 'No', Clonal_status = 'Monoclonal')) %>% 
  arrange(SPN) %>% 
  dplyr::rename(Samples = Clonal_status)

list_ht <- list()
for (spn in c('SPN01', 'SPN02','SPN03','SPN04','SPN05','SPN06','SPN07')) { 
  print(spn)

  samples_levels <-  get_sample_names(spn)
  samples <- get_sample_names(spn)

  ann_data_tmp = ann_data %>% filter(SPN == spn)
  
  cna_data_all <- lapply(samples, FUN = function(s){
    readRDS(get_process_cna(spn = spn, sample = s)) %>% mutate(sample = s)
  }) %>% bind_rows()
  
  x <- cna_data_all %>%
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
  
  out <- segment_fixed_windows_all(x, window_size = 1e6) %>%
    mutate(sample_id = factor(sample_id, levels = samples_levels)) %>%
    arrange(sample_id, desc(ratio)) %>%
    mutate(sample_id = paste(sample_id, ratio, sep = ':')) %>%
    select(-segment_id) %>%
    distinct()
  
  
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
    group_by(chr, from, to, sample_id) %>% 
    summarise(CN = mean(CN), .groups = 'drop') %>% 
    mutate(CN = as.integer(CN)) %>% 
    pivot_wider(names_from = sample_id, values_from = CN) %>%  #, values_fn = max
    as.data.frame() %>% 
    replace(is.na(.), 'NA')%>% 
    filter(chr %in% paste0("chr", 1:22)) 
  
  chr_level = paste0("chr", 1:22)
  wide_df_kar <- wide_df_kar %>% mutate(chr = factor(chr, levels = chr_level)) 
  
  chr <- wide_df_kar$chr
  subgroup = tibble(sample_id = colnames(wide_df_kar %>% select(-chr, -from, -to)))
  #subgroup <- tibble(sample_id = colnames(tmp))
  
  subgroup = subgroup %>%
    tidyr::separate(col = sample_id, into = c('sample_id', 'type'), sep = ':', convert = T) %>% 
    tidyr::separate(col = sample_id, into = c('spn', 'tmp'), sep = '_', remove = F) %>% 
    mutate(type = factor(type, levels = levels_type)) %>% 
    select(-tmp) %>% 
    arrange(sample_id, type)
  
  subgroup = subgroup %>% left_join(fga_fgs %>% dplyr::rename(sample_id = spn)) 
  subgroup = subgroup %>% left_join(ann_data_tmp %>% dplyr::rename(sample_id = Sample))
    
  # Create a named vector to assign colors to each chromosome level
  colors <- rep(c("gray", "gray30"), length.out = length(levels(chr)))
  color_map <- setNames(colors, levels(chr))
  chr_colors <- color_map[chr]
  column_ha = HeatmapAnnotation(chr = chr, 
                                col=list(chr=chr_colors),
                                # label = anno_mark(at = at, labels = labels),
                                show_legend = F,
                                annotation_name_side = "left",
                                annotation_width = unit(3, "mm"))
  
  my_palette <- c(
    "0"="#3590d2",
    "1"="#7fbbe3",
    "2" = "gainsboro", #"#d4cfcd",
    "3" = "#fdc38dff", ##fdc38dff
    "4" = "#fca16cff",
    "5" = "#f67b51ff",
    "6" = "#e7533aff",
    "7" = "#cf2518ff",
    "8" = "#ad0000ff",
    "9" = "#c4198aff",
    "NA" = 'white'
  )
  
  tmp_fga_fgs <- fga_fgs_ploidy %>% filter(spn %in% samples)
  fga_ann = HeatmapAnnotation(FGA = anno_barplot(tmp_fga_fgs$fga, beside = T, border = F, bar_width = 1, gp = gpar(fill = 'snow3', col = 'white'), height = unit(1, "cm")), 
                              which = 'row',
                              #annotation_name_side = 'left', 
                              annotation_name_rot = 0)
  fgs_ann = HeatmapAnnotation(FGS = anno_barplot(tmp_fga_fgs$fgs, beside = T, border = F, bar_width = 1, gp = gpar(fill = 'snow3', col = 'white'), height = unit(1, "cm")), 
                              which = 'row',
                              #annotation_name_side = 'left', 
                              annotation_name_rot = 0)
  ploidy_ann = HeatmapAnnotation(Ploidy = anno_barplot(tmp_fga_fgs$true_ploidy, beside = T, border = F, bar_width = 1, gp = gpar(fill = 'snow3', col = 'white'), height = unit(1, "cm")),
                                 which = 'row',
                                 #annotation_name_side = 'left',
                                 annotation_name_rot = 0)
  #left_ann = c(row_ha, fga_ann, fgs_ann, ploidy_ann, gap = unit(2, "mm")) 

  tmp <- wide_df_kar %>% select(-chr, -from, -to)
  
  col_proportion <- c('clonal' = 'seagreen4', 'sub-clonal' = '#D0F0C0')
  col_fga = c('High FGA' = 'indianred2', "Low FGA"="dodgerblue3")
  col_wgd = c('Yes' =  '#555555','No' = 'gray95')
  levels_type <- c("clonal", "sub-clonal")
  
  subgroup$type <- factor(
    #subgroup$type,
    sub(".*:", "", names(tmp)),
    levels = levels_type
  )
  #assemble final plot
  row_ha = rowAnnotation(`CN profile` = subgroup$type,
                         WGD = subgroup$WGD,
                         `FGA class` = subgroup$fga_class,
                         #spn = subgroup$spn,
                         col=list(`CN profile` = col_proportion,
                                  `FGA class` = col_fga,
                                  WGD = col_wgd
                                  #sample = col_sample,
                                  #spn = col_spn,
                         ), show_annotation_name = F
  )
  
  if (spn == 'SPN01'){
    ht <- Heatmap(t(tmp),
                name = "CN",
                col = my_palette,
                column_split = chr,
                #top_annotation = column_ha,
                row_split = sub(".*_", "", sub(":.*", "", names(tmp))), #subgroup$sample_id,
                right_annotation = row_ha,
                cluster_columns = T,

                row_gap = unit(0.5, "mm"),  # Slightly increased gap for a cleaner look
                column_gap = unit(0.5, "mm"),
                
                column_title_side = "top", # This moves the titles to the bottom
                
                row_title = unique(sub(".*_", "", sub(":.*", "", names(tmp)))), 
                row_labels = sub(".*_", "", sub(":.*", "", names(tmp))), 
                column_title = 1:22, 
                show_row_names = F,
                heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop"), 
            
                
                # --- Aesthetics for a 'Lighter' Theme ---
                border = gpar(col = "grey50", lwd = 0.5),                
                #rect_gp = gpar(col = "white", lwd = 0.5), # White grid lines between cells
                
                column_title_gp = gpar(fontsize = 8, col = "grey30"), # Soften title color
                row_title_gp = gpar(fontsize = 7, col = "grey30"),
                use_raster = F,
                row_title_rot = 0)    # Soften title color
  } else {
    ht <- Heatmap(t(tmp),
                  name = "CN",
                  col = my_palette,
                  column_split = chr,
                  #top_annotation = column_ha,
                  row_split = sub(".*_", "", sub(":.*", "", names(tmp))), #subgroup$sample_id,
                  right_annotation = row_ha,
                  cluster_columns = T,
                  
                  row_gap = unit(0.5, "mm"),  # Slightly increased gap for a cleaner look
                  column_gap = unit(0.3, "mm"),
                  
                  column_title_side = "bottom", # This moves the titles to the bottom
                  column_title = NULL,         # This removes the "1:22" titles
                  show_column_names = FALSE,
                  
                  row_title = unique(sub(".*_", "", sub(":.*", "", names(tmp)))), 
                  row_labels = sub(".*_", "", sub(":.*", "", names(tmp))), 
                  show_row_names = F,
                  heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop"), 
                  
                  # --- Aesthetics for a 'Lighter' Theme ---
                  border = gpar(col = "grey50", lwd = 0.5),
                  #rect_gp = gpar(col = "white", lwd = 0.5), # White grid lines between cells
                  
                  column_title_gp = gpar(fontsize = 8, col = "grey30"), # Soften title color
                  row_title_gp = gpar(fontsize = 7, col = "grey30"),
                  use_raster = F,
                  row_title_rot = 0)    # Soften title color
  }
  list_ht[[spn]] <- ht
}


h_final <- list_ht$SPN01 %v% list_ht$SPN02 %v% list_ht$SPN03 %v% list_ht$SPN04 %v% list_ht$SPN05 %v% list_ht$SPN06 %v% list_ht$SPN07 
pdf('cn_heatmap_all_last.pdf', width = 8.5, height = 5.5)
draw(
  h_final,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE   # merges multiple legends into one row if possible
)
dev.off()

