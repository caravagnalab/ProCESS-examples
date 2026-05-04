library(dplyr)
library(ggplot2)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

info <- tibble(spn = paste0('SPN0',1:7),
               class = c('Other', 'Hypermutant', 'Other', 'Other', 'Hypermutant', 'Other', 'Hypermutant'))

data <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/drivers/performance_driver.rds') %>% 
  left_join(info) %>% 
  mutate(tool = paste(vcf_caller, cna_caller, sep = '-')) 

f1_data <- data %>% 
  mutate(id=paste(spn, tool, class, sep =':')) %>% 
  mutate(comb=paste(coverage,purity,sep=":")) %>% 
  dplyr::select(comb, id, F1) %>%
  tidyr::pivot_wider(names_from = id, values_from = F1) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()


spn_ids <- sapply(strsplit(colnames(f1_data), ":"), `[`, 1)
tools <- sapply(strsplit(colnames(f1_data), ":"), `[`, 2)
class <- sapply(strsplit(colnames(f1_data), ":"), `[`, 3)
coverages <- as.numeric(sapply(strsplit(rownames(f1_data), ":"), `[`, 1))
purities <- as.numeric(sapply(strsplit(rownames(f1_data), ":"), `[`, 2))

column_ha <- HeatmapAnnotation(
  spn = spn_ids, 
  Hypermutant = class,
  col = list(spn=SPN_colors, Hypermutant = c('Other' = 'gray', 'Hypermutant' = 'gray30')),
  annotation_label = c("SPN", 'Hypermutant')
)

row_ha <- rowAnnotation(
  coverage = coverages,
  purity = purities,
  col = list(coverage=coverage_colors, purity=purity_colors),
  show_annotation_name = F
)

col_fun = circlize::colorRamp2(c(min(f1_data,na.rm=TRUE), 
                                 max(f1_data,na.rm=TRUE)), 
                               c("gainsboro", "steelblue"))



f1 <- ComplexHeatmap::Heatmap(f1_data,
                             cluster_rows = F,
                             cluster_columns = F,
                              top_annotation = column_ha,
                              # bottom_annotation = column_bottom_ha,
                              left_annotation = row_ha,
                              col=col_fun,
                              #row_title = "Recall\nbreakpoint",
                              row_title_side = "right",
                              row_title_gp = gpar(fontsize = 12, lineheight = 0.8),
                              row_title_rot = 0,
                              show_column_names = F,
                              show_row_names = F,
                              rect_gp = gpar(col = "white", lwd = 1),
                              column_split = tools,
                              name = "F1"
)


precision_data <- data %>% 
  mutate(id=paste(spn, tool, sep =':')) %>% 
  mutate(comb=paste(coverage,purity,sep=":")) %>% 
  dplyr::select(comb, id, Precision) %>%
  tidyr::pivot_wider(names_from = id, values_from = Precision) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()

col_fun = circlize::colorRamp2(c(min(precision_data,na.rm=TRUE), 
                                 max(precision_data,na.rm=TRUE)), 
                               c("gainsboro", "indianred"))



precision <- ComplexHeatmap::Heatmap(precision_data,
                              cluster_rows = F,
                              cluster_columns = F,
                              #top_annotation = column_ha,
                              # bottom_annotation = column_bottom_ha,
                              left_annotation = row_ha,
                              col=col_fun,
                              #row_title = "Recall\nbreakpoint",
                              row_title_side = "right",
                              row_title_gp = gpar(fontsize = 12, lineheight = 0.8),
                              row_title_rot = 0,
                              show_column_names = F,
                              show_row_names = F,
                              rect_gp = gpar(col = "white", lwd = 1),
                              column_split = tools,
                              name = "Precision"
)

recall_data <- data %>% 
  mutate(id=paste(spn, tool, sep =':')) %>% 
  mutate(comb=paste(coverage,purity,sep=":")) %>% 
  dplyr::select(comb, id, Recall) %>%
  tidyr::pivot_wider(names_from = id, values_from = Recall) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()

col_fun = circlize::colorRamp2(c(0, 
                                 max(recall_data,na.rm=TRUE)), 
                               c("gainsboro", "goldenrod"))



recall <- ComplexHeatmap::Heatmap(recall_data,
                                     cluster_rows = F,
                                     cluster_columns = F,
                                     #top_annotation = column_ha,
                                     # bottom_annotation = column_bottom_ha,
                                     left_annotation = row_ha,
                                     col=col_fun,
                                     #row_title = "Recall\nbreakpoint",
                                     row_title_side = "right",
                                     row_title_gp = gpar(fontsize = 12, lineheight = 0.8),
                                     row_title_rot = 0,
                                     show_column_names = F,
                                     show_row_names = F,
                                     rect_gp = gpar(col = "white", lwd = 1),
                                     column_split = tools,
                                     name = "Recall"
)

final <- f1 %v% precision %v% recall
pdf("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/supp_figures/Final_Driver_SCOUT_Validation.pdf",
    width = 10,height = 8)
draw(
  final,
  #heatmap_legend_side = "bottom",
  #annotation_legend_side = "bottom",
  merge_legends = TRUE   # merges multiple legends into one row if possible
)
dev.off()
