library(dplyr)
library(ggplot2)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/multivariate.R')
plt_subclone <- results %>% 
  pivot_longer(cols = c(viber, pyclone)) %>% 
  ggplot() +
  geom_boxplot(aes(x = as.factor(purity), y = value, col = type, fill = type), alpha =.2) +
  my_ggplot_theme() + 
  xlab('Purity') +
  scale_fill_manual('SPN', values = c('darkseagreen4', 'orangered3')) + 
  scale_color_manual('SPN', values = c('darkseagreen4', 'orangered3')) + 
  ylab('Relative number of inferred cluster\nClusters / Samples') + 
  ggh4x::facet_nested( vcf_caller + cna_caller ~ name + coverage)

ggsave("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/supp_figures/Final_SD_Multivariate_SCOUT_Validation.pdf",
       width = 9, height = 5)


mono_data <- results %>% 
  pivot_longer(cols = c('viber', 'pyclone')) %>% 
  mutate(tool = paste(vcf_caller, cna_caller, sep = '-')) %>% 
  mutate(id=paste(spn, tool, type,  sep =':')) %>% 
  mutate(comb=paste(coverage, purity, name, sep=":")) %>% 
  dplyr::select(comb, id, value) %>%
  tidyr::pivot_wider(names_from = id, values_from = value) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()

spn_ids <- sapply(strsplit(colnames(mono_data), ":"), `[`, 1)
comb <- sapply(strsplit(colnames(mono_data), ":"), `[`, 2)
class <- sapply(strsplit(colnames(mono_data), ":"), `[`, 3)
coverages <- as.numeric(sapply(strsplit(rownames(mono_data), ":"), `[`, 1))
purities <- as.numeric(sapply(strsplit(rownames(mono_data), ":"), `[`, 2))
tools <- sapply(strsplit(rownames(mono_data), ":"), `[`, 3)


column_ha <- HeatmapAnnotation(
  spn = spn_ids, 
  class = class,
  col = list(spn=SPN_colors, class = c("No Bias" = 'darkseagreen4', "Sampling Bias" = 'orangered3')),
  annotation_label = c("SPN", 'Class')
)

row_ha <- rowAnnotation(
  coverage = coverages,
  purity = purities,
  tool = tools, 
  col = list(coverage=coverage_colors, purity=purity_colors, tool = c("pyclone" = "chocolate", "viber" = "darkgoldenrod2")),
  show_annotation_name = F
)

col_fun = circlize::colorRamp2(c(min(mono_data,na.rm=TRUE), 
                                 max(mono_data,na.rm=TRUE)), 
                               c("gainsboro", "dodgerblue4"))
ht <- ComplexHeatmap::Heatmap(mono_data,
                              cluster_rows = F,
                              cluster_columns = F,
                              top_annotation = column_ha,
                              # bottom_annotation = column_bottom_ha,
                              left_annotation = row_ha,
                              col=col_fun,
                              #row_title = "Recall\nbreakpoint",
                              row_title_side = "right",
                              row_title_gp = gpar(fontsize = 0, lineheight = 0.8),
                              row_title_rot = 0,
                              show_column_names = F,
                              show_row_names = F,
                              rect_gp = gpar(col = "white", lwd = 1),
                              column_split = comb,
                              row_split = tools, 
                              name = "Clusters/Samples"
)


pdf("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/supp_figures/Final_SD_Multivariate_SCOUT_Validation_2.pdf",
    width = 10,height = 4)
draw(
  ht,
  #heatmap_legend_side = "bottom",
  #annotation_legend_side = "bottom",
  merge_legends = TRUE   # merges multiple legends into one row if possible
)
dev.off()
  
