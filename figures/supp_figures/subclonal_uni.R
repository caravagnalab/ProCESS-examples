library(dplyr)
library(ggplot2)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/subclone_mobster.R')
plt_subclone <- results %>% 
  filter(!is.na(gt)) %>% 
  mutate(TP = gt & pred,
         FP = !gt & pred,
         FN = gt & !pred,
         TN = !gt & !pred) %>% 
  group_by(subclass, purity, coverage, vcf_caller, cna_caller) %>% 
  summarise(Polyclonal = sum(TP) / (sum(TP) + sum(FN)),  #recall
            Monoclonal = sum(TN) / (sum(TN) + sum(FP))) %>%  #Specificity
  pivot_longer(cols = c(Polyclonal, Monoclonal),
               names_to = "metric",
               values_to = "value") %>%
  ggplot(aes(x = as.factor(purity), 
             y = value, 
             color = subclass, 
             group = subclass)) +
  geom_line() +
  geom_point(size = 2) +
  scale_color_manual('Subclass', values = c('Polyclonal' = 'goldenrod', 
                                            'Monoclonal'='#645394', 
                                            'Weak evidence' = 'khaki3',
                                            'Sampling Bias' = 'plum')) +
  ylab("Sample composition detection\nTPR") +
  xlab("Purity") +
  my_ggplot_theme() +
  ylim(0,1) + 
  ggh4x::facet_nested(vcf_caller + cna_caller~metric + coverage )

ggsave("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/supp_figures/Final_SD_Univariate_SCOUT_Validation.pdf",
    width = 9, height = 5)


data <- results %>% 
  filter(!is.na(gt)) %>% 
  mutate(TP = gt & pred,
         FP = !gt & pred,
         FN = gt & !pred,
         TN = !gt & !pred) %>% 
  group_by(sample_id, subclass, purity, coverage, vcf_caller, cna_caller) %>% 
  summarise(Polyclonal = sum(TP) / (sum(TP) + sum(FN)),  #recall
            Monoclonal = sum(TN) / (sum(TN) + sum(FP))) %>% 
  ungroup()

mono_data <- data %>% 
  rowwise() %>% 
  mutate(detection = sum(Polyclonal, Monoclonal, na.rm = T)) %>% 
  mutate(tool = paste(vcf_caller, cna_caller, sep = '-')) %>% 
  mutate(id=paste(sample_id, tool, subclass, sep =':')) %>% 
  mutate(comb=paste(coverage,purity,sep=":")) %>% 
  dplyr::select(comb, id, detection) %>%
  tidyr::pivot_wider(names_from = id, values_from = detection) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()

samples_id <- sapply(strsplit(colnames(mono_data), ":"), `[`, 1)
spn_ids <- sub("_.*", "", samples_id)
tools <- sapply(strsplit(colnames(mono_data), ":"), `[`, 2)
class <- sapply(strsplit(colnames(mono_data), ":"), `[`, 3)
coverages <- as.numeric(sapply(strsplit(rownames(mono_data), ":"), `[`, 1))
purities <- as.numeric(sapply(strsplit(rownames(mono_data), ":"), `[`, 2))

class_levels <- c('Monoclonal',  'Sampling Bias', 'Polyclonal', 'Weak evidence')
col_order <- order(factor(class, levels = class_levels))

# Reorder matrix columns and annotation vectors
mono_data <- mono_data[, col_order]
spn_ids   <- spn_ids[col_order]
tools     <- tools[col_order]
class     <- class[col_order]

column_ha <- HeatmapAnnotation(
  spn = spn_ids, 
  class = class,
  col = list(spn=SPN_colors, class = c('Polyclonal' = 'goldenrod', 
                                       'Monoclonal'='#645394', 
                                       'Weak evidence' = 'khaki3',
                                       'Sampling Bias' = 'plum')),
  annotation_label = c("SPN", 'Class')
)

row_ha <- rowAnnotation(
  coverage = coverages,
  purity = purities,
  col = list(coverage=coverage_colors, purity=purity_colors),
  show_annotation_name = F
)


ht <- ComplexHeatmap::Heatmap(mono_data,
                              cluster_rows = F,
                              cluster_columns = F,
                              top_annotation = column_ha,
                              # bottom_annotation = column_bottom_ha,
                              left_annotation = row_ha,
                              col=c("1" = 'cadetblue4', '0' = 'seashell2'),
                              #row_title = "Recall\nbreakpoint",
                              row_title_side = "right",
                              row_title_gp = gpar(fontsize = 12, lineheight = 0.8),
                              row_title_rot = 0,
                              show_column_names = F,
                              show_row_names = F,
                              rect_gp = gpar(col = "white", lwd = 1),
                              column_split = tools,
                              name = "Detected"
)


pdf("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/supp_figures/Final_SD_Univariate_SCOUT_Validation_2.pdf",
    width = 10,height = 4)
draw(
  ht,
  #heatmap_legend_side = "bottom",
  #annotation_legend_side = "bottom",
  merge_legends = TRUE   # merges multiple legends into one row if possible
)
dev.off()
