library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ComplexHeatmap)

source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/germline/germline.R')
df_all_SPN_germline_orig_n <- df_all_SPN_germline_orig %>% 
  mutate(spn = ifelse(spn == 'SPN07_last', 'SPN07', spn)) %>% 
  dplyr::select(spn, F1_Score, tool) %>%
  tidyr::pivot_wider(names_from = spn, values_from = F1_Score) %>% 
  tibble::column_to_rownames("tool") %>% 
  as.matrix()


df_mut_germline_n <- df_mut_germline %>% 
  group_by(spn, tool) %>% 
  summarise(mean_RMSE = mean(RMSE),
            mean_Correlation = mean(cor_coeff))  %>% 
  ungroup() %>% 
  mutate(spn = ifelse(spn == 'SPN07_last', 'SPN07', spn))


df_mut_germline_n2 <- df_mut_germline_n %>% 
  dplyr::select(spn, mean_RMSE, tool) %>%
  tidyr::pivot_wider(names_from = spn, values_from = mean_RMSE) %>% 
  tibble::column_to_rownames("tool") %>% 
  as.matrix()

row_ha <- rowAnnotation(
  tool = rownames(df_mut_germline_n2),
  col = list(tool=col_germline_tools),
  show_annotation_name = F
)

column_ha <- HeatmapAnnotation(
  spn = colnames(df_mut_germline_n2), 
  col = list(spn=SPN_colors))

col_fun = circlize::colorRamp2(c(min(df_mut_germline_n2,na.rm=TRUE), 
                                 max(df_mut_germline_n2,na.rm=TRUE)), c("white","skyblue4"))
h_mse <- ComplexHeatmap::Heatmap(df_mut_germline_n2,
                        cluster_rows = F,
                        cluster_columns = F,
                        col=col_fun,
                        top_annotation = column_ha,
                        left_annotation = row_ha,
                        show_column_names = T,
                        row_title = "MSE",
                        row_title_side = "right",
                        row_title_gp = gpar(fontsize = 12, lineheight = 0.8),
                        row_title_rot = 0,
                        show_row_names = F,
                        rect_gp = gpar(col = "white", lwd = 1),
                        name = "MSE"
)


col_fun = circlize::colorRamp2(c(min(df_all_SPN_germline_orig_n,na.rm=TRUE), 
                                 max(df_all_SPN_germline_orig_n,na.rm=TRUE)), c("white","mediumorchid4"))
h_f1 <- ComplexHeatmap::Heatmap(df_all_SPN_germline_orig_n,
                                 cluster_rows = F,
                                 cluster_columns = F,
                                 col=col_fun,
                                 #top_annotation = column_ha,
                                 left_annotation = row_ha,
                                 show_column_names = T,
                                 row_title = "F1 score",
                                 row_title_side = "right",
                                 row_title_gp = gpar(fontsize = 12, lineheight = 0.8),
                                 row_title_rot = 0,
                                 show_row_names = F,
                                 rect_gp = gpar(col = "white", lwd = 1),
                                 name = "F1 score"
)


h_final_germline <- h_mse %v% h_f1

pdf("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/supp_figures/Final_Germline_SCOUT_Validation.pdf",
    width =6,height = 5)
draw(
  h_final_germline,
  #heatmap_legend_side = "bottom",
  #annotation_legend_side = "bottom",
  merge_legends = TRUE   # merges multiple legends into one row if possible
)
dev.off()
