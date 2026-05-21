

# # Loading stuff ####
# forest = load_sample_forest("/orfeo/scratch/cdslab/shared/SCOUT/SPN07/process/sample_forest.sff")
# sub_forest = forest$get_subforest_for(c('SPN07_1.1', 'SPN07_2.1'))
# muts_path = get_mutations("SPN07", coverage = 100, purity = 0.9, type = "tumour")
# 
# treatment_info_df = dplyr::tibble(s=158, e=173)
# 
# tree_plot = my_plot_forest(sub_forest, 
#                            color_map = c(
#                              "Clone 1" = "gainsboro",
#                              "Clone 2" = "gainsboro",
#                              "Clone 3" = "gainsboro",
#                              "Clone 4" = "gainsboro",
#                              "Clone 5" = "gainsboro"
#                            ), color_sample = c('SPN07_1.1' = 'plum4', 'SPN07_2.1' = 'cadetblue')) %>%
#   my_annotate_forest(forest = sub_forest, samples = T, MRCAs = T, drivers = TRUE, 
#                      color_sample = c('    SPN07_1.1' = 'plum4', '    SPN07_2.1' = 'cadetblue')) +
#   geom_rect(
#     data = treatment_info_df,
#     aes(ymin = s, ymax = e, xmin = -Inf, xmax = Inf),
#     inherit.aes = FALSE,
#     alpha = 0.2,
#     fill = "coral"
#   )  + 
#   guides(shape = "none") 
# 
# tree_plot_spn07 = tree_plot + theme(legend.direction = "horizontal")
# 
# # Loading stuff ####
# forest = load_sample_forest("/orfeo/scratch/cdslab/shared/SCOUT/SPN06/process/sample_forest.sff")
# sub_forest = forest$get_subforest_for(c('SPN06_1.1', 'SPN06_2.1', 'SPN06_3.1'))
# treatment_info = readRDS("/orfeo/scratch/cdslab/shared/SCOUT/SPN06/process/chemo_timing.rds")
# treatment_info_df = dplyr::tibble(s=c(treatment_info$chemo1_start, treatment_info$chemo2_start), 
#                                   e=c(treatment_info$chemo1_end, treatment_info$chemo2_end))
# 
# tree_plot = my_plot_forest(sub_forest, 
#                            color_map = c(
#                              "Clone 1" = "gainsboro",
#                              "Clone 2" = "gainsboro",
#                              "Clone 3" = "gainsboro",
#                              "Clone 4" = "gainsboro",
#                              "Clone 5" = "gainsboro"
#                            ), color_sample = c('SPN06_1.1' = 'goldenrod', 'SPN06_2.1' = 'palevioletred', 'SPN06_3.1' = 'cornflowerblue')) %>%
#   my_annotate_forest(forest = sub_forest, samples = T, MRCAs = T, drivers = TRUE, 
#                      color_sample = c('    SPN06_1.1' = 'goldenrod', '    SPN06_2.1' = 'palevioletred', '    SPN06_3.1' = 'cornflowerblue')) +
#   geom_rect(
#     data = treatment_info_df,
#     aes(ymin = s, ymax = e, xmin = -Inf, xmax = Inf),
#     inherit.aes = FALSE,
#     alpha = 0.2,
#     fill = "coral"
#   )  + 
#   guides(shape = "none") 
# 
# tree_plot_spn06 = tree_plot + theme(legend.direction = "horizontal")
# 
# p <- tree_plot_spn04 + ggtitle('SPN04') + 
#   tree_plot_spn06 + ggtitle('SPN06')  +
#   tree_plot_spn07 + ggtitle('SPN07') + 
#   plot_layout(nrow = 1)

# ggsave(filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/mrca/forest_mrca.png',
#        plot = p, width = 11, height = 6, units = 'in', dpi = 600)
