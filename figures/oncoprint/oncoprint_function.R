# preparing oncoprint functions (as much generalised as possible)

create_annotation = function(x, ann_colors, position, show_legend = T,ann_name = TRUE, pos, annotation_legend_param = list(nrow = 3, width = 12, by_row = T)) {
  
  ComplexHeatmap::HeatmapAnnotation(df = x, 
                                    col = ann_colors, 
                                    which = position,
                                    show_annotation_name = ann_name, 
                                    na_col = '#EEEEEE',
                                    annotation_name_side = pos, 
                                    show_legend = show_legend,
                                    annotation_legend_param = annotation_legend_param)
  
}



create_cohort_oncoprint = function(data, 
                                   top_ann_list, 
                                   bottom_ann_list, 
                                   right_ann_list, 
                                   left_ann_list, 
                                   colors_and_pch, 
                                   ht_title, 
                                   ncol_ann,
                                   column_tt, 
                                   row_tt,
                                   split_columns = NULL, 
                                   split_rows = NULL, 
                                   ht_name, 
                                   ht_path) {
  
  # define the alteration function
  vars_fun = lapply(colors_and_pch$ht_colors, function(x) {
    alter_graphic("rect", width = 0.8, height = 0.8, fill = unname(x))
  })
  vars_fun = c(vars_fun, 'background' = alter_graphic("rect", width = 0.8, height = 0.8, fill = "#CCCCCC"))
  # test_alter_fun(vars_fun)
  
  cn_fun = lapply(colors_and_pch$CNA, function(s) {
    function(x, y, w, h) 
      grid.points(x,y,w,h, pch = unname(s), size = unit(3, "mm"))
  })
  
  cn_fun = c(cn_fun, "NA" = function(x, y, w, h) {
    grid.points(x,y,w,h, pch = "", size = unit(3, "mm"), gpar(col = NA, alpha = 0, fill = NA))
  })
  
  
  alter_fun = c(vars_fun, cn_fun)
  # test_alter_fun(alter_fun)
  
  
  # annotations
  
  heatmap_legend_param = list(title = ht_title, ncol = ncol_ann)
  
  ht = ComplexHeatmap::oncoPrint(data,
                                 col = colors_and_pch$ht_colors, 
                                 alter_fun = alter_fun, 
                                 show_row_names = T, 
                                 show_column_names = T, 
                                 # column_title = column_tt,
                                 heatmap_legend_param = heatmap_legend_param,
                                 # show_pct = F, 
                                 # column_names_gp = gpar(fontsize = 7),
                                 row_names_gp = gpar(fontsize = 8), 
                                 name = ht_name, 
                                 bottom_annotation = bottom_ann_list, 
                                 # top_annotation = top_ann_list,
                                 # right_annotation = right_ann_list, 
                                 # left_annotation = left_ann_list,
                                 remove_empty_rows = F, 
                                 column_split = split_columns, 
                                 row_split = split_rows
  )  
  # pdf(file = ht_path, width = 20, height = 15, bg = "white")
  draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  # graphics.off()
}

