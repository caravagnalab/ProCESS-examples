plot_scatter_process_single = function(table_wide, s1,s2, color_palette, driver){
  if(driver==F){
    ggplot()+
      geom_point(data=table_wide,aes(x=eval(parse(text = s1)),
                                     y = eval(parse(text = s2)),
                                     color = cluster_id_process),
                 size=0.5) +
      scale_color_manual(values=color_palette)+
      theme_minimal() +
      labs(
        color = "Cluster",
        x = s1,
        y = s2
      )+
      xlim(0, 1)+
      ylim(0, 1)
  }else{
    ggplot() + 
      geom_point(data = table_wide, 
                 aes(x = .data[[s1]], y =.data[[s2]], col = .data$cluster_id_process), 
                 size = .5, alpha = 1) +
      scale_color_manual(values=color_palette)+
      theme_minimal() +
      labs(
        color = "Cluster",
        x = s1,
        y = s2
      )+
      xlim(0, 1)+
      ylim(0, 1)+
      ggtitle('') +
      xlab(s1) +
      ylab(s2) + 
      geom_point(
        data = subset(table_wide, is_driver_process == TRUE),
        aes(x = .data[[s1]], y = .data[[s2]]),
        color = "black", size = 1, shape = 15
      ) +
      ggrepel::geom_label_repel(
        data = subset(table_wide, is_driver_process == TRUE),
        aes(x = .data[[s1]], y = .data[[s2]], label =.data$gene, color = .data$cluster_id_process),
        #color = 'black',
        size = 3,
        # nudge_y = 0,
        # nudge_x = 0,
        show.legend = FALSE,
        max.overlaps = Inf,
        nudge_y = 0.5,
        nudge_x = 0.5
      )  +
      guides(
        color = guide_legend(
          ncol = 1,
          override.aes = list(size = 3, alpha = 1)
        )
      )
  }
}

plot_scatter_tool_single = function(table_wide, s1, s2, color_palette, type){
  if(type =='original'){
    ggplot() + 
      geom_point(data = table_wide, 
                 aes(x = .data[[s1]], y =.data[[s2]], col = .data$cluster_id_tool), 
                 size = .5, alpha = 1) +
      scale_color_manual(values=color_palette)+
      theme_minimal() +
      labs(
        color = "Cluster",
        x = s1,
        y = s2
      )+
      xlim(0, 1)+
      ylim(0, 1)+
      guides(
        color = guide_legend(
          ncol = 1,
          override.aes = list(size = 3, alpha = 1)
        )
      )
  }else if(type=='interpreted'){
    ggplot() + 
      geom_point(data = table_wide, 
                 aes(x = .data[[s1]], y =.data[[s2]], col = .data$cluster_id_tool_interpreted), 
                 size = .5, alpha = 1) +
      scale_color_manual(values=color_palette)+
      theme_minimal() +
      labs(
        color = "Cluster",
        x = s1,
        y = s2
      )+
      xlim(0, 1)+
      ylim(0, 1)+
      ggtitle('') +
      xlab(s1) +
      ylab(s2) + 
      geom_point(
        data = subset(table_wide, is_driver_tool == TRUE),
        aes(x = .data[[s1]], y = .data[[s2]]),
        color = "black", size = 1, shape = 15
      ) +
      ggrepel::geom_label_repel(
        data = subset(table_wide, is_driver_tool == TRUE),
        aes(x = .data[[s1]], y = .data[[s2]], label =.data$gene, color = .data$cluster_id_tool_interpreted),
        #color = 'black',
        size = 3,
        nudge_y = 0.5,
        nudge_x = 0.5,
        show.legend = FALSE,
        max.overlaps = Inf
      )  +
      guides(
        color = guide_legend(
          ncol = 1,
          override.aes = list(size = 3, alpha = 1)
        )
      )
  }else if(type=='interpreted_driver'){
    ggplot() + 
      geom_point(data = table_wide, 
                 aes(x = .data[[s1]], y =.data[[s2]], col = .data$cluster_id_tool_interpreted_driver), 
                 size = .5, alpha = 1) +
      scale_color_manual(values=color_palette)+
      theme_minimal() +
      labs(
        color = "Cluster",
        x = s1,
        y = s2
      )+
      xlim(0, 1)+
      ylim(0, 1)+
      ggtitle('') +
      xlab(s1) +
      ylab(s2) + 
      geom_point(
        data = subset(table_wide, is_driver_process == TRUE),
        aes(x = .data[[s1]], y = .data[[s2]]),
        color = "black", size = 1, shape = 15
      ) +
      ggrepel::geom_label_repel(
        data = subset(table_wide, is_driver_process == TRUE),
        aes(x = .data[[s1]], y = .data[[s2]], label =.data$gene, color = .data$cluster_id_tool_interpreted_driver),
        #color = 'black',
        size = 3,
        nudge_y = 0.5,
        nudge_x = 0.5,
        show.legend = FALSE,
        max.overlaps = Inf
      )  +
      guides(
        color = guide_legend(
          ncol = 1,
          override.aes = list(size = 3, alpha = 1)
        )
      )
    
  }
}

get_plots_path = function(save_path, tool, spn, simulation_id, plot_name) {
  file.path(save_path, "plots", paste0(plot_name, "_", tool, "_", spn, "_", simulation_id, ".png"))
}

get_plots_path_process = function(save_path, tool, spn, simulation_id, plot_name) {
  file.path(save_path, "plots", paste0(plot_name, "_", spn, "_", simulation_id, ".png"))
}

get_plots_path_shared = function(save_path, tool, spn, simulation_id, plot_name) {
  file.path(save_path, "subclonal/plots", paste0(plot_name, "_", tool, "_", spn, "_", simulation_id, ".png"))
}

get_table_path = function(save_path, tool, spn, simulation_id) {
  file.path(save_path, "tables", paste0("table_", tool, "_", spn, "_", simulation_id, ".rds"))
}

plot_scatter_process = function(table, sample_names, color_palette_process, driver=T, vertical = T){
  # table = join_table_process
  if(spn=='SPN07'){
    table = table %>%
      distinct(patient_id, sample_id, mutation_id, .keep_all = TRUE)
  }
  table_wide = table %>%
    mutate(gene=driver_label_process) %>% 
    select(patient_id, sample_id, mutation_id, cluster_id_process, vaf_process, is_driver_process,gene) %>%
    pivot_wider(values_from="vaf_process", names_from="sample_id")
  
  table_wide <- table_wide %>%
    filter(!is.na(cluster_id_process))
  
  # table_wide[is.na(table_wide)] = 0.0
  table_wide <- table_wide %>%
    mutate(across(starts_with("Spn"), ~replace_na(., 0.0)))
  
  sample_names = as.character(sample_names)
  cm = combn(sample_names, 2)
  
  # s1 = cm[[1]]
  # s2 = cm[[2]]
  plots <- apply(
    cm,
    2,
    function(w) plot_scatter_process_single(table_wide, s1 = w[1], s2 = w[2], color_palette=color_palette_process, driver)
  )
  
  if(vertical == F){
    if(cm %>% ncol() == 1){
      nrows = 1
      ncols = 1
      
    }else{
      num_pairs = cm %>% ncol()
      ncols = min(3, num_pairs) # max 3 cols
      nrows = ceiling(num_pairs / ncols)
    }
  }else{
    ncols = 1
    nrows = length(plots)
  }
  
  plot_to_save = ggpubr::ggarrange(
    plotlist = plots,
    ncol = ncols,
    nrow = nrows,
    common.legend = T,
    legend = "right")
  
  # wrap_plots(plots, guides = 'collect')
  
  plot_to_save
}


# plot_scatter_tool(final_table, color_palette=color_palette_tool, sample_names, type ='interpreted', vertical = vertical)
plot_scatter_tool = function(final_table, color_palette, sample_names, 
                             type ='original', vertical = T){
  
  if(type =='original'){
    if(spn=='SPN07'){
      final_table = final_table %>%
        distinct(patient_id, sample_id, mutation_id, .keep_all = TRUE)
    }
    table_wide = final_table %>%
      select(patient_id, sample_id, mutation_id, cluster_id_tool, vaf_tool) %>%
      pivot_wider(values_from="vaf_tool", names_from="sample_id")
    
    table_wide <- table_wide %>%
      filter(!is.na(cluster_id_tool))
    
    table_wide[is.na(table_wide)] = 0.0
    
    table_wide$cluster_id_tool <- factor(
      table_wide$cluster_id_tool,
      levels = sort(unique(table_wide$cluster_id_tool))
    )
    
  }else if(type=='interpreted'){
    if(spn=='SPN07'){
      final_table = final_table %>%
        distinct(patient_id, sample_id, mutation_id, .keep_all = TRUE)
    }
    table_wide = final_table %>%
      select(patient_id, sample_id, mutation_id, cluster_id_tool_interpreted, vaf_tool, is_driver_tool, driver_label_tool) %>%
      mutate(gene=driver_label_tool) %>% 
      pivot_wider(values_from="vaf_tool", names_from="sample_id")
    
    table_wide <- table_wide %>%
      filter(!is.na(cluster_id_tool_interpreted))
    
    table_wide <- table_wide %>%
      mutate(across(starts_with("Spn"), ~replace_na(., 0.0)))
    
    table_wide$cluster_id_tool_interpreted <- factor(
      table_wide$cluster_id_tool_interpreted,
      levels = sort(unique(table_wide$cluster_id_tool_interpreted))
    )
  }else if(type=='interpreted_driver'){
    
    if(spn=='SPN07'){
      final_table = final_table %>%
        distinct(patient_id, sample_id, mutation_id, .keep_all = TRUE)
    }
    table_wide = final_table %>%
      select(patient_id, sample_id, mutation_id, cluster_id_tool_interpreted_driver, 
             vaf_tool, is_driver_tool, is_driver_process, driver_label_tool,driver_label_process) %>%
      mutate(gene=driver_label_process) %>% 
      pivot_wider(values_from="vaf_tool", names_from="sample_id")
    
    table_wide <- table_wide %>%
      filter(!is.na(cluster_id_tool_interpreted_driver))
    
    table_wide <- table_wide %>%
      mutate(across(starts_with("Spn"), ~replace_na(., 0.0)))
    
    table_wide$cluster_id_tool_interpreted_driver <- factor(
      table_wide$cluster_id_tool_interpreted_driver,
      levels = sort(unique(table_wide$cluster_id_tool_interpreted_driver))
    )
  }
  
  sample_names = as.character(sample_names)
  cm = combn(sample_names, 2)
  # s1 = cm[[1]]
  # s2 = cm[[2]]
  plots <- apply(
    cm,
    2,
    function(w) plot_scatter_tool_single(table_wide, s1 = w[1], s2 = w[2], color_palette, type=type)
  )
  
  if(vertical == F){
    if(cm %>% ncol() == 1){
      nrows = 1
      ncols = 1
      
    }else{
      num_pairs = cm %>% ncol()
      ncols = min(3, num_pairs) # max 3 cols
      nrows = ceiling(num_pairs / ncols)
    }
  }else{
    ncols = 1
    nrows = length(plots)
  }
  
  plot_to_save = ggpubr::ggarrange(
    plotlist = plots,
    ncol = ncols,
    nrow = nrows,
    common.legend = T,
    legend = "right")
  
  plot_to_save = plot_to_save + plot_layout(guides = 'collect')
  plot_to_save
}


plot_mutations_on_tree = function(table_tool, 
                                  process_seq, 
                                  sample_forest,
                                  phylo_forest,
                                  color_palette_tool, 
                                  color_palette_process,
                                  plot_tool =TRUE){
  
  # Here I need to do a join between the process table of this simulation and the tool table
  join_table = table_tool %>%
    inner_join(process_seq) %>% 
    select(chr, chr_pos, ref, alt, mutation_id, cluster_id_tool, cluster_id_process, vaf_tool)
  
  mutations_with_cell = join_table %>% 
    rowwise() %>%
    mutate(cell_id=phylo_forest$get_first_occurrences(Mutation(
      chr, chr_pos, ref, alt
    ))[[1]]) %>%
    ungroup()
  
  # tool
  cells_labels_tool = mutations_with_cell %>% 
    select(mutation_id, cell_id, cluster_id_tool) %>% 
    group_by(cell_id) %>% 
    summarise(label_list=list(cluster_id_tool)) %>% 
    rowwise() %>% 
    mutate(label=names(sort(table(label_list[[1]]), decreasing=TRUE))[1]) %>% 
    ungroup() %>% 
    select(-label_list)
  
  final_labels_tool = sample_forest$get_nodes() %>% as_tibble() %>% 
    left_join(cells_labels_tool)
  
  pl_sticks_tool = plot_sticks(sample_forest, labels=final_labels_tool, cls = color_palette_tool) %>%
    annotate_forest(sample_forest, samples=TRUE, drivers=TRUE)
  
  # Process
  cells_labels_process = mutations_with_cell %>% 
    select(mutation_id, cell_id, cluster_id_process) %>% 
    group_by(cell_id) %>% 
    summarise(label_list=list(cluster_id_process)) %>% 
    rowwise() %>% 
    mutate(label=names(sort(table(label_list[[1]]), decreasing=TRUE))[1]) %>% 
    ungroup() %>% 
    select(-label_list)
  
  final_labels_process = sample_forest$get_nodes() %>% as_tibble() %>% 
    left_join(cells_labels_process)
  
  pl_sticks_process = plot_sticks(sample_forest, labels=final_labels_process, cls = color_palette_process) %>%
    annotate_forest(sample_forest, samples=TRUE, drivers=TRUE)

  return(list(pl_sticks_tool,pl_sticks_process))
}

plot_mutations_on_tree_process = function(table_tool, 
                                  process_seq, 
                                  sample_forest,
                                  phylo_forest,
                                  color_palette_process){
  
  # Here I need to do a join between the process table of this simulation and the tool table
  join_table = table_tool %>%
    inner_join(process_seq) %>% 
    select(chr, chr_pos, ref, alt, mutation_id, cluster_id_process)
  
  mutations_with_cell = join_table %>% 
    rowwise() %>%
    mutate(cell_id=phylo_forest$get_first_occurrences(Mutation(
      chr, chr_pos, ref, alt
    ))[[1]]) %>%
    ungroup()
  
  # Process
  cells_labels_process = mutations_with_cell %>% 
    select(mutation_id, cell_id, cluster_id_process) %>% 
    group_by(cell_id) %>% 
    summarise(label_list=list(cluster_id_process)) %>% 
    rowwise() %>% 
    mutate(label=names(sort(table(label_list[[1]]), decreasing=TRUE))[1]) %>% 
    ungroup() %>% 
    select(-label_list)
  
  final_labels_process = sample_forest$get_nodes() %>% as_tibble() %>% 
    left_join(cells_labels_process)
  
  pl_sticks_process = plot_sticks(sample_forest, labels=final_labels_process, cls = color_palette_process) %>%
    annotate_forest(sample_forest, samples=TRUE, drivers=TRUE)
  
  return(pl_sticks_process)
}

get_cell_id = function(mutation_object) {
  tryCatch(
    expr = { phylo_forest$get_first_occurrences(mutation_object)[[1]] },
    error = function(e) return(NA)
  )
}
# table_tool =final_table
# process_seq
# sample_forest
# phylo_forest
# color_palette_process
plot_mutations_on_subtree = function(table_tool, 
                                    process_seq, 
                                    sample_forest,
                                    phylo_forest,
                                    color_palette_tool,
                                    color_palette_process){
  
  # Here I need to do a join between the process table of this simulation and the tool table
  join_table = table_tool %>%
    inner_join(process_seq) %>% 
    select(chr, chr_pos, ref, alt, mutation_id, cluster_id_tool, cluster_id_process, vaf_tool)
  
  mutations_with_cell = join_table %>% 
    rowwise() %>%
    mutate(cell_id=get_cell_id(Mutation(chr, chr_pos, ref, alt))) %>%
    ungroup()
  
  # tool
  cells_labels_tool = mutations_with_cell %>% 
    select(mutation_id, cell_id, cluster_id_tool) %>% 
    group_by(cell_id) %>% 
    summarise(label_list=list(cluster_id_tool)) %>% 
    rowwise() %>% 
    mutate(label=names(sort(table(label_list[[1]]), decreasing=TRUE))[1]) %>% 
    ungroup() %>% 
    select(-label_list)
  
  final_labels_tool = sample_forest$get_nodes() %>% as_tibble() %>%
    left_join(cells_labels_tool)
  
  pl_sticks_tool = plot_sticks(sample_forest, labels=final_labels_tool, cls = color_palette_tool) %>%
    annotate_forest(sample_forest, samples=TRUE, drivers=TRUE)
  
  # Process
  cells_labels_process = mutations_with_cell %>% 
    select(mutation_id, cell_id, cluster_id_process) %>% 
    group_by(cell_id) %>% 
    summarise(label_list=list(cluster_id_process)) %>% 
    rowwise() %>% 
    mutate(label=names(sort(table(label_list[[1]]), decreasing=TRUE))[1]) %>% 
    ungroup() %>% 
    select(-label_list)
  
  final_labels_process = sample_forest$get_nodes() %>% as_tibble() %>% 
    left_join(cells_labels_process)
  
  pl_sticks_process = plot_sticks(sample_forest, labels=final_labels_process, cls = color_palette_process) %>%
    annotate_forest(sample_forest, samples=TRUE, drivers=TRUE)
  
  return(list(pl_sticks_tool,pl_sticks_process))
}
