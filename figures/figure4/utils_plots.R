my_theme = theme_light(base_size=8) +
  theme(legend.position="bottom",
        legend.key.height=unit(0.1, "cm"),
        legend.key.width=unit(0.1, "cm"),
        legend.key.spacing=unit(0.1, "cm"),
        panel.background=element_rect(fill="white"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title=element_text(size=6),
        legend.text=element_text(size=6, margin=margin(l=0.1, unit="cm")),
        legend.title=element_text(size=6),
        text=element_text(size=8,family="Helvetica"),
        panel.grid = element_blank())

### Scatter ####

plot_scatter = function(table_wide, s1,s2){
  ggplot(table_wide, aes(x = eval(parse(text = s1)),
                         y = eval(parse(text = s2)),
                         color = cluster_id)) +
    geom_point(size = 0.1) +
    theme_minimal() +
    labs(
      x = s1,
      y = s2
    )+
    xlim(0,1)+
    ylim(0,1)
  
}

plot_scatter_driver = function(table_wide, s1, s2, color_palette){
  ggplot() + 
    geom_point(data = table_wide, 
               aes(x = .data[[s1]], y =.data[[s2]], col = .data$cluster_id), 
               size = .5,
               alpha = 1) +
    scale_color_manual(values=color_palette)+
    theme_minimal()+
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
      data = subset(table_wide, is_driver == TRUE),
      aes(x = .data[[s1]], y = .data[[s2]]),
      color = "black", size = 1, shape = 15
    ) +
    ggrepel::geom_label_repel(
      data = subset(table_wide, is_driver == TRUE),
      aes(x = .data[[s1]], y = .data[[s2]], label =.data$gene, color = .data$cluster_id),
      #color = 'black',
      size = 4,
      nudge_y = 0,
      nudge_x = 0,
      show.legend = FALSE,
      max.overlaps = Inf
    )  +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1)
    )+
    guides(color = guide_legend(override.aes = list(size = 3, alpha=1)))+
    theme(legend.text=element_text(size=10),
          legend.title=element_text(size=10))
}

# Devo usare i dati e i colori di raw.
# Create a column is_merged (TRUE: no driver; FALSE, driver)
# questa colonna la uso per il cerchio attorno
# Create a column is_true_driver (TRUE: alpha 1, true driver; FALSE, alpha 0.5, fake driver)
# i true saranno tutti tranne driver clone 4.
# s1 = 'Sample_1'
# s2 = 'Sample_2'
# color_palette = color_palette_tool

plot_scatter_tool_final = function(table_wide, s1, s2, color_palette){
  # n_sub <- 2000
  # 
  # table_sub <- table_wide %>%
  #   group_by(cluster_id) %>%
  #   slice_sample(n = round(n_sub * n()/nrow(table_wide))) %>%
  #   ungroup()
  # 
  # table_sub <- table_wide %>%
  #   group_by(cluster_id) %>%
  #   slice_sample(prop = 2000 / nrow(table_wide)) %>%
  #   ungroup()
  
  ggplot(table_wide) + 
    geom_point(
      aes(
        x = .data[[s1]],
        y = .data[[s2]],
        color = .data$cluster_id,       # inside of dot
        alpha = ifelse(.data$is_true_driver, 1, 0.1)  # conditional transparency
      ),
      shape = 16,  # filled circle with outline
      size = 1.5,
      stroke = 0.7
    )+
    scale_color_manual(values=color_palette)+
    theme_minimal()+
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
    scale_alpha_identity() +
    geom_point(
      data = subset(table_wide, is_driver == TRUE),
      aes(x = .data[[s1]], y = .data[[s2]]),
      color = "black", size = 1, shape = 15
    ) +
    ggrepel::geom_label_repel(
      data = subset(table_wide, is_driver == TRUE),
      aes(x = .data[[s1]], y = .data[[s2]], label =.data$gene, color = .data$cluster_id),
      #color = 'black',
      size = 3,
      nudge_y = 0,
      nudge_x = 0,
      show.legend = FALSE,
      max.overlaps = Inf
    )  #+
    guides(color = guide_legend(override.aes = list(size = 3, alpha=1)))
}

plot_scatter_tool_all_infos = function(table_wide, s1, s2, color_palette){
  table_wide$merged_status <- ifelse(table_wide$is_merged, "Without driver", "With driver")
  
  ggplot(table_wide) + 
    geom_point(
      aes(
        x = .data[[s1]],
        y = .data[[s2]],
        color = .data$cluster_id,
        shape = .data$merged_status, 
        size = .data$merged_status, 
        # alpha = ifelse(.data$is_true_driver, 1, 0.1)  # conditional transparency
      ),
      alpha=1
    ) +
    scale_color_manual(values = color_palette) +
    # scale_shape_manual(values = c("Without driver" = 16, "With driver" = 17)) +
    scale_shape_manual(values = c("Without driver" = 17, "With driver" = 16)) +
    # scale_shape_manual(values = c("Without driver" = 16, "With driver" = 1)) +
    scale_size_manual(values = c("Without driver" = 1, "With driver" = 1)) +
    theme_minimal() +
    labs(
      color = "Cluster",
      shape = "Interpreted",
      x = s1,
      y = s2
    ) +
    xlim(0, 1)+
    ylim(0, 1)+
    ggtitle('') +
    xlab(s1) +
    ylab(s2) + 
    geom_point(
      data = subset(table_wide, is_driver == TRUE),
      aes(x = .data[[s1]], y = .data[[s2]]),
      color = "black", size = 1, shape = 15
    ) +
    ggrepel::geom_label_repel(
      data = subset(table_wide, is_true_driver == TRUE),
      aes(x = .data[[s1]], y = .data[[s2]], label =.data$gene, color = .data$cluster_id),
      #color = 'black',
      size = 4,
      nudge_y = 0,
      nudge_x = 0,
      fill = "white", 
      show.legend = FALSE,
      max.overlaps = Inf,
    )  +
    ggrepel::geom_label_repel(
      data = subset(table_wide, is_true_driver == FALSE),
      aes(x = .data[[s1]], y = .data[[s2]], label =.data$gene, color = .data$cluster_id),
      size = 4,
      nudge_y = 0,
      nudge_x = 0,
      show.legend = FALSE,
      fill = "gainsboro",
      alpha = 1,
      linetype = "dashed",
      max.overlaps = Inf
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha=1), nrow=1),
           shape = guide_legend(override.aes = list(size = 3, alpha=1), nrow=1),
           size = "none")+
    theme(legend.text=element_text(size=10),
          legend.title=element_text(size=10))
  
}

plot_scatter_tool_final_old = function(table_wide, s1, s2, color_palette){
  table_wide$merged_status <- ifelse(table_wide$is_merged, "No driver cluster", "Driver cluster")
  
  ggplot(table_wide) + 
    geom_point(
      aes(
        x = .data[[s1]],
        y = .data[[s2]],
        fill = .data$cluster_id,       # inside of dot
        color = .data$merged_status,   # outline mapped to legend
        alpha = ifelse(.data$is_true_driver, 1, 0.1)  # conditional transparency
      ),
      shape = 21,  # filled circle with outline
      size = 1.5,
      stroke = 0.7
    ) +
    scale_fill_manual(values = color_palette) +
    scale_color_manual(values = c("No driver cluster" = "grey", "Driver cluster" = "black")) +
    scale_alpha_identity() +   # use alpha values literally
    theme_minimal() +
    labs(
      color = "True/False clone",
      fill = "Cluster",
      x = s1,
      y = s2
    ) +
    xlim(0, 1)+
    ylim(0, 1)+
    ggtitle('') +
    xlab(s1) +
    ylab(s2) + 
    geom_point(
      data = subset(table_wide, is_driver == TRUE),
      aes(x = .data[[s1]], y = .data[[s2]]),
      color = "black", size = 1, shape = 15
    ) +
    ggrepel::geom_label_repel(
      data = subset(table_wide, is_driver == TRUE),
      aes(x = .data[[s1]], y = .data[[s2]], label =.data$gene, color = .data$cluster_id),
      #color = 'black',
      size = 3,
      nudge_y = 0,
      nudge_x = 0,
      show.legend = FALSE,
      max.overlaps = Inf
    )  +
    guides(color = guide_legend(override.aes = list(size = 3, alpha=1)))
  
  # + 
  
  # small black driver points
  # geom_point(
  #   data = subset(table_wide, is_driver == TRUE),
  #   aes(x = .data[[s1]], y = .data[[s2]]),
  #   color = "black", size = 1, shape = 15
  # ) +
  # 
  # # ggrepel labels colored like fill
  # ggrepel::geom_label_repel(
  #   data = subset(table_wide, is_driver == TRUE),
  #   aes(
  #     x = .data[[s1]], 
  #     y = .data[[s2]], 
  #     label = .data$gene,
  #     fill = .data$cluster_id       # set fill to match dot
  #   ),
  #   color = "white",       # text color inside label
  #   size = 3,
  #   label.size = 0.2,      # thin border for the label
  #   fontface = "bold",
  #   show.legend = FALSE,
  #   max.overlaps = Inf
  # ) +
  # 
  # guides(
  #   color = guide_legend(override.aes = list(size = 3, alpha = 1)),
  #   fill = guide_legend(override.aes = list(shape = 21, size = 3))
  # )
  
  # geom_point(
  #   data = subset(table_wide, is_driver == TRUE),
  #   aes(x = .data[[s1]], y = .data[[s2]]),
  #   color = "black", size = 1, shape = 15
  # ) +
  # ggrepel::geom_label_repel(
  #   data = subset(table_wide, is_driver == TRUE),
  #   aes(x = .data[[s1]], y = .data[[s2]], label =.data$gene, color = .data$cluster_id),
  #   #color = 'black',
  #   size = 3,
  #   nudge_y = 0,
  #   nudge_x = 0,
  #   show.legend = FALSE,
  #   max.overlaps = Inf
  # )  +
  # guides(color = guide_legend(override.aes = list(size = 3, alpha=1)))
}
plot_scatter_driver_cartoon = function(table_wide, s1, s2, color_palette, psize = .3, palpha = 1){
  ggplot() + 
    geom_point(data = table_wide, 
               aes(x = .data[[s1]], y =.data[[s2]], col = .data$cluster_id), 
               size = psize,
               alpha = palpha) +
    scale_color_manual(values=color_palette)+
    theme_minimal()+
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
      data = subset(table_wide, is_driver == TRUE),
      aes(x = .data[[s1]], y = .data[[s2]]),
      color = "black", size = 0.3, shape = 15
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1)
    ) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha=1), nrow=1))+
    theme_bw() + theme(legend.position="bottom") + coord_fixed() +
    my_theme
    # theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position="bottom")
  
}

plot_scatter_driver_cartoon_marginals = function(table_wide, s1, s2, color_palette){
  p = ggplot() + 
    geom_point(data = table_wide, 
               aes(x = .data[[s1]], y =.data[[s2]], col = .data$cluster_id), 
               size = .2,
               alpha = 1) +
    scale_color_manual(values=color_palette)+
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
      data = subset(table_wide, is_driver == TRUE),
      aes(x = .data[[s1]], y = .data[[s2]]),
      color = "black", size = 1, shape = 15
    ) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha=1), nrow=2))+
    theme_bw() + theme(legend.position="bottom") + coord_fixed() + xlim(0,1) + ylim(0,1) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position="bottom")
  
  ggMarginal(p, type = "histogram",
             groupColour = TRUE,
             groupFill = TRUE,
             bins = 100,
             position="identity")
}

plot_mutations_on_tree = function(sticks_process, 
                                  sample_forest,
                                  color_palette){
  mutations_with_cell=sticks_process
  
  cells_labels_process = mutations_with_cell %>% 
    select(mutation_id, cell_id, cluster_id,sample_id) %>% 
    group_by(cell_id) %>% 
    summarise(label_list=list(cluster_id)) %>% 
    rowwise() %>% 
    mutate(label=names(sort(table(label_list[[1]]), decreasing=TRUE))[1]) %>% 
    ungroup() %>% 
    select(-label_list)
  
  final_labels = sample_forest$get_nodes() %>% as_tibble() %>% 
    left_join(cells_labels_process)
  final_labels = final_labels %>% filter(!is.na(label))
  
  pl_sticks = plot_sticks(sample_forest, labels=final_labels, cls = color_palette) %>%
    annotate_forest(sample_forest, samples=TRUE, drivers=TRUE)
  
  return(pl_sticks)
}

my_plot_forest = function (forest, sticks_process, highlight_sample = NULL, color_map = NULL, horizontal = TRUE) {
  stopifnot(inherits(forest, "Rcpp_SampleForest"))
  nodes <- forest$get_nodes() # if it has an ancestor it is a leaf node
  if (nrow(nodes) == 0) {
    warning("The forest does not contain any node")
    return(ggplot2::ggplot())
  } else {
    forest_data <- forest$get_nodes()
    forest_data[nrow(forest_data) + 1, ] <- c(NA, NA, NA, NA, NA, 0) # add the last row as NA with birth time = 0
    forest_data <- forest_data %>%
      dplyr::as_tibble() %>%
      dplyr::rename(from = .data$ancestor, to = .data$cell_id) %>%
      dplyr::select(.data$from, .data$to, .data$mutant, .data$epistate, .data$sample, .data$birth_time) %>%
      dplyr::mutate(from = ifelse(is.na(.data$from), "WT", .data$from), # WT if root
                    to   = ifelse(is.na(.data$to), "WT", .data$to), # WT the last row
                    species = paste0(.data$mutant, .data$epistate),
                    sample  = ifelse(is.na(.data$sample), "N/A", .data$sample),
                    highlight = FALSE)
    
    if (!is.null(highlight_sample)) {
      highlight <- ProCESS:::paths_to_sample(forest_data, highlight_sample)
      forest_data$highlight <- forest_data$to %in% highlight$to
    }
    
    edges <- forest_data %>% dplyr::select("from", "to", "highlight")
    graph <- tidygraph::as_tbl_graph(edges, directed = TRUE)
    graph <- graph %>%
      tidygraph::activate("nodes") %>%
      dplyr::left_join(forest_data %>%
                         dplyr::rename(name = .data$to) %>%
                         dplyr::mutate(name = as.character(.data$name)),
                       by = "name")
    
    layout <- ggraph::create_layout(graph, layout = "tree", root = "WT")
    max_Y <- max(layout$birth_time, na.rm = TRUE)
    
    layout$y <- layout$birth_time
    
    if (is.null(color_map)) {
      color_map <- ProCESS:::get_species_colors(forest$get_species_info())
    }
    
    nsamples <- forest$get_samples_info() %>% nrow()
    labels_every <- max_Y / 10
    point_size <- c(0.5, rep(2, nsamples))
    names(point_size) <- c("N/A", forest$get_samples_info() %>% dplyr::pull(.data$name))
    
    group_name <- ProCESS:::get_group_cell_name(forest)
    
    
    mutations_with_cell=sticks_process
    cells_labels_process = sticks_process %>% 
      select(mutation_id, cell_id, cluster_id,sample_id) %>% 
      group_by(cell_id) %>% 
      summarise(label_list=list(cluster_id)) %>% 
      rowwise() %>% 
      mutate(label=names(sort(table(label_list[[1]]), decreasing=TRUE))[1]) %>% 
      ungroup() %>% 
      select(-label_list) %>% 
      dplyr::rename(name=cell_id)
    
    cells_labels_process$name <- as.character(cells_labels_process$name)
    
    # mutations which are not sample have a NA in label (e.g. leaves nodes)
    my_layout = layout %>% left_join(cells_labels_process, by='name')
    
    # --- Base tree edges
    graph_plot <- ggraph::ggraph(my_layout, "tree") +
      ggraph::geom_edge_link(edge_width = 0.1,
                             ggplot2::aes(edge_color = ifelse(highlight, "indianred3", "black")))
    
    # --- Base nodes: colored by species
    point_size <- c(2, rep(2, nsamples)) # 2 1 1 1
    
    # p <- graph_plot +
    #   ggraph::geom_node_point(
    #     data = function(d) d[!(is.na(d$label) & (is.na(d$sample) | d$sample == "N/A")), ],
    #     ggplot2::aes(
    #       color = .data$label,
    #       shape = ifelse(is.na(.data$sample), "N/A", .data$sample),
    #       size  = .data$sample
    #     )
    #   )+
    #   ggplot2::scale_shape_manual(values = c(16,0,2,3)) + # https://ggplot2.tidyverse.org/reference/scale_shape.html
    #   # ggplot2::scale_shape_manual(values = c(16,15,17,18)) +
    #   ggplot2::scale_size_manual(values = point_size) +
    #   ggplot2::scale_color_manual(values = color_map) +
    #   ggplot2::theme_minimal() +
    #   ggplot2::theme(legend.position = "bottom") +
    #   ggplot2::labs(color = 'Cluster', shape = "Sample", x = NULL, y = "Time") +
    #   ggplot2::guides(size = "none",
    #                   shape = ggplot2::guide_legend("Sample"),
    #                   color = ggplot2::guide_legend('Cluster')) +
    #   ggplot2::scale_size_manual(values = point_size) +
    #   ggplot2::scale_y_continuous(labels = seq(0, max_Y, labels_every) %>% round,
    #                               breaks = seq(0, max_Y, labels_every) %>% round) +
    #   ggplot2::theme(axis.line.x = ggplot2::element_blank(),
    #                  axis.text.x = ggplot2::element_blank(),
    #                  axis.ticks.x = ggplot2::element_blank())
    
    p <- graph_plot +
      ggraph::geom_node_point(
        # data = function(d) d[!(is.na(d$label) & (is.na(d$sample) | d$sample == "N/A")), ],
        data = dplyr::filter(my_layout, (sample == "N/A") & !(is.na(label))),
        ggplot2::aes(
          color = .data$label,
          size  = .data$sample
        )
      )+
      ggplot2::scale_size_manual(values = 2) +
      ggplot2::scale_color_manual(values = color_map) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(color = 'Cluster', x = NULL, y = "Time") +
      ggplot2::guides(size = "none",
                      color = ggplot2::guide_legend('Cluster',override.aes = list(size = 3, alpha=1), nrow=2)) +
      ggplot2::scale_y_continuous(labels = seq(0, max_Y, labels_every) %>% round,
                                  breaks = seq(0, max_Y, labels_every) %>% round) +
      ggplot2::theme(axis.line.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
    
    p <- p +
      ggraph::geom_node_point(
        data = dplyr::filter(my_layout, sample != "N/A"),
        # inherit.aes = FALSE,
        ggplot2::aes(x = x, y = y, shape = sample),
        size = 2,
        # stroke = 0, # border thickness
        color='black'
      ) +
      ggplot2::scale_shape_manual(values = c(0,2,3))+
      # ggplot2::scale_shape_manual(values = c(15,17,18))+
      ggplot2::labs(shape = 'Sample', x = NULL, y = "Time") +
      ggplot2::guides(size = "none",
                      shape = ggplot2::guide_legend("Sample",override.aes = list(size = 3, alpha=1), nrow=2))
    
    
    if (horizontal) {
      p = p + scale_y_reverse()
    } else {
      p = p +
        coord_flip() +
        theme_minimal() +
        theme(
          axis.ticks.y = element_blank(),
          axis.text.y  = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = "bottom"
        )
    }
  }
  return(p)
}

my_plot_forest_cartoon = function(forest, sticks_process, highlight_sample = NULL, color_map = NULL, horizontal = TRUE) {
  stopifnot(inherits(forest, "Rcpp_SampleForest"))
  nodes <- forest$get_nodes() # if it has an ancestor it is a leaf node
  if (nrow(nodes) == 0) {
    warning("The forest does not contain any node")
    return(ggplot2::ggplot())
  } else {
    forest_data <- forest$get_nodes()
    forest_data[nrow(forest_data) + 1, ] <- c(NA, NA, NA, NA, NA, 0) # add the last row as NA with birth time = 0
    forest_data <- forest_data %>%
      dplyr::as_tibble() %>%
      dplyr::rename(from = .data$ancestor, to = .data$cell_id) %>%
      dplyr::select(.data$from, .data$to, .data$mutant, .data$epistate, .data$sample, .data$birth_time) %>%
      dplyr::mutate(from = ifelse(is.na(.data$from), "WT", .data$from), # WT if root
                    to   = ifelse(is.na(.data$to), "WT", .data$to), # WT the last row
                    species = paste0(.data$mutant, .data$epistate),
                    sample  = ifelse(is.na(.data$sample), "N/A", .data$sample),
                    highlight = FALSE)
    
    if (!is.null(highlight_sample)) {
      highlight <- ProCESS:::paths_to_sample(forest_data, highlight_sample)
      forest_data$highlight <- forest_data$to %in% highlight$to
    }
    
    edges <- forest_data %>% dplyr::select("from", "to", "highlight")
    graph <- tidygraph::as_tbl_graph(edges, directed = TRUE)
    graph <- graph %>%
      tidygraph::activate("nodes") %>%
      dplyr::left_join(forest_data %>%
                         dplyr::rename(name = .data$to) %>%
                         dplyr::mutate(name = as.character(.data$name)),
                       by = "name")
    
    layout <- ggraph::create_layout(graph, layout = "tree", root = "WT")
    max_Y <- max(layout$birth_time, na.rm = TRUE)
    
    layout$y <- layout$birth_time
    
    if (is.null(color_map)) {
      color_map <- ProCESS:::get_species_colors(forest$get_species_info())
    }
    
    nsamples <- forest$get_samples_info() %>% nrow()
    labels_every <- max_Y / 10
    point_size <- c(0.5, rep(2, nsamples))
    names(point_size) <- c("N/A", forest$get_samples_info() %>% dplyr::pull(.data$name))
    
    group_name <- ProCESS:::get_group_cell_name(forest)
    
    
    mutations_with_cell=sticks_process
    cells_labels_process = sticks_process %>% 
      select(mutation_id, cell_id, cluster_id,sample_id) %>% 
      group_by(cell_id) %>% 
      summarise(label_list=list(cluster_id)) %>% 
      rowwise() %>% 
      mutate(label=names(sort(table(label_list[[1]]), decreasing=TRUE))[1]) %>% 
      ungroup() %>% 
      select(-label_list) %>% 
      dplyr::rename(name=cell_id)
    
    cells_labels_process$name <- as.character(cells_labels_process$name)
    
    # mutations which are not sample have a NA in label (e.g. leaves nodes)
    my_layout = layout %>% left_join(cells_labels_process, by='name')
    
    # Subset the leaves
    # p <- 0.30
    # 
    # my_layout <- my_layout %>%
    #   filter(!is.na(label)) %>%
    #   bind_rows(
    #     my_layout %>%
    #       filter(is.na(label)) %>%
    #       slice_sample(prop = p)
    #   )
    
    # --- Base tree edges
    graph_plot <- ggraph::ggraph(my_layout, "tree") +
      ggraph::geom_edge_link(edge_width = 0.1,
                             ggplot2::aes(edge_color = ifelse(highlight, "indianred3", "black")))
    
    # --- Base nodes: colored by species
    point_size <- c(2, rep(2, nsamples)) # 2 1 1 1
    
    # p <- graph_plot +
    #   ggraph::geom_node_point(
    #     data = function(d) d[!(is.na(d$label) & (is.na(d$sample) | d$sample == "N/A")), ],
    #     ggplot2::aes(
    #       color = .data$label,
    #       shape = ifelse(is.na(.data$sample), "N/A", .data$sample),
    #       size  = .data$sample
    #     )
    #   )+
    #   ggplot2::scale_shape_manual(values = c(16,0,2,3)) + # https://ggplot2.tidyverse.org/reference/scale_shape.html
    #   # ggplot2::scale_shape_manual(values = c(16,15,17,18)) +
    #   ggplot2::scale_size_manual(values = point_size) +
    #   ggplot2::scale_color_manual(values = color_map) +
    #   ggplot2::theme_minimal() +
    #   ggplot2::theme(legend.position = "bottom") +
    #   ggplot2::labs(color = 'Cluster', shape = "Sample", x = NULL, y = "Time") +
    #   ggplot2::guides(size = "none",
    #                   shape = ggplot2::guide_legend("Sample"),
    #                   color = ggplot2::guide_legend('Cluster')) +
    #   ggplot2::scale_size_manual(values = point_size) +
    #   ggplot2::scale_y_continuous(labels = seq(0, max_Y, labels_every) %>% round,
    #                               breaks = seq(0, max_Y, labels_every) %>% round) +
    #   ggplot2::theme(axis.line.x = ggplot2::element_blank(),
    #                  axis.text.x = ggplot2::element_blank(),
    #                  axis.ticks.x = ggplot2::element_blank())
    
    p <- graph_plot +
      ggraph::geom_node_point(
        # data = function(d) d[!(is.na(d$label) & (is.na(d$sample) | d$sample == "N/A")), ],
        data = dplyr::filter(my_layout, (sample == "N/A") & !(is.na(label))),
        ggplot2::aes(
          color = .data$label,
          size  = .data$sample
        )
      )+
      ggplot2::scale_size_manual(values = 2) +
      ggplot2::scale_color_manual(values = color_map) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(color = 'Cluster', x = NULL, y = "Time") +
      ggplot2::guides(size = "none",
                      color = ggplot2::guide_legend('Cluster',override.aes = list(size = 3, alpha=1), nrow=2)) +
      ggplot2::scale_y_continuous(labels = seq(0, max_Y, labels_every) %>% round,
                                  breaks = seq(0, max_Y, labels_every) %>% round) +
      ggplot2::theme(axis.line.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
    
    p <- p +
      ggraph::geom_node_point(
        data = dplyr::filter(my_layout, sample != "N/A"),
        # inherit.aes = FALSE,
        ggplot2::aes(x = x, y = y, shape = sample),
        size = 2,
        # stroke = 0, # border thickness
        color='black'
      ) +
      ggplot2::scale_shape_manual(values = c(0,2,3))+
      # ggplot2::scale_shape_manual(values = c(15,17,18))+
      ggplot2::labs(shape = 'Sample', x = NULL, y = "Time") +
      ggplot2::guides(size = "none",
                      shape = ggplot2::guide_legend("Sample",override.aes = list(size = 3, alpha=1), nrow=2))
    
    
    if (horizontal) {
      p = p + scale_y_reverse()
    } else {
      p = p +
        coord_flip() +
        theme_minimal() +
        theme(
          axis.ticks.y = element_blank(),
          axis.text.y  = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = "bottom"
        )
    }
  }
  return(p)
}


# forest=sample_forest
# labels=sticks_process 
# cls = color_palette_process

my_plot_sticks = function(forest, labels, cls = NULL) {
  stopifnot(inherits(forest, "Rcpp_SampleForest"))
  nodes <- forest$get_nodes()
  if (nrow(nodes) == 0) {
    warning("The forest does not contain any node")
    return(ggplot2::ggplot())
  } else {
    labels <- labels %>%
      dplyr::mutate(
        cell_id = as.character(.data$cell_id),
        ancestor = ifelse(is.na(.data$ancestor), "WT",
                          as.character(.data$ancestor))
      ) %>%
      dplyr::add_row(cell_id="WT", ancestor="WT", mutant=NA, epistate="",
                     sample=NA, birth_time=0, label="Truncal")
    
    forest_data <- forest$get_nodes()
    
    forest_data[nrow(forest_data) + 1, ] <- c(NA, NA, NA, NA, NA, 0)
    
    forest_data <- forest_data %>% dplyr::as_tibble()  %>%
      dplyr::rename(from = .data$ancestor, to = .data$cell_id) %>%
      dplyr::select(.data$from,
                    .data$to,
                    .data$mutant,
                    .data$epistate,
                    .data$sample,
                    .data$birth_time) %>%
      dplyr::mutate(
        from = ifelse(is.na(.data$from), "WT", .data$from),
        to = ifelse(is.na(.data$to), "WT", as.character(.data$to)),
        species = paste0(.data$mutant, .data$epistate),
        sample = ifelse(is.na(.data$sample), "N/A", .data$sample),
        highlight = FALSE
      ) %>%
      dplyr::full_join(labels %>%
                         dplyr::rename(to = cell_id) %>%
                         dplyr::select(c("to", "label")),
                       by = "to") %>%
      dplyr::mutate(label = ifelse(is.na(label), "Subclonal", label))
    
    edges <- forest_data %>% dplyr::select("from", "to",
                                           "highlight")
    graph <- tidygraph::as_tbl_graph(edges, directed = TRUE)
    graph <- graph %>% tidygraph::activate("nodes") %>%
      dplyr::left_join(forest_data %>% dplyr::rename(name = .data$to) %>%
                         dplyr::mutate(name = as.character(.data$name)),
                       by = "name")
    layout <- ggraph::create_layout(graph, layout = "tree",
                                    root = "WT")
    max_Y <- max(layout$birth_time, na.rm = TRUE)
    layout$reversed_btime <- max_Y - layout$birth_time
    layout$y <- layout$reversed_btime
    
    if (is.null(cls)) {
      cls <- get_colors_for(forest_data %>% dplyr::arrange(desc(to)) %>%
                              dplyr::pull(label) %>% unique())
      cls["Subclonal"] <- "gainsboro"
    }
    nsamples <- forest$get_samples_info() %>% nrow()
    labels_every <- max_Y/10
    not_subclonal <- forest_data$label[forest_data$label != "Subclonal"] %>%
      unique()
    point_size <- c(0,rep(0.5,length(not_subclonal)))
    names(point_size) <- c("Subclonal", not_subclonal)
    
    p = ggraph::ggraph(layout, "tree") +
      ggraph::geom_edge_link(edge_width = 0.5,
                             ggplot2::aes(edge_color = ifelse(highlight,
                                                              "indianred3",
                                                              "gainsboro"))) +
      
      ggraph::geom_node_point(
        ## added 
        data = dplyr::filter(layout, (sample == "N/A") & !(is.na(label))),
        ##
        ggplot2::aes(
        color = .data$label,
        size = .data$label
      )) +
      ggplot2::scale_color_manual(values = cls, na.translate=FALSE)  +
      ggplot2::theme_minimal() + 
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(
        shape = "Sample",
        x = NULL,
        y = "Birth time"
      ) + ggplot2::guides(
        size = "none",
        shape = ggplot2::guide_legend("Sample"),
        fill = ggplot2::guide_legend("Sample")
      ) +
      ggplot2::scale_size_manual(values = point_size) +
      ggplot2::scale_y_continuous(labels = seq(0, max_Y, labels_every) %>%
                                    round %>% rev,
                                  breaks = seq(0, max_Y, labels_every) %>%
                                    round ) +
      # ggplot2::scale_y_continuous(
      #   breaks = scales::pretty_breaks(n = 2),
      #   labels = scales::pretty_breaks(n = 2)
      # )+
      
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),# remove this
        panel.grid.minor = ggplot2::element_blank(),# and this to add the grid
        axis.line.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
    
    p <- p +
      ggraph::geom_node_point(
        data = dplyr::filter(layout, .data$sample  != "N/A"),
        # inherit.aes = FALSE,
        ggplot2::aes(x = x, y = y, shape = sample),
        size = 1,
        # stroke = 0, # border thickness
        color='black'
      ) +
      ggplot2::scale_shape_manual(values = c(0,2,3))+
      # ggplot2::scale_shape_manual(values = c(15,17,18))+
      ggplot2::labs(shape = 'Sample', x = NULL, y = "Time") +
      ggplot2::guides(size = "none",
                      shape = ggplot2::guide_legend("Sample",override.aes = list(size = 2, alpha=1), nrow=2),
                      color = ggplot2::guide_legend("Population",override.aes = list(size = 1, alpha=1), nrow=2))+
      my_theme+
      theme(
        panel.border = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x  = element_blank(),
        axis.title.x = element_blank()
      )
      # theme(legend.key.size = unit(0.1, 'cm'), #change legend key size
      #       axis.title = element_text(size=8),
      #       axis.text = element_text(size=6),
      #       legend.key.height = unit(0.1, 'cm'), #change legend key height
      #       legend.key.width = unit(0.1, 'cm'), #change legend key width
      #       legend.key.spacing=unit(0.1, "cm"),
      #       legend.title = element_text(size=8), #change legend title font size
      #       legend.text = element_text(size=8)) #change legend text font size
    p
  }
}
