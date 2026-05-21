library(tidyverse)
source("/orfeo/cephfs/scratch/cdslab/gsantacatterina/ProCESS-examples/getters/process_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/report/plotting/utils.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

my_plot_forest = function (forest, highlight_sample = NULL, color_map = NULL, horizontal = T, 
                           color_sample = c('SPN04_1.1' = 'darkseagreen', 'SPN04_2.1' = 'sienna3')) {
  stopifnot(inherits(forest, "Rcpp_SampleForest"))
  nodes <- forest$get_nodes()
  if (nrow(nodes) == 0) {
    warning("The forest does not contain any node")
    return(ggplot2::ggplot())
  } else {
    forest_data <- forest$get_nodes()
    forest_data[nrow(forest_data) + 1, ] <- c(NA, NA, NA, NA, NA, 0)
    forest_data <- forest_data %>%
      dplyr::as_tibble() %>%
      dplyr::rename(from = .data$ancestor, to = .data$cell_id) %>%
      dplyr::select(.data$from, .data$to, .data$mutant, .data$epistate, .data$sample, .data$birth_time) %>%
      dplyr::mutate(from = ifelse(is.na(.data$from), "WT", .data$from),
                    to   = ifelse(is.na(.data$to), "WT", .data$to),
                    species = paste0(.data$mutant, .data$epistate),
                    sample  = ifelse(is.na(.data$sample), "N/A", .data$sample),
                    highlight = FALSE)
    
    if (!is.null(highlight_sample)) {
      highlight <- paths_to_sample(forest_data, highlight_sample)
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
    point_size <- c(0.1, rep(1, nsamples))
    names(point_size) <- c("N/A", forest$get_samples_info() %>% dplyr::pull(.data$name))
    
    group_name <- ProCESS:::get_group_cell_name(forest)
    
    # --- Base tree edges
    graph_plot <- ggraph::ggraph(layout, "tree") +
      ggraph::geom_edge_link(edge_width = 0.1,
                             ggplot2::aes(edge_color = ifelse(highlight, "indianred3", "black")))
    
    # --- Base nodes: colored by species
    p <- graph_plot +
      ggraph::geom_node_point(
        ggplot2::aes(color = .data$species,
                     shape = ifelse(is.na(.data$sample), "N/A", .data$sample),
                     size  = .data$sample), show.legend = F
      ) +
      ggplot2::scale_shape_manual(values = c(0:nsamples + 1)) +
      ggplot2::scale_color_manual(values = color_map) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(color = group_name, shape = "Sample", x = NULL, y = "Time") +
      ggplot2::guides(size = "none",
                      shape = ggplot2::guide_legend("Sample"),
                      color = ggplot2::guide_legend(group_name)) +
      ggplot2::scale_size_manual(values = point_size) +
      ggplot2::scale_y_continuous(labels = seq(0, max_Y, labels_every) %>% round,
                                  breaks = seq(0, max_Y, labels_every) %>% round) +
      ggplot2::theme(axis.line.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
    
    # --- NEW LAYER: outline circles for Sample != "N/A"
    p <- p +
      ggraph::geom_node_point(
        data = dplyr::filter(layout, sample != "N/A"),
        inherit.aes = FALSE,
        ggplot2::aes(x = x, y = y, fill = sample, color = sample),  # <-- map sample to fill
        pch = 21,          # circle with fill
        #stroke = 0,      # border thickness
        #color = "white",   # border color stays black
        size = 1,show.legend = F,
        alpha = .5# fixed circle size
      ) +
      ggplot2::scale_fill_manual(
        values = color_sample
      ) + 
      ggplot2::scale_color_manual(
        values = color_sample
      )
    
    p = p +
      labs(color = "Clone", fill = "Sample")
    
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
    
    return(p)
  }
}

my_annotate_forest = function (tree_plot, 
                               forest,
                               phylo = NULL, 
                               samples = TRUE, MRCAs = TRUE, exposures = FALSE,
                               facet_signatures = TRUE,
                               drivers = TRUE, add_driver_label = TRUE,
                               color_sample = c('    SPN04_1.1' = 'darkseagreen', '    SPN04_2.1' = 'sienna3')) {
  
  samples_info <- forest$get_samples_info() %>% filter(name %in% names(color_sample))
  if (samples) {
    #max_Y <- max(tree_plot$data$y, na.rm = TRUE)
    tree_plot <- tree_plot + ggplot2::geom_hline(yintercept = samples_info$time,
                                                 color = "gray30",
                                                 linetype = "dashed",
                                                 linewidth = 0.5)
  }
  if (MRCAs) {
    sample_names <- names(color_sample)#samples_info %>% dplyr::pull(.data$name) 
    MRCAs_cells <- lapply(sample_names, function(s) {
      forest$get_coalescent_cells(forest$get_nodes() %>%
                                    dplyr::filter(sample %in% s) %>%
                                    dplyr::pull(.data$cell_id)) %>%
        dplyr::mutate(sample = s)
    }) %>% Reduce(f = dplyr::bind_rows) %>% dplyr::group_by(.data$cell_id) %>%
      dplyr::mutate(cell_id = paste(.data$cell_id)) %>%
      dplyr::summarise(label = paste0("    ", .data$sample,
                                      collapse = "\n"))
    
    layout <- tree_plot$data %>% dplyr::select(.data$x, .data$y,
                                               .data$name, sample) %>%
      dplyr::mutate(cell_id = paste(.data$name)) %>%
      dplyr::filter(.data$name %in% MRCAs_cells$cell_id) %>%
      dplyr::left_join(MRCAs_cells, by = "cell_id")
    layout$label <- trimws(layout$label)
    
    # tree_plot <- tree_plot +
    #   ggplot2::geom_point(data = layout, ggplot2::aes(x = .data$x, y = .data$y), color = "purple3",size = 3, pch = 21) +
    #   ggplot2::geom_text(data = layout, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label), color = "purple3", size = 3, hjust = 0, vjust = 1, check_overlap = TRUE)
    
    tree_plot <- tree_plot +
      ggplot2::geom_point(
        data = layout,
        aes(x = x, y = y, color = label),
        #color = "purple3",
        size = 3,
        pch = 21,
        show.legend = F
      ) +
      ggplot2::scale_color_manual(
        values = color_sample) + 
      #ggrepel::geom_text_repel(
      ggrepel::geom_label_repel(
        data = layout,
        aes(x = x, y = y, label = label, color = label),
        show.legend = F,
        #color = "purple3",
        size = 3,
        max.overlaps = Inf,    # allow flexible placement
        box.padding = 0.3,     # padding around text
        point.padding = 0.2,   # padding around the plotted point
        force = 2,             # repelling force
        min.segment.length = 0 # always draw connecting lines
      ) +
      ggplot2::scale_color_manual(
        values = color_sample) 
  }
  
  if (drivers) {
    if (inherits(phylo, "Rcpp_PhylogeneticForest")) {
      drivers_mutations = drivers_CNAs = data.frame()
      try(expr = {
        drivers_mutations <- phylo$get_sampled_cell_mutations() %>%
          dplyr::filter(class == "driver") %>%
          dplyr::mutate(driver_id = paste0(chr,":", chr_pos, ":", ref, ">", alt), driver_type = type) %>%
          dplyr::select(cell_id, driver_id, driver_type)
      })
      try(expr = {
        drivers_CNAs <- phylo$get_sampled_cell_CNAs() %>%
          dplyr::mutate(driver_id = paste0(chr, ":", begin, "-", end, ":", allele), driver_type = "CNA") %>%
          dplyr::select(cell_id, driver_id, driver_type)
      })
      drivers <- dplyr::bind_rows(drivers_mutations, drivers_CNAs)
      drivers_start_nodes <- lapply(unique(drivers$driver_id),
                                    function(d) {
                                      nodes_with_driver = drivers %>%
                                        dplyr::filter(driver_id ==d) %>%
                                        dplyr::pull(cell_id)
                                      d_type = drivers %>%
                                        dplyr::filter(driver_id == d) %>%
                                        dplyr::pull(driver_type) %>%
                                        unique()
                                      
                                      forest$get_coalescent_cells(nodes_with_driver) %>%
                                        dplyr::mutate(driver_id = d, driver_type = d_type)
                                    }) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(cell_id = as.character(cell_id)) %>%
        dplyr::group_by(cell_id) %>%
        dplyr::summarise(driver_id = paste0(driver_id, collapse = "\n"))
      
      layout <- tree_plot$data %>%
        dplyr::select(x, y, name) %>%
        dplyr::mutate(cell_id = paste(name), has_driver = TRUE) %>%
        dplyr::filter(name %in% drivers_start_nodes$cell_id) %>%
        dplyr::left_join(drivers_start_nodes, by = "cell_id") 
      
      if (add_driver_label) {
        if ('SPN04_1.1' %in% samples_info$name){
          
          tree_plot <- tree_plot + ggplot2::geom_point(data = layout,
                                                       ggplot2::aes(x = .data$x, y = .data$y), fill = "#d68910",
                                                       color = "#FF000000", size = 2, pch = 21)
          
          layout  = layout %>% mutate(driver_label = case_when(
            cell_id == 24633 ~ 'IDH1 R123C',
            cell_id == 24946 ~ 'MYC AMP',
            cell_id == 168798 ~ 'NRAS Q61K'
          ))
          nudge_x = (max(tree_plot$data$x) - min(tree_plot$data$x)) *
            0.15
          tree_plot <- tree_plot + ggrepel::geom_label_repel(data = layout,
                                                             ggplot2::aes(x = .data$x, y = .data$y, label = .data$driver_label),
                                                             color = "#d68910", 
                                                             size = 3, 
                                                             max.overlaps = Inf,    # allow flexible placement
                                                             box.padding = 0.3,     # padding around text
                                                             point.padding = 0.2,   # padding around the plotted point
                                                             force = 2,             # repelling force
                                                             min.segment.length = 0, # always draw connecting lines
                                                             direction = "x")  
        } else if ('SPN06_1.1' %in% samples_info$name){
          layout <- layout %>% 
              filter(cell_id %in% c(3188, 29887, 296224, 207734))
          
          tree_plot_n <- tree_plot + ggplot2::geom_point(data = layout,
                                                       ggplot2::aes(x = .data$x, y = .data$y), fill = "maroon3",
                                                       color = "#FF000000", size = 2, pch = 21)
          
          layout  = layout %>% mutate(driver_label = case_when(
            cell_id == 3188 ~ 'TP53 R248W',
            cell_id == 29887 ~ 'STK11 LOH',
            cell_id == 207734 ~ 'EGFR AMP',
            cell_id == 296224 ~ 'KRAS G12D'
          ))
          nudge_x = (max(tree_plot$data$x) - min(tree_plot$data$x)) *
            0.15
          tree_plot_n <- tree_plot_n + ggrepel::geom_label_repel(data = layout,
                                                             ggplot2::aes(x = .data$x, y = .data$y, label = .data$driver_label),
                                                             color = "maroon3", 
                                                             size = 3, 
                                                             max.overlaps = Inf,    # allow flexible placement
                                                             box.padding = 0.3,     # padding around text
                                                             point.padding = 0.2,   # padding around the plotted point
                                                             force = 2,             # repelling force
                                                             min.segment.length = 0, # always draw connecting lines
                                                             direction = "x")  
        }
      }
    }
  }
  tree_plot_n
}


forest = load_sample_forest("/orfeo/scratch/cdslab/shared/SCOUT/SPN06/process/sample_forest.sff")
phylo = load_phylogenetic_forest("/orfeo/scratch/cdslab/shared/SCOUT/SPN06/process/phylo_forest.sff")
sub_forest = forest$get_subforest_for(c('SPN06_1.1', 'SPN06_2.1'))
treatment_info = readRDS("/orfeo/scratch/cdslab/shared/SCOUT/SPN06/process/chemo_timing.rds")
treatment_info_df = dplyr::tibble(s=c(treatment_info$chemo1_start),
                                  e=c(treatment_info$chemo1_end))

tree_plot = my_plot_forest(sub_forest,
                           color_map = c(
                             "Clone 1" = "gainsboro",
                             "Clone 2" = "gainsboro",
                             "Clone 3" = "gainsboro",
                             "Clone 4" = "gainsboro",
                             "Clone 5" = "gainsboro"
                           ), color_sample = c('SPN06_1.1' = 'goldenrod', 'SPN06_2.1' = 'cornflowerblue')) %>%
  my_annotate_forest(forest = forest, 
                     phylo = phylo, 
                     samples = T, 
                     MRCAs = T, 
                     drivers = T, 
                     add_driver_label = T,
                     color_sample = c('SPN06_1.1' = 'goldenrod', 'SPN06_2.1' = 'cornflowerblue')) +
  geom_rect(
    data = treatment_info_df,
    aes(ymin = s, ymax = e, xmin = -Inf, xmax = Inf),
    inherit.aes = FALSE,
    alpha = 0.2,
    fill = "coral"
  )  +
  guides(shape = "none")


base = '/orfeo/cephfs/scratch/cdslab/shared/SCOUT/validation_subclonal_new/tables/'
gt_f =  readRDS(paste0(base, 'table_process_univariate_w_private_', 'SPN06', '_', '150', 'x_', '0.9', 'p_mutect2_ascat.rds')) %>% 
  filter(!stringr::str_detect(causes, "errors"))  %>% 
  select(sample_id, mutation_id, cell_id, causes, is_driver_process, cluster_id_process, vaf_process, ccf_process) %>% 
  filter(sample_id == 'SPN06_SPN06_2.1') %>% 
  group_by(cluster_id_process, sample_id) %>%
  mutate(is_clonal_process=replace(FALSE, ccf_process > 0.95, TRUE)) %>% 
  ungroup() %>%
  mutate(cluster_id_process_full = cluster_id_process) %>%
  mutate(cluster_id_process = replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal')) %>% 
  tidyr::separate(cluster_id_process, sep = '_', into = c('cluster_id_process', 'tmp')) %>% 
  select(-tmp) 

gt_f <- gt_f %>% 
  mutate(cluster_id_process = ifelse(cluster_id_process == 'Clonal', 'Clonal', 'Other'))

filt_gt <- gt_f %>% 
  filter(!grepl(":(X|Y):", mutation_id))

treatment_cell = filt_gt %>% filter(causes == 'SBS31') %>% pull(cell_id) 
treatment_cell_df = tibble(cell_id = paste(treatment_cell), tretment_label = 'SBS31') %>% distinct()

filt_gt <- filt_gt %>% mutate(cell_id = paste(cell_id))


s = 'SPN06_1.1'
cells <- phylo$get_nodes()
cell_in_sample <- phylo$get_nodes() %>% filter(sample == 'SPN06_1.1')  %>% pull(cell_id)
mrca_id <- forest$get_coalescent_cells(
  forest$get_nodes() %>%
    dplyr::filter(sample %in% s) %>%
    dplyr::pull(.data$cell_id)
) %>%
  dplyr::mutate(sample = s)  %>% 
  pull(cell_id)

ancestor_map <- setNames(cells$ancestor, cells$cell_id)

trace_to_mrca <- function(start_id, ancestor_map, mrca_id, max_steps = 10000) {
  path <- start_id
  current <- start_id
  
  for (i in seq_len(max_steps)) {
    if (current == mrca_id) break
    
    parent <- ancestor_map[as.character(current)]
    
    if (is.na(parent)) {
      warning(paste("Reached root without finding MRCA from cell", start_id))
      break
    }
    
    current <- parent
    path <- c(path, current)
  }
  
  return(path)
}

all_paths <- lapply(cell_samples, 
                    trace_to_mrca,
                    ancestor_map = ancestor_map,
                    mrca_id = mrca_id)

paths_cell_in_sample <- all_paths %>% unlist() %>% unique()


s = 'SPN06_2.1'
cell_id <- forest$get_coalescent_cells(
  forest$get_nodes() %>%
    dplyr::filter(sample %in% s) %>%
    dplyr::pull(.data$cell_id)
) %>%
  dplyr::mutate(sample = s)  %>% 
  pull(cell_id)


all_paths_s2 <- lapply(cell_id, 
                    trace_to_mrca,
                    ancestor_map = ancestor_map,
                    mrca_id = mrca_id) %>% unlist()

private_s2 <- all_paths_s2[!all_paths_s2 %in% paths_cell_in_sample]



plt_vaf_treatment <- filt_gt %>% 
  mutate(new_label = ifelse(causes == 'SBS31', causes, 'Other'))  %>% 
  ggplot() + 
  geom_histogram(aes(x = vaf_process, fill = new_label), alpha=0.9, binwidth = 0.01, position = "identity") +
  my_ggplot_theme() +
  xlim(0.02, 1) +
  scale_fill_manual('Signature',values = c(sbs_colors, 'Other' = 'gainsboro')) +
  xlab('VAF SPN06 2.1')+
  ylab('')+
  guides(fill = guide_legend(override.aes = list(alpha = 1)))


plt_vaf <- filt_gt %>% 
  mutate(label = ifelse(cell_id %in% private_s2, 'Private SPN06 2.1', cluster_id_process)) %>% 
  #mutate(new_label = ifelse(causes == 'SBS31', causes, 'Other'))  %>% 
  ggplot() + 
  geom_histogram(aes(x = vaf_process, fill = label), alpha=0.9, binwidth = 0.01, position = "identity") +
  my_ggplot_theme() +
  xlim(0.02, 1) +
  scale_fill_manual('Signature',values =  c('Clonal' = 'palegreen4', 'Other' = 'gainsboro', 'Private SPN06 2.1' = 'palegreen3')) +
  xlab('VAF SPN06 2.1')+
  ylab('')+
  guides(fill = guide_legend(override.aes = list(alpha = 1)))

private_cell_df = tibble(cell_id = paste(private_s2), label_cell = 'Private SPN06 2.1')

layout <- tree_plot$data %>% 
  dplyr::select(.data$x, .data$y,.data$name, sample) %>%
  dplyr::mutate(cell_id = paste(.data$name)) %>%
  dplyr::filter(.data$name %in% c(treatment_cell_df$cell_id, private_cell_df$cell_id)) %>%
  left_join(treatment_cell_df) %>% 
  left_join(private_cell_df) %>% 
  mutate(label = ifelse(is.na(tretment_label), label_cell, tretment_label))
layout$label <- trimws(layout$label)


sbs_colors <- c(sbs_colors, c('SPN06_1.1' = 'goldenrod', 'SPN06_2.1' = 'cornflowerblue', 'Private SPN06 2.1' = 'palegreen3'))

tree_plot_signature <- tree_plot +
  ggplot2::geom_point(
    data = layout,
    aes(x = x, y = y, color = as.factor(label)),
    size = 2,
    alpha = .6
  ) +
  ggplot2::scale_color_manual(
    values = sbs_colors)  +
  guides(color = guide_legend(override.aes = list(alpha = 1)))



p <- tree_plot_signature +  labs(tag ='A')  + 
  plt_vaf + theme(legend.position = 'bottom') +   labs(tag ='B') +
  plt_vaf_treatment + theme(legend.position = 'bottom') +   labs(tag ='C') +
  plot_layout(design = 'AAABB\nAAABB\nAAACC\nAAACC')

ggsave(filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/mrca/forest_mrca_spn06.pdf',
       plot = p, width = 11, height = 8, units = 'in')
