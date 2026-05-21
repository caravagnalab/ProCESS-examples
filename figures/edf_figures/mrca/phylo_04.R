rm(list = ls())

library(ProCESS)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(patchwork)
source("/orfeo/cephfs/scratch/cdslab/gsantacatterina/ProCESS-examples/getters/process_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/report/plotting/utils.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

# Plot tissue sampling
library(ProCESS)
library(dplyr)

seed = 777
set.seed(seed)

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
  
  samples_info <- forest$get_samples_info()
  if (samples) {
    #max_Y <- max(tree_plot$data$y, na.rm = TRUE)
    tree_plot <- tree_plot + ggplot2::geom_hline(yintercept = samples_info$time,
                                                 color = "gray30",
                                                 linetype = "dashed",
                                                 linewidth = 0.5)
  }
  if (MRCAs) {
    sample_names <- samples_info %>% dplyr::pull(.data$name)
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
  
  if (exposures) {
    if (inherits(forest, "Rcpp_PhylogeneticForest")) {
      #max_Y <- max(tree_plot$data$y, na.rm = TRUE)
      exposures <- forest$get_exposures()
      exposure_colors <- get_colors_for(exposures %>% dplyr::pull(signature) %>%
                                          unique)
      times <- exposures$time %>% unique() %>% sort()
      exposures <- exposures %>% dplyr::rowwise() %>% dplyr::mutate(t_end = dplyr::case_when(time ==
                                                                                               max(times) ~ Inf, .default = min(times[times >=
                                                                                                                                        time]))) %>% dplyr::mutate(signature = factor(signature,
                                                                                                                                                                                      levels = exposures %>% dplyr::arrange(time) %>%
                                                                                                                                                                                        dplyr::pull(signature) %>% unique()))
      tree_plot <- tree_plot + ggplot2::geom_rect(data = exposures,
                                                  ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = ifelse(is.infinite(t_end),
                                                                                                      0, max_Y - t_end), ymax = max_Y - time, fill = signature,
                                                               alpha = exposure)) + ggplot2::scale_fill_manual(values = exposure_colors) +
        ggplot2::scale_alpha_continuous(range = c(0.25,
                                                  0.75), breaks = sort(unique(exposures$exposure))) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "Signature"),
                        alpha = ggplot2::guide_legend(title = "Exposure"))
      if (facet_signatures)
        tree_plot <- tree_plot + ggplot2::facet_wrap(~signature)
      layers_new <- list(tree_plot$layers[[length(tree_plot$layers)]])
      layers_new <- c(layers_new, tree_plot$layers[1:length(tree_plot$layers) -
                                                     1])
      tree_plot$layers <- layers_new
    }
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
      
      tree_plot <- tree_plot + ggplot2::geom_point(data = layout,
                                                   ggplot2::aes(x = .data$x, y = .data$y), fill = "#d68910",
                                                   color = "#FF000000", size = 2, pch = 21)
      if (add_driver_label) {
        if ('SPN04_1.1' %in% samples_info$name){
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
        }
      }
    }
  }
  tree_plot
}

# Loading stuff ####
forest = load_sample_forest("/orfeo/scratch/cdslab/shared/SCOUT/SPN04/process/sample_forest.sff")
phylo = load_phylogenetic_forest("/orfeo/scratch/cdslab/shared/SCOUT/SPN04/process/phylo_forest.sff")
muts_path = get_mutations("SPN04", coverage = 100, purity = 0.9, type = "tumour")
muts = readRDS(muts_path)
sim = recover_simulation("/orfeo/scratch/cdslab/shared/SCOUT/SPN04/process/SPN04/")
treatment_info = readRDS("/orfeo/scratch/cdslab/shared/SCOUT/SPN04/process/treatment_info.rds")

# Plotting
library(ggplot2)

treatment_info_df = dplyr::tibble(s=treatment_info$treatment_start, e=treatment_info$treatment_end)

tree_plot = my_plot_forest(forest = forest, 
                           color_map = c(
                             "Clone 1" = "gainsboro",
                              "Clone 2" = "gainsboro",
                              "Clone 3" = "gainsboro"), 
                           color_sample = c('SPN04_1.1' = 'darkseagreen', 'SPN04_2.1' = 'cadetblue')) %>%
  my_annotate_forest(forest = forest, 
                     phylo = phylo, 
                     samples = T, 
                     MRCAs = T, 
                     drivers = T, 
                     add_driver_label = T,
                     color_sample = c('SPN04_1.1' = 'darkseagreen', 'SPN04_2.1' = 'cadetblue')) +
  guides(color = guide_legend(override.aes = list(alpha = 1)),
         fill = guide_legend(override.aes = list(alpha = 1))) + 
  geom_rect(
    data = treatment_info_df,
    aes(ymin = s, ymax = e, xmin = -Inf, xmax = Inf),
    inherit.aes = FALSE,
    alpha = 0.2,
    fill = "coral"
  )  + 
  guides(shape = "none") 

tree_plot_spn04 = tree_plot + theme(legend.direction = "horizontal")
# ggsave(filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/mrca/forest_mrca.png',
#        plot = tree_plot_spn04, width = 7, height = 8, units = 'in', dpi = 600)
# ggsave(filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/mrca/forest_mrca.pdf',
#        plot = tree_plot_spn04, width = 7, height = 8, units = 'in')


cell_theraphy <- phylo$get_nodes() %>% 
  filter(birth_time > 302.4552, birth_time < 344.6101) %>% 
  filter(is.na(sample) | sample == 'SPN04_2.1') 

cell_surviving <- cell_theraphy %>% filter(!ancestor %in% cell_theraphy$cell_id)
cell_ancestor <- cell_surviving$ancestor

s = 'SPN04_2.1'
forest = phylo
mrca <- forest$get_coalescent_cells(
  forest$get_nodes() %>%
    dplyr::filter(sample %in% s) %>%
    dplyr::pull(.data$cell_id)
) %>%
  dplyr::mutate(sample = s) 

sampled_cell <- phylo$get_nodes() %>% select(cell_id, ancestor) 
head(sampled_cell)


# get path
# Build a named lookup: cell_id -> ancestor
ancestor_map <- setNames(sampled_cell$ancestor, sampled_cell$cell_id)

mrca_id <- mrca$cell_id  # 29208

# Function: walk one cell up to the MRCA, returning the full path
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

# Trace all 9 start cells
all_paths <- lapply(cell_ancestor, 
                    trace_to_mrca,
                    ancestor_map = ancestor_map,
                    mrca_id = mrca_id)

names(all_paths) <- as.character(cell_ancestor)

path_table <- do.call(rbind, lapply(names(all_paths), function(start_id) {
  path <- all_paths[[start_id]]
  data.frame(
    cell_id = path,
    label   = start_id,
    stringsAsFactors = FALSE
  )
}))

path_table <- path_table %>% filter(cell_id != mrca$cell_id)


rownames(path_table) <- NULL
head(path_table)
path_table_unique <- path_table[!duplicated(path_table$cell_id), ]

path_table_unique <- path_table_unique %>% mutate(cell_id = paste(cell_id),
                                                  label = paste0('cell_', label))



base = '/orfeo/cephfs/scratch/cdslab/shared/SCOUT/validation_subclonal_new/tables/'
gt_f =  readRDS(paste0(base, 'table_process_univariate_w_private_', 'SPN04', '_', '150', 'x_', '0.9', 'p_mutect2_ascat.rds')) %>% 
  filter(!stringr::str_detect(causes, "errors"))  %>% 
  select(sample_id, mutation_id, cell_id, causes, is_driver_process, cluster_id_process, vaf_process, ccf_process) %>% 
  filter(sample_id == 'SPN04_SPN04_2.1') %>% 
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
  
treatment_cell = filt_gt %>% filter(causes == 'SBS25') %>% pull(cell_id) 
treatment_cell_df = tibble(cell_id = paste(treatment_cell), tretment_label = 'SBS25')


layout <- tree_plot$data %>% 
  dplyr::select(.data$x, .data$y,.data$name, sample) %>%
  dplyr::mutate(cell_id = paste(.data$name)) %>%
  dplyr::filter(.data$name %in% c(path_table_unique$cell_id, treatment_cell_df$cell_id)) %>%
  dplyr::left_join(path_table_unique, by = "cell_id") %>% 
  left_join(treatment_cell_df) %>% 
  mutate(final_label = ifelse(is.na(tretment_label), label, tretment_label))
layout$final_label <- trimws(layout$final_label)



labels_color = c('indianred2', 'steelblue', 'goldenrod', 
                 'darkorange', 'mediumpurple', 'seagreen',
                 'darkolivegreen2', 'SBS25' = 'maroon3')
names(labels_color) = unique(layout$final_label)
labels_color = c(labels_color, c('SPN04_1.1' = 'darkseagreen', 'SPN04_2.1' = 'cadetblue'))

tree_plot_path <- tree_plot +
  ggplot2::geom_point(
    data = layout,
    aes(x = x, y = y, color = as.factor(final_label)),
    size = 2,
    alpha = .5
  ) +
  ggplot2::scale_color_manual(
    values = labels_color)  +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

filt_gt <- filt_gt %>% mutate(cell_id = paste(cell_id))

plt_vaf <- filt_gt %>% 
  left_join(path_table_unique) %>% 
  mutate(new_label = ifelse(!is.na(label), label, cluster_id_process))  %>% 
  ggplot()+
  geom_histogram(aes(x = vaf_process, fill = new_label), alpha=0.6, binwidth = 0.007, position = "identity") +
  my_ggplot_theme() +
  xlim(0.02, 0.7) +
  scale_fill_manual('Lineages', values = c('Clonal' = 'antiquewhite3', 'Other' = 'gainsboro', labels_color))  +
  xlab('VAF SPN04 2.1') +
  ylab('')+
  guides(fill = guide_legend(override.aes = list(alpha = 1)))
  #theme(legend.position = 'none') 

plt_vaf_treatment <- filt_gt %>% 
  mutate(new_label = ifelse(causes == 'SBS25', causes, 'Other'))  %>% 
  ggplot() + 
  geom_histogram(aes(x = vaf_process, fill = new_label), alpha=0.9, binwidth = 0.01, position = "identity") +
  my_ggplot_theme() +
  xlim(0.02, 0.7) +
  scale_fill_manual('Signature',values = c(labels_color, 'Other' = 'gainsboro')) +
  xlab('VAF SPN04 2.1')+
  ylab('')+
  guides(fill = guide_legend(override.aes = list(alpha = 1)))


p <- pp + labs(tag = 'A') +
  ggplot() + labs(tag = 'B') +
  tree_plot_path +  labs(tag ='C')  + 
  plt_vaf + theme(legend.position = 'bottom') +   labs(tag ='D') +
  plt_vaf_treatment + theme(legend.position = 'bottom') +   labs(tag ='E') +
  plot_layout(design = 'AABBB\nAABBB\nCCCDD\nCCCDD\nCCCEE\nCCCEE')
   
ggsave(filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/mrca/forest_mrca.pdf',
       plot = p, width = 11, height = 11, units = 'in')




  