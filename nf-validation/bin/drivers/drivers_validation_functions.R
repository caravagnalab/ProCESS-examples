# validation of driver calls
get_drivers_results = function(path = '/orfeo/cephfs/scratch/cdslab/shared/SCOUT', 
                               spn, 
                               samples, 
                               purity, 
                               coverage, 
                               callers, 
                               cohort) {
  
  drivers_path = paste0(path, '/', spn, '/tumourevo/', coverage, 'x_', purity, 'p_', callers, '/driver_annotation/annotate_driver/', cohort, '/', spn)
  drivers_path = lapply(samples, function(x) {
    list.files(drivers_path, pattern = x, full.names = T, recursive = T)
  })
  
  drivers = lapply(drivers_path, readRDS)
  names(drivers) = samples
  return(drivers)
}

# get process sequencing
get_process_seq = function(path = '/orfeo/cephfs/scratch/cdslab/shared/SCOUT', 
                           spn, 
                           coverage, 
                           purity
) {
  seq_path = paste0(path, '/', spn, '/sequencing/tumour/purity_', purity, '/data/mutations/seq_results_muts_merged_coverage_', coverage, 'x.rds')
  seq_obj = readRDS(seq_path)
  return(seq_obj)
}

##########################################################################################################################################################################
# heatmap of process drivers found in the tumourevo results 

# 1. preprocessing of data: selecting process drivers

# get_process_drivers_ids = function(process_seq_res) {
#   process_drivers = process_seq_res %>% 
#     dplyr::filter(classes == 'driver') %>% 
#     dplyr::mutate(chr = paste0('chr', chr)) %>% 
#     # dplyr::mutate(mut_id = paste0(chr, ':', chr_pos, '_', ref, '>', alt))
#     dplyr::mutate(mut_id = paste0(chr, ':', chr_pos))
#   
#   process_drivers_ids = process_drivers$mut_id
#   return(process_drivers_ids)
# }

get_process_drivers_codes = function(phylo_forest) {
  phylo_forest$get_driver_mutations() %>% 
    dplyr::filter(type == 'SID') %>% 
    dplyr::mutate(code = gsub(' ', '_p.', code)) %>% 
    dplyr::pull(code)
}


# this will be needed later --> filter the process results to keep only drivers -- get both mut ids and codes
get_process_drivers = function(process_seq_res, phylo_forest) {
  process_drivers = process_seq_res %>% 
    dplyr::filter(classes == 'driver') 
    # dplyr::mutate(chr = paste0('chr', chr)) %>% 
    # dplyr::mutate(mut_id = paste0(chr, ':', chr_pos, '_', ref, '>', alt))
    # dplyr::mutate(mut_id = paste0(chr, ':', chr_pos))
  dr = phylo_forest$get_driver_mutations() %>% 
    dplyr::filter(type == 'SID') %>% 
    dplyr::mutate(code = gsub(' ', '_p.', code))
  
  process_drivers = full_join(process_drivers, dr, by = join_by(
    'chr' == 'chr', 
    'chr_pos' == 'start', 
    'ref' == 'ref', 
    'alt' == 'alt'
  ))
  
  return(process_drivers)
}

# select tumourevo drivers 

# select from tumourevo the driver mutations from process
get_process_drivers_in_tumourevo = function(tumourevo_mutations, process_drivers_ids) {
  lapply(names(tumourevo_mutations), function(x) {
    tumourevo_mutations[[x]][[1]]$mutations %>% 
      # dplyr::mutate(mut_id = paste0(chr, ':', from, '_', ref, '>', alt)) %>% 
      # dplyr::mutate(mut_id = paste0(chr, ':', from)) %>% 
      # dplyr::filter(mut_id %in% process_drivers_ids) %>% 
      dplyr::filter(driver_label %in% process_drivers_ids) %>% 
      dplyr::mutate(sample = x) %>% 
      dplyr::mutate(origin = 'tumourevo')
  }) %>% 
    dplyr::bind_rows()
} %>% 
  dplyr::select(driver_label, SYMBOL, VAF, sample, origin) %>%
  tidyr::pivot_wider(names_from = c(sample, origin), names_sep = '-', values_from = VAF)


# creating the input for the heatmap plot

create_true_drivers_table = function(process_drivers,
                                     tumourevo_mutations, 
                                     process_drivers_ids){
  
  process_drivers = process_drivers %>% 
    dplyr::filter(!grepl("errors",causes)) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(code, dplyr::ends_with('VAF')) %>% 
    dplyr::rename_with(~ gsub('.VAF', '-process', .x))
  
  process_in_tumourevo = get_process_drivers_in_tumourevo(tumourevo_mutations, process_drivers_ids)
  
  # function to create the heatmap
  true_drivers = full_join(process_in_tumourevo, process_drivers, by = join_by('driver_label' == 'code')) %>% 
    tibble::column_to_rownames('driver_label')
  return(true_drivers)
  
}

# create annotation for the heatmap
create_annotation = function(x, # data 
                             ann_colors, 
                             position) {
  
  ComplexHeatmap::HeatmapAnnotation(df = x, 
                                    col = ann_colors, 
                                    which = position, 
                                    annotation_legend_param = list(nrow = 3, width = 12, by_row = T))
  
}

# create the heatmap 
col= c('white',ggsci::pal_bs5('cyan', alpha = 0.7)(10))
cols = circlize::colorRamp2(colors = col, seq(0,1, 0.1))

plot_drivers_heatmap = function(true_drivers_table, 
                                samples_palette = awtools::ppalette, # palette from which subsampling
                                origin_palette = c('#FFC107', '#347433'), # palette from which subsampling
                                # genes_palette, # palette from which subsampling
                                heatmap_cols = cols) {
  
  options(digits = 4)
  # prepare the data
  true_drivers_table = true_drivers_table %>% 
    dplyr::mutate(SYMBOL = ifelse(is.na(SYMBOL), 'Unknown driver', SYMBOL))
  
  genes = true_drivers_table$SYMBOL
  
  true_drivers = true_drivers_table %>% 
    dplyr::select(-SYMBOL) %>% 
    as.matrix()
  
  # create the heatmap annotations
  
  # rows - genes 
  row_colors = list(
    genes = setNames(nm = genes, # still a vector
                     object = viridis::viridis(length(genes)))
  )
  # convert the vector of genes to a tibble 
  genes = tibble(genes = genes) 
  
  # generate annotation
  rows_ann = create_annotation(genes, ann_colors = row_colors, position = 'row')
  
  # column - samples and origin 
  ann_info_columns = tibble(
    x = colnames(true_drivers)
  ) %>% 
    tidyr::separate(x, into = c('sample', 'origin'), sep = '-')
  
  # colors 
  column_colors = list(
    sample = setNames(nm = ann_info_columns$sample %>% unique, 
                      object = samples_palette[1:length(ann_info_columns$sample %>% unique)]), 
    origin = setNames(nm = c('process', 'tumourevo'), 
                      origin_palette)
  )
  # generate annotation
  column_ann = create_annotation(ann_info_columns, ann_colors = column_colors, position = 'column')
  
  drivers_heatmap = Heatmap(true_drivers, 
                            bottom_annotation = column_ann, 
                            left_annotation = rows_ann, 
                            name = 'VAF', 
                            row_title = 'driver mutations', 
                            column_title = 'samples', 
                            row_dend_side = 'right', 
                            row_names_side = 'left', 
                            col = heatmap_cols, 
                            heatmap_legend_param = list(
                              legend_direction = "horizontal")
  )
  
  draw(drivers_heatmap, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
}


# plot comparison 

# plot_comparison_drivers = function(x) {
#   x = x %>% 
#     dplyr::mutate(causes = ifelse(is.na(causes), 'Tumourevo annotation', causes)) %>% 
#     dplyr::select(-c(to)) %>% 
#     dplyr::mutate(across(matches('NV|DP|VAF'), ~ ifelse(is.na(.), 0, .)))
#   
#   
#   x %>% 
#     ggplot(aes(x = tumourevo_VAF, y = process_VAF, colour =  classes, label = driver_label)) + 
#     geom_point() + 
#     facet_wrap(vars(causes), scales = 'free') + 
#     theme_bw() + 
#     xlim(c(-0.01, 1.01)) + 
#     ylim(c(-0.01, 1.01)) + 
#     ggrepel::geom_text_repel(max.overlaps = 1)
#   
# }

##########################################################################################################################################################################

# select tumourevo drivers ids

# get_tumourevo_drivers_ids = function(tumourevo_mutations) {
#   lapply(names(tumourevo_mutations), function(x) {
#     tumourevo_mutations[[x]][[1]]$mutations %>% 
#       dplyr::mutate(mut_id = paste0(chr, ':', from, '_', ref, '>', alt)) %>% 
#       dplyr::filter(is_driver, VAF > 0) %>% 
#       dplyr::mutate(sample = x) %>% 
#       dplyr::select(mut_id, driver_label, sample, is_driver, VAF) %>% 
#       dplyr::rename(is_driver_tumourevo = is_driver)
#   }) %>% 
#     dplyr::bind_rows() %>% 
#     dplyr::pull(mut_id) %>% 
#     unique
# }

get_tumourevo_drivers_codes = function(tumourevo_mutations) {
  lapply(names(tumourevo_mutations), function(x) {
    tumourevo_mutations[[x]][[1]]$mutations %>% 
      # dplyr::mutate(mut_id = paste0(chr, ':', from, '_', ref, '>', alt)) %>% 
      dplyr::filter(is_driver, VAF > 0) %>% 
      dplyr::mutate(sample = x) %>% 
      dplyr::select(driver_label, sample, is_driver, VAF) %>% 
      dplyr::rename(is_driver_tumourevo = is_driver)
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::pull(driver_label) %>% 
    unique
}


get_all_drivers_process = function(process_drivers, all_drivers) {
  process_drivers_all = process_drivers %>% 
    dplyr::mutate(chr = paste0('chr', chr)) %>% 
    # dplyr::mutate(mut_id = paste0(chr, ':', chr_pos, '_', ref, '>', alt)) %>% 
    dplyr::filter(code %in% all_drivers) %>%
    dplyr::ungroup() %>% 
    dplyr::select(code, classes, ends_with('VAF')) %>% 
    reshape2::melt() %>% 
    dplyr::rename(sample = variable) %>% 
    dplyr::rename(VAF = value) %>% 
    dplyr::mutate(sample = gsub('.VAF', '', sample)) %>% 
    dplyr::mutate(classes = ifelse((classes == 'driver' & VAF > 0), TRUE, FALSE)) %>% 
    dplyr::rename(is_driver_process = classes)
  
  return(process_drivers_all)
}

get_all_drivers_tumourevo = function(tumourevo_mutations, all_drivers) {
  lapply(names(tumourevo_mutations), function(x) {
    tumourevo_mutations[[x]][[1]]$mutations %>% 
      # dplyr::mutate(mut_id = paste0(chr, ':', from, '_', ref, '>', alt)) %>% 
      dplyr::filter(driver_label %in% all_drivers) %>% 
      dplyr::mutate(sample = x) %>% 
      # dplyr::filter(VAF > 0) %>%
      dplyr::select(driver_label, sample, is_driver, VAF) %>% 
      dplyr::rename(is_driver_tumourevo = is_driver)
  }) %>% 
    dplyr::bind_rows()
}

merge_drivers = function(all_dr_process, all_drivers_tumourevo) {
  all_drivers = dplyr::full_join(all_dr_process, all_drivers_tumourevo, 
                                 suffix = c('_process', '_tumourevo'), 
                                 by = join_by('code' == 'driver_label', 
                                              'sample' == 'sample')) %>% 
    dplyr::rename(driver_label = code) %>% 
    dplyr::mutate(is_driver_process = ifelse((VAF_process == 0 | is.na(VAF_process)), FALSE, is_driver_process), 
                  is_driver_tumourevo = ifelse(VAF_tumourevo == 0, FALSE, is_driver_tumourevo)) %>% 
    dplyr::mutate(driver_class = 
                    dplyr::case_when(
                      (is_driver_process == FALSE & 
                         is_driver_tumourevo == FALSE) ~ 'Not a driver in sample', 
                      (is_driver_process == TRUE &
                         is_driver_tumourevo == FALSE) ~ 'Process True - Tumourevo False', 
                      (is_driver_process == FALSE &
                         is_driver_tumourevo == TRUE) ~ 'Process False - Tumourevo True',
                      (is_driver_process == TRUE &
                         is_driver_tumourevo == TRUE) ~ 'Process True - Tumourevo True', 
                    )) 
  
  return(all_drivers)
  
}

colors = setNames(nm = c("Process False - Tumourevo True", 
                         "Process True - Tumourevo False", 
                         "Process True - Tumourevo True", 
                         'Not a driver in sample'),
                  object = c('goldenrod', 
                             'indianred4', 
                             'forestgreen', 
                             '#393E46'))


plot_drivers = function(x, 
                        colors) {
  
  plt = x %>% 
    ggplot2::ggplot(
      ggplot2::aes(
        x = driver_label, 
        y = sample, 
        fill = driver_class, alpha = VAF_tumourevo
      )) + 
    ggplot2::geom_tile(height = 0.8, width = 0.8) +
    ggplot2::theme_light() +
    ggplot2::scale_fill_manual(values = colors) + 
    # ggplot2::facet_grid(vars(sample), scales = 'free') + 
    ggplot2::theme(legend.position = 'bottom', 
                   axis.text.x = ggplot2::element_text(angle = 90, 
                                                       hjust = 1, 
                                                       vjust = 0.5)) +
    ggplot2::guides(fill = ggplot2::guide_legend(title.position = "top",
                                                 title.hjust =0.5,
                                                 nrow = 2))
  
  return(plt)
  
}


discovery_rate_plot = function(x, colors) {
  x %>% 
    dplyr::filter(driver_class != 'Not a driver in sample') %>%
    ggplot2::ggplot(ggplot2::aes(
      sample,
      fill = driver_class
    )) + 
  ggplot2::geom_histogram(stat = 'count', position = 'dodge') +
  ggplot2::theme_light() + 
  ggplot2::scale_fill_manual(values = colors)
  
}
  

