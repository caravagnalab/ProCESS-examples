get_drivers = function(process_seq_res, drivers_annotated) {
  
  # select process drivers
  process_drivers_ids = process_seq_res %>% 
    dplyr::filter(classes == 'driver') %>% 
    dplyr::mutate(chr = paste0('chr', chr)) %>% 
    dplyr::mutate(mut_id = paste0(chr, ':', chr_pos, '_', ref, '>', alt)) %>% 
    dplyr::pull(mut_id) %>% 
    unique
  
  # select tumourevo drivers
  tumourevo_drivers_ids = lapply(names(drivers_annotated), function(x) {
    drivers_annotated[[x]][[1]]$mutations %>% 
      dplyr::mutate(mut_id = paste0(chr, ':', from, '_', ref, '>', alt)) %>% 
      dplyr::filter(is_driver, VAF > 0) %>% 
      dplyr::mutate(sample = x) %>% 
      dplyr::select(mut_id, driver_label, sample, is_driver, VAF) %>% 
      dplyr::rename(is_driver_tumourevo = is_driver)
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::pull(mut_id) %>% 
    unique
  
  all_drivers = c(process_drivers_ids, tumourevo_drivers_ids) %>% unique
  
  # filter the mutation tables to get only mutation ids that are labelled as driver in one of the two
  process_drivers_all = process_seq_res %>% 
    dplyr::mutate(chr = paste0('chr', chr)) %>% 
    dplyr::mutate(mut_id = paste0(chr, ':', chr_pos, '_', ref, '>', alt)) %>% 
    dplyr::filter(mut_id %in% all_drivers)
  
  tumourevo_drivers_ids = lapply(names(drivers_annotated), function(x) {
    drivers_annotated[[x]][[1]]$mutations %>% 
      dplyr::mutate(mut_id = paste0(chr, ':', from, '_', ref, '>', alt)) %>% 
      dplyr::filter(mut_id %in% all_drivers) %>% 
      dplyr::mutate(sample = x) %>% 
      dplyr::select(mut_id, driver_label, sample, is_driver, VAF) %>% 
      dplyr::rename(is_driver_tumourevo = is_driver)
  }) %>% 
    dplyr::bind_rows()
  
  
}


