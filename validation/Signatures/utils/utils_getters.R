### get and reshape ProCESS exposure data ###

samples_table <- function(snapshot, sample_forest) {
  sim <- ProCESS::recover_simulation(snapshot)
  info = sim$get_samples_info() ## requested from either the simulation recovery or as saved table
  
  nodes = sample_forest$get_nodes()
  clones = nodes %>% 
    dplyr::filter(!is.na(sample)) %>% 
    dplyr::group_by(sample, mutant) %>% 
    dplyr::pull(mutant) %>% 
    unique()
  clones_of_origin = nodes %>%
    dplyr::filter(!is.na(sample)) %>% 
    dplyr::group_by(sample, mutant) %>% 
    dplyr::count(mutant) %>% 
    tidyr::pivot_wider(values_from = n, names_from = mutant, values_fill = 0) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(Sample_Type = sum(c_across(clones) == 0)) %>% 
    dplyr::mutate(Sample_Type = ifelse(Sample_Type == (length(clones)-1), "Monoclonal", "Polyclonal")) %>% 
    dplyr::mutate(Total_Cells = sum(c_across(clones), na.rm = T)) %>% 
    dplyr::mutate(across(all_of(clones), ~ round(.x/Total_Cells,2), .names = "{.col} proportion"))
  
  samples_tb = dplyr::full_join(info, clones_of_origin, by = join_by("name" == "sample")) %>% 
    dplyr::select(!c("xmin","xmax","ymin","ymax","id","tumour_cells","tumour_cells_in_bbox")) %>% 
    dplyr::rename(Sample_ID=name) %>% 
    dplyr::rename(Samping_Time=time) %>% 
    dplyr::mutate(Samping_Time=round(Samping_Time,2)) %>% 
    dplyr::arrange(Sample_ID)
  
  return(samples_tb)
}

get_process_exposures <- function(spn,coverage,purity){
  samples <- get_sample_names(spn = spn)
  muts <- readRDS(get_mutations(spn = spn,coverage=coverage,purity = purity,type = "tumour"))
  muts_som <- muts %>% 
    filter(classes!="germinal") %>% 
    seq_to_long()
  exposure_sbs <- list()
  exposure_ids <- list()
  for (s in samples){
    exposure_sbs[[s]] <- muts_som %>% 
      filter(sample_name==s & VAF!=0) %>% 
      filter(!grepl("errors", causes)) %>% 
      filter(grepl("SBS", causes)) %>% 
      group_by(causes) %>% 
      summarise(count = n(), .groups = 'drop') %>%
      mutate(exposure = count / sum(count)) %>% 
      mutate(Sample_ID=s)
    
    exposure_ids[[s]] <- muts_som %>% 
      filter(sample_name==s & VAF!=0) %>% 
      filter(!grepl("errors", causes)) %>% 
      filter(grepl("ID", causes)) %>% 
      group_by(causes) %>% 
      summarise(count = n(), .groups = 'drop') %>%
      mutate(exposure = count / sum(count)) %>% 
      mutate(Sample_ID=s)
  }
  
  exposure_sbs_all <- do.call("rbind",exposure_sbs) %>% 
    select(Sample_ID, causes, exposure) %>%
    pivot_wider(
      names_from = causes,
      values_from = exposure
    )
  exposure_ids_all <- do.call("rbind",exposure_ids) %>% 
    select(Sample_ID, causes, exposure) %>%
    pivot_wider(
      names_from = causes,
      values_from = exposure
    )
  return(list("SBS"=exposure_sbs_all,
              "ID"=exposure_ids_all))
}


# get_exposures_by_context <- function(spn, base_path, context = "SBS") {
#   
#   phylo_forest <- get_phylo_forest(spn = spn, base_path = base_path)
#   sample_forest <- get_sample_forest(spn = spn, base_path = base_path)
#   snapshot_path <- file.path(base_path, spn, "process_old", spn)
#   
#   # Load phyloforest object 
#   phylo_forest <- ProCESS::load_phylogenetic_forest(phylo_forest)
#   # Load sample forest object
#   samples_forest <- ProCESS::load_samples_forest(sample_forest)
#   
#   # Generate sample-level data
#   samples_data <- samples_table(snapshot = snapshot_path, sample_forest = samples_forest)
#   
#   sample_ids <- samples_data %>%
#     dplyr::pull(Sample_ID)
#   
#   # Get exposures
#   exposure_prop <- phylo_forest$get_exposures()
#   
#   # Filter out time = 0 exposures if others exist
#   if (any(exposure_prop$time != 0)) {
#     exposure_prop_filtered <- exposure_prop %>% dplyr::filter(time != 0)
#   } else {
#     exposure_prop_filtered <- exposure_prop
#   }
#   
#   # Match sample IDs to timepoints
#   assign_sample_id <- function(t) {
#     if (t == 0) {
#       return(sample_ids)
#     } else {
#       match_idx <- which.min(abs(samples_data$Samping_Time - t))
#       return(samples_data$Sample_ID[match_idx])
#     }
#   }
#   
#   # Expand each exposure row to a sample ID
#   expanded_rows <- do.call(rbind, lapply(1:nrow(exposure_prop_filtered), function(i) {
#     t <- exposure_prop_filtered$time[i]
#     sids <- assign_sample_id(t)
#     do.call(rbind, lapply(sids, function(sid) {
#       data.frame(
#         Sample_ID = sid,
#         signature = exposure_prop_filtered$signature[i],
#         exposure = exposure_prop_filtered$exposure[i],
#         stringsAsFactors = FALSE
#       )
#     }))
#   }))
#   
#   # Aggregate and reshape
#   avg_exposure <- expanded_rows %>%
#     dplyr::group_by(Sample_ID, signature) %>%
#     dplyr::summarise(mean_exposure = mean(exposure), .groups = "drop")
#   
#   exposure_matrix <- avg_exposure %>%
#     tidyr::pivot_wider(
#       names_from = signature,
#       values_from = mean_exposure,
#       values_fill = list(mean_exposure = 0)
#     ) %>%
#     dplyr::arrange(Sample_ID)
#   
#   # Keep only signatures matching the context (e.g., "SBS", "ID", etc.)
#   exposure_filtered <- exposure_matrix %>%
#     dplyr::select(Sample_ID, tidyselect::matches(paste0("^", context)))
#   
#   return(exposure_filtered)
# }
# 

load_signature_data <- function(paths) {
  lapply(paths, function(p) {
    if (!file.exists(p)) {
      warning("File not found: ", p)
      return(NULL)
    }

    if (grepl("\\.rds$", p)) {
      return(readRDS(p))
    } else if (grepl("\\.txt$", p)) {
      return(read.delim(p, header = TRUE, sep = "\t"))
    } else {
      warning(paste("Unknown file type:", p))
      return(NULL)
    }
  })
}
