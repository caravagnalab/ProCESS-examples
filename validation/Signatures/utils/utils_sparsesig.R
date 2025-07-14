map_sparsesig_to_cosmic <- function(sparsesig_out, mut_counts, cosmic_path, threshold = 0.7) {
  
  cosmic_signatures <- read.delim(cosmic_path) %>% 
    tibble::column_to_rownames("Type") %>%
    as.matrix()
  
  # Extract SparseSignature de novo signatures and exposures
  de_novo_signatures <- t(sparsesig_out[["beta"]]) %>% as.matrix()
  de_novo_exposures <- sparsesig_out[["alpha"]]
  
  # Calculate cosine similarity matrix
  similarity_matrix <- MutationalPatterns::cos_sim_matrix(de_novo_signatures, cosmic_signatures)
  
  # Separate background
  similarity_matrix_noBackground <- similarity_matrix[rownames(similarity_matrix) != "Background", , drop = FALSE]
  
  # Find best COSMIC matches above threshold
  best_matches <- apply(similarity_matrix_noBackground, 1, function(similarities) {
    above_idx <- which(similarities >= threshold)
    if (length(above_idx) == 0) {
      return(NA)  
    } else {
      matches <- similarities[above_idx]
      matches <- sort(matches, decreasing = TRUE)  
      return(matches)
    }
  })
  
  best_matches_denovo <- t(best_matches) %>% as.data.frame()
  
  # Handle "Background" separately (assign to SBS5)
  background_row <- as.data.frame(similarity_matrix["Background", , drop = FALSE])
  max_col <- which.max(background_row[1, ])
  background_max <- background_row[, max_col, drop = FALSE]
  
  sim_matrix_all <- list(background_max, best_matches_denovo)
  names(sim_matrix_all) <- rownames(similarity_matrix)
  
  # Build remapped exposure matrix
  cosmic_sigs <- unique(unlist(lapply(sim_matrix_all, names)))
  samples <- rownames(de_novo_exposures)
  remapped_exposures <- matrix(0, nrow = length(samples), ncol = length(cosmic_sigs),
                               dimnames = list(samples, cosmic_sigs))
  
  for (de_novo_sig in colnames(de_novo_exposures)) {
    if (de_novo_sig == "Background") {
      for (sample in samples) {
        remapped_exposures[sample, "SBS5"] <- remapped_exposures[sample, "SBS5"] +
          de_novo_exposures[sample, de_novo_sig]
      }
    } else {
      sim_weights <- sim_matrix_all[[de_novo_sig]]
      
      if (!is.null(sim_weights) && length(sim_weights) > 0 && sum(sim_weights) > 0) {
        sim_weights <- as.numeric(sim_weights)
        names(sim_weights) <- names(sim_matrix_all[[de_novo_sig]])
        sim_weights <- sim_weights / sum(sim_weights)
        
        for (sample in samples) {
          exposure <- de_novo_exposures[sample, de_novo_sig]
          remapped_exposures[sample, names(sim_weights)] <- 
            remapped_exposures[sample, names(sim_weights)] +
            exposure * sim_weights
        }
      }
    }
  }
  
  # Normalize exposures
  row_totals <- rowSums(remapped_exposures)
  remapped_exposures_prop <- remapped_exposures / row_totals
  
  return(remapped_exposures_prop)
}


