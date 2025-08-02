map_sparsesig_to_cosmic <- function(sparsesig_out, cosmic_path, threshold = 0.7) {
  # Load COSMIC reference
  cosmic_signatures <- read.delim(cosmic_path) %>% 
    tibble::column_to_rownames("Type") %>%
    as.matrix()
  
  # Extract SparseSignature de novo signatures and exposures
  de_novo_signatures <- t(sparsesig_out[["beta"]]) %>% as.matrix()
  de_novo_signatures <- de_novo_signatures[rownames(cosmic_signatures), ]
  de_novo_exposures <- sparsesig_out[["alpha"]]
  
  # Compute cosine similarity matrix
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
  
  # Convert best_matches into a proper named list
  best_matches_list <- vector("list", length = nrow(best_matches))
  names(best_matches_list) <- rownames(best_matches)
  for (i in seq_along(best_matches)) {
    best_matches_list[[i]] <- best_matches[[i]]
  }
  
  # Handle "Background" separately (assign to SBS5)
  background_row <- similarity_matrix["Background", , drop = FALSE]
  max_col <- which.max(background_row[1, ])
  match_value <- background_row[1, max_col, drop = FALSE]
  col_name <- colnames(background_row)[max_col]
  
  # Create a match structure consistent with best_matches_list
  background_match <- setNames(list(match_value), col_name)
  
  # Combine into full match list
  sim_matrix_all <- c(list(background_match, best_matches_list))
  sim_matrix_all <- c(sim_matrix_all[[1]], sim_matrix_all[[2]])
  
  # Build remapped exposure matrix
  cosmic_sigs <- names(sim_matrix_all)
  samples <- rownames(de_novo_exposures)

  remapped_exposures <- matrix(0, nrow = length(samples), ncol = length(cosmic_sigs),
                               dimnames = list(samples, cosmic_sigs))
  
  # Background handled separately (already assigned to SBS5)
  for (de_novo_sig in colnames(de_novo_exposures)) {
    if (de_novo_sig == "Background") {
      for (sample in samples) {
        remapped_exposures[sample, "SBS5"] <- remapped_exposures[sample, "SBS5"] +
          de_novo_exposures[sample, de_novo_sig]
      }
    } else if (de_novo_sig == colnames(best_matches)) {
      # Extract and normalize weights from best_matches
      sim_weights <- best_matches[, de_novo_sig]
      sim_weights <- sim_weights / sum(sim_weights, na.rm = TRUE)
      
      for (sample in samples) {
        exposure <- de_novo_exposures[sample, de_novo_sig]
        for (sig in names(sim_weights)) {
          remapped_exposures[sample, sig] <- remapped_exposures[sample, sig] +
            exposure * sim_weights[sig]
        }
      }
    } else {
      warning(paste("No mapping for de novo sig:", de_novo_sig))
    }
  }
  
  # Normalize exposures
  row_totals <- rowSums(remapped_exposures)
  remapped_exposures_prop <- remapped_exposures / row_totals
  
  return(remapped_exposures_prop)
}


