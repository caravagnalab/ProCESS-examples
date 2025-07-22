### Validation of exposures ###

compute_mse <- function(inferred_mat, ground_truth_mat) {
  inferred_mat <- as.matrix(inferred_mat)
  ground_truth_mat <- as.matrix(ground_truth_mat)

  # Align columns
  common_cols <- intersect(colnames(inferred_mat), colnames(ground_truth_mat))
  inferred_mat <- inferred_mat[, common_cols, drop = FALSE]
  ground_truth_mat <- ground_truth_mat[, common_cols, drop = FALSE]

  # Compute squared error matrix
  mse_matrix <- (inferred_mat - ground_truth_mat)^2

  # MSE per signature (column-wise mean)
  mse_per_signature <- colMeans(mse_matrix, na.rm = TRUE)
  names(mse_per_signature) <- common_cols  

  mse_overall <- mean(mse_matrix, na.rm = TRUE)

  return(list(
    per_signature = mse_per_signature,
    overall = mse_overall
  ))
}



cosine_similarity_exposures <- function(gt, tool) {
  # Ensure samples order
  common_samples <- intersect(rownames(gt), rownames(tool))
  gt <- gt[common_samples, , drop = FALSE]
  tool <- tool[common_samples, , drop = FALSE]

  # Normalize rows (sample vectors) to length 1
  norm_rows <- function(mat) mat / sqrt(rowSums(mat^2))
  gt_norm <- norm_rows(gt)
  tool_norm <- norm_rows(tool)

  # Element-wise rowwise dot product = cosine similarity per sample
  cos_sim <- rowSums(gt_norm * tool_norm)
  cos_sim[is.na(cos_sim)] <- 0
  cos_sim
}


mse_per_signature <- function(gt, tool) {
  # Ensure same samples & signatures order
  common_samples <- intersect(rownames(gt), rownames(tool))
  gt <- gt[common_samples, , drop = FALSE]
  tool <- tool[common_samples, , drop = FALSE]

  # Align columns (signatures)
  common_sigs <- intersect(colnames(gt), colnames(tool))
  gt <- gt[, common_sigs, drop = FALSE]
  tool <- tool[, common_sigs, drop = FALSE]

  # Calculate column-wise MSE
  mse_vec <- colMeans((gt - tool)^2, na.rm = TRUE)
  return(mse_vec)
}



compute_exposure_metrics <- function(aligned_exposures) {
  cosine_results <- list()
  mse_results <- list()
  cosine_idx <- 1
  mse_idx <- 1
  
  for (method in names(aligned_exposures)) {
    keys <- names(aligned_exposures[[method]])
    
    for (key in keys) {
      # Extract matrices
      gt_mat <- aligned_exposures[[method]][[key]]$gt
      tool_mat <- aligned_exposures[[method]][[key]]$tool
      
      # Skip if ground truth is empty
      if (nrow(gt_mat) == 0 || ncol(gt_mat) == 0) next
      
      parts <- str_split(key, "_", simplify = TRUE)
      SPN <- parts[1]
      Coverage <- as.numeric(parts[3])
      Purity <- as.numeric(parts[5])
      
      # COSINE SIMILARITY
      # Align sample rows
      common_samples <- intersect(rownames(gt_mat), rownames(tool_mat))
      if (length(common_samples) == 0) next
      
      gt_common <- as.matrix(gt_mat[common_samples, , drop = FALSE])
      tool_common <- as.matrix(tool_mat[common_samples, , drop = FALSE])
      
      # Compute cosine similarity
      sims <- cosine_similarity_exposures(gt_common, tool_common)
      
      cosine_results[[cosine_idx]] <- tibble(
        Sample = common_samples,
        SPN = SPN,
        Coverage = Coverage,
        Purity = Purity,
        Tool = method,
        CosineSimilarity = sims
      )
      cosine_idx <- cosine_idx + 1
      
      # Align signatures
      common_sigs <- intersect(colnames(gt_mat), colnames(tool_mat))
      if (length(common_sigs) == 0) next
      
      gt_sub <- as.matrix(gt_mat[, common_sigs, drop = FALSE])
      tool_sub <- as.matrix(tool_mat[, common_sigs, drop = FALSE])
      
      # Compute MSE
      mse_vec <- colMeans((gt_sub - tool_sub)^2, na.rm = TRUE)
      
      mse_results[[mse_idx]] <- tibble(
        SPN = SPN,
        Coverage = Coverage,
        Purity = Purity,
        Tool = method,
        Signature = names(mse_vec),
        MSE = mse_vec
      )
      mse_idx <- mse_idx + 1
    }
  }
  
  # Combine results into data frames
  list(
    cosine = bind_rows(cosine_results),
    mse = bind_rows(mse_results)
  )
}



