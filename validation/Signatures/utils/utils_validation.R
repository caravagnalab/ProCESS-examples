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







