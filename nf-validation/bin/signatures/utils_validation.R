# Align columns for cosine/MSE
align_matrices <- function(inferred, simulated) {
  all_sigs <- union(colnames(inferred), colnames(simulated))
  
  for (sig in setdiff(all_sigs, colnames(inferred))) {
    inferred[[sig]] <- rep(0, nrow(inferred))
  }
  for (sig in setdiff(all_sigs, colnames(simulated))) {
    simulated[[sig]] <- rep(0, nrow(inferred))
  }
  
  # Fill missing columns with 0
  inferred_full <- inferred[, all_sigs, drop = FALSE]
  simulated_full <- simulated[, all_sigs, drop = FALSE]
  
  
  # Order columns the same way
  inferred_full <- inferred_full[, sort(colnames(inferred_full))]
  simulated_full <- simulated_full[, sort(colnames(simulated_full))]
  simulated_full[is.na(simulated_full)] <- 0
  inferred_full[is.na(inferred_full)] <- 0
  return(list(inferred = as.matrix(inferred_full),
              simulated = as.matrix(simulated_full)))
}

compute_cosine_mse <- function(inferred, simulated) {
  aligned <- align_matrices(inferred, simulated)
  inf_mat <- aligned$inferred
  sim_mat <- aligned$simulated
  
  # Compute per-row cosine similarity and MSE
  cosine <- mapply(function(x, y) cosine(x, y), 
                   split(inf_mat, row(inf_mat)), 
                   split(sim_mat, row(sim_mat)))
  mse <- rowMeans((inf_mat - sim_mat)^2)
  
  return(data.frame(
    sample = rownames(inf_mat),
    cosine = cosine,
    mse = mse
  ))
}


### Validation of exposures ###
per_sample_metrics <- function(inferred_df, simulated_df) {
  # Get union of all signatures
  all_sigs <- union(colnames(inferred_df), colnames(simulated_df))
  
  # Order columns the same
  inferred_df <- inferred_df[, sort(colnames(inferred_df))]
  simulated_df <- simulated_df[, sort(colnames(simulated_df))]
  
  # Per sample metrics
  results <- data.frame()
  
  for (sample in rownames(simulated_df)) {
    inferred_active <- names(inferred_df[which(inferred_df[sample, ] > 0)])
    simulated_active <- names(simulated_df[which(simulated_df[sample, ] > 0)])
    
    TP <- length(intersect(inferred_active, simulated_active))
    FP <- length(setdiff(inferred_active, simulated_active))
    FN <- length(setdiff(simulated_active, inferred_active))
    
    precision <-  TP / (TP + FP)
    recall <- TP / (TP + FN)
    f1 <- ifelse(!is.na(precision) & !is.na(recall) & (precision + recall) > 0,
                 2 * precision * recall / (precision + recall), NA)
    
    results <- rbind(results, data.frame(
      sample = sample,
      TP = TP,
      FP = FP,
      FN = FN,
      precision = precision,
      recall = recall,
      f1 = f1
    ))
  }
  
  return(results)
}

summary_stats <- function(df) {
  data.frame(
    metric = c("precision", "recall", "f1"),
    mean = c(mean(df$precision, na.rm = TRUE),
             mean(df$recall, na.rm = TRUE),
             mean(df$f1, na.rm = TRUE)),
    sd = c(sd(df$precision, na.rm = TRUE),
           sd(df$recall, na.rm = TRUE),
           sd(df$f1, na.rm = TRUE))
  )
}


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



