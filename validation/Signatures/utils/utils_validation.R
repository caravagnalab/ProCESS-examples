### Validation of Signatures ###

evaluate_signatures <- function(true_signatures, predicted_signatures) {
  
  # Union of all signatures involved
  all_signatures <- union(true_signatures, predicted_signatures)
  
  # Binary vectors
  true_labels <- ifelse(all_signatures %in% true_signatures, 1, 0)
  pred_labels <- ifelse(all_signatures %in% predicted_signatures, 1, 0)
  
  # Confusion matrix
  cm <- table(True = true_labels, Predicted = pred_labels)
  print("Confusion Matrix:")
  print(cm)
  
  # Compute metrics
  TP <- sum(true_labels == 1 & pred_labels == 1)
  FP <- sum(true_labels == 0 & pred_labels == 1)
  FN <- sum(true_labels == 1 & pred_labels == 0)
  TN <- sum(true_labels == 0 & pred_labels == 0)
  
  accuracy <- (TP + TN) / (TP + FP + FN + TN)
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  recall <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  f1_score <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA
  }
  
  return(list(
    ConfusionMatrix = cm,
    Accuracy = accuracy,
    Precision = precision,
    Recall = recall,
    F1_Score = f1_score,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN
  ))
}


evaluate_all_combined <- function(ground_truth_list, predicted_list, threshold = 0) {
  results <- data.frame()

  for (spn in names(ground_truth_list)) {
    for (coverage in names(ground_truth_list[[spn]])) {
      for (purity in names(ground_truth_list[[spn]][[coverage]])) {

        gt_mat <- ground_truth_list[[spn]][[coverage]][[purity]]
        pred_mat <- predicted_list[[spn]][[coverage]][[purity]]

        if (is.null(gt_mat) || is.null(pred_mat)) next

        # Align samples and signatures
        common_samples <- intersect(rownames(gt_mat), rownames(pred_mat))
        gt_mat <- gt_mat[common_samples, , drop = FALSE]
        pred_mat <- pred_mat[common_samples, , drop = FALSE]

        all_sigs <- union(colnames(gt_mat), colnames(pred_mat))

        # Convert to data.frame for easier column handling
        gt_mat <- as.data.frame(gt_mat)
        pred_mat <- as.data.frame(pred_mat)

        # Add missing columns with zeros
        for (sig in all_sigs) {
          if (!(sig %in% colnames(gt_mat))) gt_mat[[sig]] <- 0
          if (!(sig %in% colnames(pred_mat))) pred_mat[[sig]] <- 0
        }

        # Reorder columns
        gt_mat <- gt_mat[, all_sigs, drop = FALSE]
        pred_mat <- pred_mat[, all_sigs, drop = FALSE]

        # Flatten matrices into vectors: presence across all samples and all signatures
        gt_presence <- as.vector(gt_mat > threshold)
        pred_presence <- as.vector(pred_mat > threshold)

        # Get indices where present (=1)
        true_sigs <- which(gt_presence)
        pred_sigs <- which(pred_presence)

        # Calculate metrics
        metrics <- evaluate_signatures(true_sigs, pred_sigs)

        # Collect results
        results <- rbind(results, data.frame(
          SPN = spn,
          Coverage = coverage,
          Purity = purity,
          Accuracy = metrics$Accuracy,
          Precision = metrics$Precision,
          Recall = metrics$Recall,
          F1_Score = metrics$F1_Score,
          TP = metrics$TP,
          FP = metrics$FP,
          FN = metrics$FN,
          TN = metrics$TN
        ))
      }
    }
  }
  rownames(results) <- NULL
  return(results)
}


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







