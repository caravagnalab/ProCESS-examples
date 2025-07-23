library(caret)


safe_cor_test <- function(x, y, min_points = 3) {
  # Remove NA values
  valid_idx <- !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
  x_clean <- x[valid_idx]
  y_clean <- y[valid_idx]
  
  # Check if we have enough data points
  if (length(x_clean) < min_points) {
    return(NA_real_)
  }
  
  # Check for zero variance (constant values)
  if (stats::var(x_clean) == 0 || stats::var(y_clean) == 0) {
    return(NA_real_)
  }
  
  # Try correlation test with error handling
  tryCatch({
    cor_result <- stats::cor.test(x_clean, y_clean)
    return(cor_result$estimate)
  }, error = function(e) {
    return(NA_real_)
  })
}

# Function to safely calculate RMSE
safe_rmse <- function(actual, predicted) {
  # Remove NA values
  valid_idx <- !is.na(actual) & !is.na(predicted) & is.finite(actual) & is.finite(predicted)
  actual_clean <- actual[valid_idx]
  predicted_clean <- predicted[valid_idx]
  
  # Check if we have data points
  if (length(actual_clean) == 0) {
    return(NA_real_)
  }
  
  # Calculate RMSE
  sqrt(mean((actual_clean - predicted_clean)^2))
}


plot_density_comparison_multi <- function(vec_list, labels, colors, x_label, n) {
  if (length(vec_list) != length(labels) || length(labels) != length(colors)) {
    stop("vec_list, labels, and colors must have the same length.")
  }
  
  df <- purrr::map2_dfr(vec_list, labels, ~ data.frame(value = .x, group = .y))
  
  ks_tests <- purrr::map2(vec_list[-1], labels[-1], ~ stats::ks.test(vec_list[[1]], .x))
  p_values <- purrr::map_chr(ks_tests, ~ sprintf("%.4g", .x$p.value))
  comp_labels <- paste0(labels[1], " vs ", labels[-1], ": p = ", p_values)
  
  ns <- purrr::map_chr(vec_list, ~ format(length(.x), big.mark = "'"))
  group_info <- paste0('sample ', n, ' mutations')
  
  ggplot2::ggplot(df, ggplot2::aes(x = value, color = group)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::scale_color_manual('', values = colors) +
    ggplot2::labs(
      title = paste0(x_label, " Comparison"),
      subtitle = group_info,
      x = x_label,
      y = "Density",
      fill = ""
    ) +
    ggplot2::theme_bw()
}

plot_ecdf_comparison_multi <- function(vec_list, labels, colors, x_label, n) {
  if (length(vec_list) != length(labels) || length(labels) != length(colors)) {
    stop("vec_list, labels, and colors must have the same length.")
  }
  
  df <- purrr::map2_dfr(vec_list, labels, ~ data.frame(value = .x, group = .y))
  
  ks_tests <- purrr::map2(vec_list[-1], labels[-1], ~ stats::ks.test(vec_list[[1]], .x))
  p_values <- purrr::map_chr(ks_tests, ~ sprintf("%.4g", .x$p.value))
  ns <- purrr::map_chr(vec_list, ~ format(length(.x), big.mark = ""))
  group_info <- paste0('sample ', n, ' mutations')
  
  ggplot2::ggplot(df, ggplot2::aes(x = value, color = group)) +
    ggplot2::stat_ecdf(geom = "step", size = 1) +
    ggplot2::scale_color_manual('', values = colors) +
    ggplot2::labs(
      title = paste0(x_label, " ECDF Comparison"),
      subtitle = group_info,
      x = x_label,
      y = "ECDF",
      color = ""
    ) +
    ggplot2::theme_bw()
}


merge_datasets <- function(snp_caller, ground_truth) {
  df <- snp_caller %>%
    full_join(ground_truth, by = c("mutationID"), suffix = c(".caller", ".races")) %>%
    mutate(
      positive_truth = !is.na(BAF.races) & BAF.races > 0,
      positive_call = !is.na(BAF.caller) & BAF.caller > 0,
      true_positive = positive_truth & positive_call,
      false_positive = !positive_truth & positive_call,
      false_negative = positive_truth & !positive_call,
      true_negative = !positive_truth & !positive_call
    ) %>% 
    dplyr::mutate(BAF.races = ifelse(is.na(BAF.races), 0, BAF.races)) %>% 
    dplyr::mutate(BAF.caller = ifelse(is.na(BAF.caller), 0, BAF.caller)) %>% 
    dplyr::mutate(DP.races = ifelse(is.na(DP.races), 0, DP.races)) %>% 
    dplyr::mutate(DP.caller = ifelse(is.na(DP.caller), 0, DP.caller))
  
  return(df)
}

compute_metrics <- function(actual, predicted) {
  cm <- table(Actual = actual, Predicted = predicted)
  confusion_matrix <- confusionMatrix(as.factor(predicted), as.factor(actual), positive = "1")
  
  metrics <- dplyr::tibble(
    Accuracy = confusion_matrix$overall["Accuracy"],
    Sensitivity = confusion_matrix$byClass["Sensitivity"],
    Precision = confusion_matrix$byClass["Precision"],
    Recall = confusion_matrix$byClass["Recall"],
    F1_Score = confusion_matrix$byClass["F1"]
  )
  
  return(metrics)
}
