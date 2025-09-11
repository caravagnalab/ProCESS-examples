library(caret)

method_colors = c(
  "ProCESS" = "gray80",
  "mutect2" = "lightsteelblue",
  "mutect2 (all)" = "steelblue4",  # darker steelblue
  "strelka" = "coral",
  "strelka (all)" = "coral4",     # darker coral
  "freebayes" = "#8FBC8B",
  "freebayes (all)" = "#228B22"   # darker green (ForestGreen)
)

# Function to compute confusion matrix and performance metrics
compute_metrics <- function(actual, predicted) {
  cm <- table(Actual = actual, Predicted = predicted)
  confusion_matrix <- caret::confusionMatrix(as.factor(predicted), as.factor(actual), positive = "1")
  
  metrics <- dplyr::tibble(
    Accuracy = confusion_matrix$overall["Accuracy"],
    Sensitivity = confusion_matrix$byClass["Sensitivity"],
    Precision = confusion_matrix$byClass["Precision"],
    Recall = confusion_matrix$byClass["Recall"],
    F1_Score = confusion_matrix$byClass["F1"]
  )
  
  return(metrics)
}

plot_venn_diagram = function(merged_df, caller_name) {
  id_muts_races = merged_df %>% dplyr::filter(positive_truth) %>% dplyr::pull(mutationID)
  id_muts_caller = merged_df %>% dplyr::filter(positive_call) %>% dplyr::pull(mutationID)
  x = list("races"=id_muts_races, caller_name=id_muts_caller)
  names(x) = c("races", caller_name)
  
  ggVennDiagram::ggVennDiagram(x) +
    ggplot2::scale_fill_gradient2(low = "white", high = "#4981BF", mid = "white", midpoint=0) +
    ggplot2::coord_flip() +
    ggplot2::theme(legend.position = "none")
}

# Function to plot the confusion matrix
plot_confusion_matrix <- function(actual, predicted) {
  # Create the confusion matrix
  cm <- table(Actual = actual, Predicted = predicted)
  
  # Convert to data frame for ggplot
  cm_df <- as.data.frame(cm)
  
  # Calculate global percentages (out of total observations)
  total <- sum(cm)
  cm_df$Percentage <- (cm_df$Freq / total) * 100
  
  # Create labels with both count and percentage
  cm_df$Label <- paste0(cm_df$Freq, "\n(", round(cm_df$Percentage, 1), "%)")
  
  # Add a column to identify diagonal vs. off-diagonal elements
  cm_df$Diagonal <- ifelse(cm_df$Actual == cm_df$Predicted, "Correct", "Incorrect")
  
  # Plot
  ggplot2::ggplot(cm_df, ggplot2::aes(x = Predicted, y = Actual, fill = Diagonal)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = Label), color = "black", size = 4) +
    ggplot2::scale_fill_manual(values = c("Correct" = alpha("#2E8B57", .9), "Incorrect" = alpha("indianred", .9))) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Confusion Matrix", x = "Predicted", y = "Actual")
}

plot_scatter_with_corr <- function(data, x_var, y_var, col_var = "FILTER", title = NULL) {
  data[[x_var]][is.na(data[[x_var]])] <- 0
  data[[y_var]][is.na(data[[y_var]])] <- 0
  
  if (nrow(data) < 3) {
    p <- ggplot2::ggplot(data, ggplot2::aes_string(x = x_var, y = y_var, col = col_var)) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
      ggplot2::theme_bw() +
      ggplot2::labs(x = x_var, y = y_var)
    return(p)
  }
  
  corr_test <- stats::cor.test(data[[x_var]], data[[y_var]], use = "complete.obs")
  corr_coef <- corr_test$estimate
  p_val <- corr_test$p.value
  
  if (is.na(p_val)) {
    p_text <- "NA"
    corr_coef <- NA
  } else {
    p_text <- format.pval(p_val, digits = 1)
  }
  
  corr_text <- paste0("r = ", round(corr_coef, 2), ", ", p_text)
  
  p <- ggplot2::ggplot(data, ggplot2::aes_string(x = x_var, y = y_var, col = col_var)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    ggplot2::annotate("label", x = Inf, y = -Inf, label = corr_text,
                      hjust = 1.1, vjust = -0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = x_var, y = y_var)
  p
}

plot_filter_distribution <- function(merged_df, colors, log_scale = TRUE) {
  p <- merged_df %>%
    dplyr::group_by(FILTER) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    tidyr::drop_na() %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(FILTER, +n), y = n, fill = FILTER)) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::geom_col() +
    ggplot2::theme_bw() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "FILTER", y = "count") +
    ggplot2::theme(legend.position = 'none')
  
  if (log_scale) {
    p <- p + ggplot2::scale_y_continuous(trans = "log10")
  }
  
  p
}

plot_precision_recall <- function(df) {
  pr <- df %>%
    dplyr::arrange(dplyr::desc(VAF_caller)) %>%
    dplyr::mutate(
      TP = cumsum(true_positive),
      FP = cumsum(false_positive),
      FN = sum(true_positive) - TP,
      precision = TP / (TP + FP),
      recall = TP / (TP + FN)
    )
  
  vaf_threshold_point <- df %>%
    dplyr::filter(!is.na(VAF_caller)) %>%
    dplyr::slice(which.min(abs(VAF_caller - 0.1)))
  
  threshold_recall <- pr$recall[which.min(abs(df$VAF_caller - 0.1))]
  threshold_precision <- pr$precision[which.min(abs(df$VAF_caller - 0.1))]
  
  ggplot2::ggplot(pr, ggplot2::aes(x = recall, y = precision)) +
    ggplot2::geom_line(color = "blue") +
    ggplot2::geom_point(ggplot2::aes(x = threshold_recall, y = threshold_precision),
                        color = "red", size = 3) +
    ggplot2::annotate("text", x = threshold_recall, y = threshold_precision,
                      label = "VAF = 0.1", vjust = -1, color = "red", size = 4) +
    ggplot2::labs(title = "Precision-Recall Curve", x = "Recall", y = "Precision") +
    ggplot2::theme_bw()
}

plot_roc_curve <- function(df, vaf_thresh = 0.1) {
  roc_data <- df %>%
    dplyr::arrange(dplyr::desc(VAF_caller)) %>%
    dplyr::mutate(
      TP = cumsum(true_positive),
      FP = cumsum(false_positive),
      TN = sum(true_negative),
      FN = sum(false_negative),
      TPR = TP / (TP + FN),
      FPR = FP / (FP + TN)
    )
  
  vaf_threshold_point <- df %>%
    dplyr::filter(!is.na(VAF_caller)) %>%
    dplyr::slice(which.min(abs(VAF_caller - vaf_thresh)))
  
  threshold_tpr <- roc_data$TPR[which.min(abs(df$VAF_caller - vaf_thresh))]
  threshold_fpr <- roc_data$FPR[which.min(abs(df$VAF_caller - vaf_thresh))]
  
  ggplot2::ggplot(roc_data, ggplot2::aes(x = FPR, y = TPR)) +
    ggplot2::geom_line(color = "blue") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    ggplot2::geom_point(ggplot2::aes(x = threshold_fpr, y = threshold_tpr),
                        color = "red", size = 3) +
    ggplot2::annotate("text", x = threshold_fpr, y = threshold_tpr,
                      label = paste0("VAF = ", vaf_thresh), vjust = -1, color = "red", size = 4) +
    ggplot2::labs(title = "ROC Curve", x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)") +
    ggplot2::theme_bw()
}

# Function to plot Histogram of VAF Differences

plot_vaf_difference <- function(df) {
  df <- df %>%
    dplyr::mutate(VAF_diff = VAF_caller - VAF_truth)
  
  ggplot2::ggplot(df, mapping = ggplot2::aes(x = FILTER, y = VAF_diff)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(title = "Boxplot of VAF differences", 
                  y = "VAF Difference (Caller - Truth)", x = "") +
    ggplot2::theme_bw()
}

plot_cov_difference <- function(df) {
  df <- df %>%
    dplyr::mutate(DP_diff = DP_caller - DP_truth)
  
  ggplot2::ggplot(df, mapping = ggplot2::aes(x = FILTER, y = DP_diff)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(title = "Boxplot of DP differences", 
                  y = "DP Difference (Caller - Truth)", x = "") +
    ggplot2::theme_bw()
}

plot_depth_vs_performance <- function(df) {
  df %>%
    dplyr::group_by(DP_caller) %>%
    dplyr::summarize(
      precision = mean(true_positive / (true_positive + false_positive), na.rm = TRUE),
      recall = mean(true_positive / (true_positive + false_negative), na.rm = TRUE)
    ) %>%
    tidyr::pivot_longer(cols = c(precision, recall), names_to = "metric", values_to = "value") %>%
    ggplot2::ggplot(ggplot2::aes(x = DP_caller, y = value, color = metric)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Read Depth vs. Precision/Recall", 
                  x = "Read Depth", y = "Performance Metric") +
    ggplot2::theme_bw()
}

plot_chromosomal_errors <- function(df) {
  df %>%
    dplyr::group_by(chr_caller) %>%
    dplyr::summarize(
      false_positive = sum(false_positive), 
      false_negative = sum(false_negative)
    ) %>%
    tidyr::pivot_longer(cols = c(false_positive, false_negative), 
                        names_to = "error_type", values_to = "count") %>%
    ggplot2::ggplot(ggplot2::aes(x = chr_caller, y = count, fill = error_type)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(title = "Chromosomal Distribution of Errors", 
                  x = "Chromosome", y = "Count") +
    ggplot2::theme_bw()
}

plot_depth_comparison <- function(df) {
  df <- df %>%
    dplyr::filter(!is.na(DP_caller) & !is.na(DP_truth))
  
  ggplot2::ggplot(df, ggplot2::aes(x = DP_truth, y = DP_caller)) +
    ggplot2::geom_point(alpha = 0.5, color = "blue") +
    ggplot2::geom_abline(slope = 1, intercept = 0, 
                         linetype = "dashed", color = "red") +
    ggplot2::labs(title = "Comparison of Read Depth (Caller vs. Ground Truth)", 
                  x = "Ground Truth Depth", y = "Caller Depth") +
    ggplot2::theme_bw()
}

plot_races_coverage <- function(seq_res_long) {
  seq_res_long %>% 
    ggplot2::ggplot(ggplot2::aes(x = DP)) +
    ggplot2::geom_histogram(binwidth = 1) +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(ggplot2::aes(xintercept = stats::median(DP)), color = "indianred") +
    ggplot2::labs(
      title = "rRACES coverage",
      subtitle = paste0(nrow(seq_res_long), " mutations; Median coverage = ", 
                        stats::median(seq_res_long$DP)),
      color = ""
    )
}

plot_races_vaf <- function(seq_res_long) {
  seq_res_long %>%
    ggplot2::ggplot(ggplot2::aes(x = VAF)) +
    ggplot2::xlim(c(0, 1)) +
    ggplot2::geom_histogram(binwidth = 0.01) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "rRACES VAF",
      subtitle = paste0(nrow(seq_res_long), " mutations"),
      color = ""
    )
}

plot_caller_coverage <- function(caller_res_pass, sample_info) {
  caller_res_pass %>%
    ggplot2::ggplot(ggplot2::aes(x = DP)) +
    ggplot2::geom_histogram(binwidth = 1) +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(ggplot2::aes(xintercept = stats::median(DP)), color = "indianred") +
    ggplot2::labs(
      title = paste0(sample_info$caller_name, " coverage"),
      subtitle = paste0(nrow(caller_res_pass), 
                        " PASS mutations; Median coverage = ", 
                        stats::median(caller_res_pass$DP)),
      color = ""
    )
}

# Function to plot caller VAF histogram
plot_caller_vaf <- function(caller_res_pass, sample_info) {
  caller_res_pass %>% 
    ggplot2::ggplot(mapping = ggplot2::aes(x = VAF)) +
    ggplot2::geom_histogram(binwidth = 0.01) +
    ggplot2::theme_bw() +
    ggplot2::xlim(c(0, 1)) +
    ggplot2::labs(
      title = paste0(sample_info$caller_name, " VAF"),
      subtitle = paste0(base::nrow(caller_res_pass), " PASS mutations"),
      color = ""
    )
}

plot_false_negative_vaf_dist <- function(merged_df0) {
  merged_df0 <- merged_df0 %>%
    dplyr::filter(false_negative)
  
  median_VAF0 <- stats::median(merged_df0$VAF_truth)
  median_CCF = stats::median(merged_df0$CCF)
  
  merged_df0 %>%
    dplyr::filter(false_negative, CCF > 0) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = CCF)) +
    ggplot2::geom_histogram(bins = 100) +
    ggplot2::geom_vline(xintercept = median_CCF, col = "indianred") +
    #ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "CCF process", 
      y = "Count", 
      title = "False negatives' Cancer Cell Fraction"
    )
  
  # merged_df0 %>%
  #   dplyr::filter(false_negative) %>%
  #   ggplot2::ggplot(mapping = ggplot2::aes(x = VAF_truth)) +
  #   ggplot2::geom_histogram(bins = 100) +
  #   ggplot2::geom_vline(xintercept = median_VAF0, col = "indianred") +
  #   ggplot2::scale_x_continuous(trans = "log10") +
  #   ggplot2::theme_bw() +
  #   ggplot2::labs(
  #     x = "VAF races", 
  #     y = "Count", 
  #     title = "False negatives' Variant allele frequency"
  #   )
}

plot_vaf_scatter_all <- function(merged_df, colors, sample_info) {
  merged_df <- merged_df %>% dplyr::filter(positive_truth | positive_call)
  
  p <- plot_scatter_with_corr(merged_df, "VAF_truth", "VAF_caller") +
    ggplot2::ggtitle("VAF correlation", subtitle = paste0(base::nrow(merged_df), " total mutations")) +
    ggplot2::xlim(c(0, 1)) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      x = "VAF races",
      y = paste0("VAF ", sample_info$caller_name)
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1)))
  
  p_with_marginal <- ggExtra::ggMarginal(p, type = "boxplot", groupFill = TRUE, groupColour = TRUE)
  ggplotify::as.ggplot(p_with_marginal)
}

plot_dp_scatter_all <- function(merged_df, colors, sample_info) {
  merged_df <- merged_df %>% dplyr::filter(positive_truth | positive_call)
  
  p <- plot_scatter_with_corr(merged_df, "DP_truth", "DP_caller") +
    ggplot2::ggtitle("DP correlation", subtitle = paste0(base::nrow(merged_df), " total mutations")) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      x = "DP races",
      y = paste0("DP ", sample_info$caller_name)
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1)))
  
  p_with_marginal <- ggExtra::ggMarginal(p, type = "boxplot", groupFill = TRUE, groupColour = TRUE)
  ggplotify::as.ggplot(p_with_marginal)
}

plot_vaf_scatter_pass <- function(merged_df, colors, sample_info) {
  merged_df <- merged_df %>% dplyr::filter(positive_truth | positive_call)
  
  p <- plot_scatter_with_corr(merged_df, "VAF_truth", "VAF_caller") +
    ggplot2::ggtitle("VAF correlation", subtitle = paste0(base::nrow(merged_df), " total mutations")) +
    ggplot2::xlim(c(0, 1)) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      x = "VAF races",
      y = paste0("VAF ", sample_info$caller_name)
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1)))
  
  p_with_marginal <- ggExtra::ggMarginal(p, type = "boxplot", groupFill = TRUE, groupColour = TRUE)
  ggplotify::as.ggplot(p_with_marginal)
}

plot_dp_scatter_pass <- function(merged_df, colors, sample_info) {
  merged_df <- merged_df %>% dplyr::filter(positive_truth | positive_call)
  
  p <- plot_scatter_with_corr(merged_df, "DP_truth", "DP_caller") +
    ggplot2::ggtitle("DP correlation", subtitle = paste0(base::nrow(merged_df), " total mutations")) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      x = "DP races",
      y = paste0("DP ", sample_info$caller_name)
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1)))
  
  p_with_marginal <- ggExtra::ggMarginal(p, type = "boxplot", groupFill = TRUE, groupColour = TRUE)
  ggplotify::as.ggplot(p_with_marginal)
}

plot_metrics_comparison <- function(metrics) {
  metrics %>% 
    ggplot2::ggplot(mapping = ggplot2::aes(x = name, y = value, fill = Mutations)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::theme_bw() +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::labs(title = "Performance Metrics Comparison", 
                  x = "Metric", y = "Value")
}

compute_metrics_from_vectors <- function(y_true, y_pred) {
  cm <- base::table(Actual = y_true, Predicted = y_pred)
  
  TP <- base::sum(y_true == 1 & y_pred == 1)
  FP <- base::sum(y_true == 0 & y_pred == 1)
  TN <- base::sum(y_true == 0 & y_pred == 0)
  FN <- base::sum(y_true == 1 & y_pred == 0)
  
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  sensitivity <- TP / (TP + FN)
  precision <- TP / (TP + FP)
  recall <- sensitivity
  f1_score <- 2 * (precision * recall) / (precision + recall)
  false_positive_rate <- FP / (FP + TN)
  
  dplyr::tibble(
    Accuracy = accuracy,
    Precision = precision,
    Recall = recall,
    F1_Score = f1_score,
    FPR = false_positive_rate
  )
}

plot_metric_over_VAF_threshold <- function(seq_res_long, caller_res, only_pass, VAF_spectrum = seq(0, 0.1, by = 0.005)) {
  if (only_pass) {
    caller_res <- caller_res %>% dplyr::filter(FILTER == "PASS")
  }
  
  dfm <- lapply(VAF_spectrum, function(min_vaf) {
    merged_df <- merge_datasets(caller_res, seq_res_long, min_vaf)
    y_true <- base::as.numeric(base::factor((merged_df$VAF_truth > min_vaf), levels = c(FALSE, TRUE))) - 1
    y_pred <- base::as.numeric(base::factor((merged_df$VAF_caller > min_vaf), levels = c(FALSE, TRUE))) - 1
    metrics <- compute_metrics_from_vectors(y_true, y_pred)
    dplyr::tibble(min_vaf = min_vaf, value = base::as.numeric(metrics), metric = base::colnames(metrics))
  }) %>% base::do.call("rbind", .)
  
  metric_colors <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"
  )
  
  dfm %>%
    ggplot2::ggplot(ggplot2::aes(x = min_vaf, y = value, col = metric)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "VAF threshold", y = "Value") +
    ggplot2::scale_color_manual(values = metric_colors) +
    ggplot2::ylim(c(0, 1))
}

plot_metric_over_CCF_threshold <- function(seq_res_long, caller_res, only_pass, min_vaf_caller, CCF_spectrum = seq(0, 0.1, by = 0.005)) {
  if (only_pass) {
    caller_res <- caller_res %>% dplyr::filter(FILTER == "PASS")
  }
  
  dfm <- lapply(CCF_spectrum, function(min_ccf) {
    merged_df <- merge_datasets(caller_res, seq_res_long, min_ccf = min_ccf, min_vaf_caller = min_vaf_caller)
    y_true <- base::as.numeric(base::factor((merged_df$CCF > min_ccf), levels = c(FALSE, TRUE))) - 1
    y_pred <- base::as.numeric(base::factor((merged_df$VAF_caller > min_vaf_caller), levels = c(FALSE, TRUE))) - 1
    metrics <- compute_metrics_from_vectors(y_true, y_pred)
    dplyr::tibble(min_ccf = min_ccf, value = base::as.numeric(metrics), metric = base::colnames(metrics))
  }) %>% base::do.call("rbind", .)
  
  metric_colors <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"
  )
  
  dfm %>%
    ggplot2::ggplot(ggplot2::aes(x = min_ccf, y = value, col = metric)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "CCF threshold", y = "Value", color = "Metric") +
    ggplot2::scale_color_manual(values = metric_colors) +
    ggplot2::ylim(c(0, 1))
}

plot_flow_of_calls <- function(merged_df) {
  merged_df$FILTER[base::is.na(merged_df$FILTER)] <- "Not called"
  
  df <- merged_df %>%
    dplyr::select(FILTER, positive_truth, positive_call) %>%
    dplyr::filter(positive_truth | positive_call) %>%
    dplyr::mutate(
      Pre_filtering = dplyr::case_when(
        positive_truth & positive_call ~ "shared",
        positive_truth ~ "races private",
        positive_call ~ "caller private",
        .default = "N/A"
      ),
      Post_filtering = dplyr::case_when(
        (FILTER != "PASS") & (!positive_truth) ~ "Correclty filtered out",
        (FILTER != "PASS") & (positive_truth) ~ "races private",
        positive_truth & positive_call ~ "shared",
        positive_truth ~ "races private",
        positive_call ~ "caller private",
        .default = "N/A"
      )
    ) %>%
    dplyr::select(FILTER, Pre_filtering, Post_filtering) %>%
    dplyr::group_by(FILTER, Pre_filtering, Post_filtering) %>%
    dplyr::summarise(Freq = dplyr::n(), .groups = "drop")
  
  df %>%
    ggplot2::ggplot(ggplot2::aes(y = Freq, axis1 = Pre_filtering, axis2 = Post_filtering)) +
    ggalluvial::geom_flow(ggplot2::aes(fill = FILTER), curve_type = "quintic") +
    ggalluvial::geom_stratum() +
    ggplot2::geom_text(stat = "stratum", ggplot2::aes(label = ggplot2::after_stat(stratum))) +
    ggplot2::scale_x_discrete(limits = c("Pre_filtering", "Post_filtering"), expand = c(.05, .05)) +
    ggplot2::scale_fill_brewer(type = "qual", palette = "Set1") +
    ggplot2::theme_bw()
  
  df <- merged_df %>%
    dplyr::select(FILTER, positive_truth, positive_call) %>%
    dplyr::filter(positive_truth | positive_call) %>%
    dplyr::mutate(
      Pre_filtering = dplyr::case_when(
        positive_truth & positive_call ~ "TP",
        positive_truth ~ "FN",
        positive_call ~ "FP",
        .default = "TN"
      ),
      Post_filtering = dplyr::case_when(
        (FILTER != "PASS") & (!positive_truth) ~ "TN",
        (FILTER != "PASS") & (positive_truth) ~ "FN",
        positive_truth & positive_call ~ "TP",
        positive_truth ~ "FN",
        positive_call ~ "FP",
        .default = "N/A"
      )
    ) %>%
    dplyr::select(FILTER, Pre_filtering, Post_filtering) %>%
    dplyr::group_by(FILTER, Pre_filtering, Post_filtering) %>%
    dplyr::summarise(Freq = dplyr::n(), .groups = "drop")
  
  df %>%
    ggplot2::ggplot(ggplot2::aes(y = Freq, axis1 = Pre_filtering, axis2 = Post_filtering)) +
    ggalluvial::geom_flow(ggplot2::aes(fill = FILTER), curve_type = "quintic") +
    ggalluvial::geom_stratum() +
    ggplot2::geom_text(stat = "stratum", ggplot2::aes(label = ggplot2::after_stat(stratum))) +
    ggplot2::scale_x_discrete(limits = c("Pre_filtering", "Post_filtering"), expand = c(.05, .05)) +
    ggplot2::scale_fill_brewer(type = "qual", palette = "Set1") +
    ggplot2::theme_bw()
}

plot_ecdf_comparison <- function(vec1, vec2, label1, label2, x_label) {
  df <- dplyr::bind_rows(
    data.frame(value = vec1, group = label1),
    data.frame(value = vec2, group = label2)
  )
  
  ks_result <- stats::ks.test(vec1, vec2)
  
  ggplot2::ggplot(df, ggplot2::aes(x = value, color = group)) +
    ggplot2::stat_ecdf() +
    ggplot2::labs(
      title = "ECDF Comparison",
      subtitle = sprintf("KS test p-value: %.4g", ks_result$p.value),
      x = x_label,
      y = "ECDF",
      color = ""
    ) +
    ggplot2::theme_minimal()
}

plot_density_comparison <- function(vec1, vec2, label1, label2, x_label) {
  df <- dplyr::bind_rows(
    data.frame(value = vec1, group = label1),
    data.frame(value = vec2, group = label2)
  )
  
  ks_result <- stats::ks.test(vec1, vec2)
  
  n1 <- base::format(base::length(vec1), big.mark = "'")
  n2 <- base::format(base::length(vec2), big.mark = "'")
  
  subtitle_text <- sprintf(
    "%s (n = %s), %s (n = %s)\nKS test p-value: %.4g",
    label1, n1, label2, n2, ks_result$p.value
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = value, fill = group)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::labs(
      title = paste0(x_label, " comparison"),
      subtitle = subtitle_text,
      x = x_label,
      y = "Density",
      fill = ""
    ) +
    ggplot2::theme_bw()
}

plot_density_comparison_multi <- function(vec_list, labels, colors, x_label) {
  if (length(vec_list) != length(labels) || length(labels) != length(colors)) {
    stop("vec_list, labels, and colors must have the same length.")
  }
  
  df <- purrr::map2_dfr(vec_list, labels, ~ data.frame(value = .x, group = .y))
  
  ks_tests <- purrr::map2(vec_list[-1], labels[-1], ~ stats::ks.test(vec_list[[1]], .x))
  p_values <- purrr::map_chr(ks_tests, ~ sprintf("%.4g", .x$p.value))
  comp_labels <- paste0(labels[1], " vs ", labels[-1], ": p = ", p_values)
  
  ns <- purrr::map_chr(vec_list, ~ format(length(.x), big.mark = "'"))
  group_info <- paste0(labels, " (n = ", ns, ")")
  subtitle_text <- paste(paste(group_info, collapse = ", "))
  
  ggplot2::ggplot(df, ggplot2::aes(x = value, fill = group)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      title = paste0(x_label, " Comparison"),
      subtitle = subtitle_text,
      x = x_label,
      y = "Density",
      fill = ""
    ) +
    ggplot2::theme_bw()
}

plot_ecdf_comparison_multi <- function(vec_list, labels, colors, x_label) {
  if (length(vec_list) != length(labels) || length(labels) != length(colors)) {
    stop("vec_list, labels, and colors must have the same length.")
  }
  
  df <- purrr::map2_dfr(vec_list, labels, ~ data.frame(value = .x, group = .y))
  
  ks_tests <- purrr::map2(vec_list[-1], labels[-1], ~ stats::ks.test(vec_list[[1]], .x))
  p_values <- purrr::map_chr(ks_tests, ~ sprintf("%.4g", .x$p.value))
  ns <- purrr::map_chr(vec_list, ~ format(length(.x), big.mark = "'"))
  group_info <- paste0(labels, " (n = ", ns, ")")
  subtitle_text <- paste(paste(group_info, collapse = ", "))
  
  ggplot2::ggplot(df, ggplot2::aes(x = value, color = group)) +
    ggplot2::stat_ecdf(geom = "step", size = 1) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      title = paste0(x_label, " ECDF Comparison"),
      subtitle = subtitle_text,
      x = x_label,
      y = "ECDF",
      color = ""
    ) +
    ggplot2::theme_bw()
}

plot_precision_recall_vaf <- function(vaf_analysis_results, 
                                      title = "Precision and Recall Across VAF Bins",
                                      text_size = 3.5,
                                      point_size = 1.5,
                                      line_size = 1) {
  performance_data <- vaf_analysis_results$performance_table
  
  required_cols <- c("VAF_bin", "precision", "sensitivity")
  if (!all(required_cols %in% colnames(performance_data))) {
    stop("Required columns missing: ", 
         paste(setdiff(required_cols, colnames(performance_data)), collapse = ", "))
  }
  
  plot_data <- performance_data %>%
    dplyr::select(VAF_bin, precision, recall = sensitivity) %>%
    dplyr::filter(!is.na(VAF_bin)) %>%
    tidyr::pivot_longer(cols = c(precision, recall),
                        names_to = "metric", values_to = "value") %>%
    dplyr::mutate(
      metric = factor(metric, levels = c("precision", "recall"), 
                      labels = c("Precision", "Recall")),
      value_pct = round(value * 100, 1),
      label_text = paste0(value_pct, "%")
    )
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = VAF_bin, y = value, color = metric, group = metric)) +
    ggplot2::geom_line(size = line_size, alpha = 0.8) +
    ggplot2::geom_point(size = point_size, alpha = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = label_text), vjust = -0.5, size = text_size,
                       show.legend = FALSE, fontface = "bold") +
    ggplot2::scale_y_continuous(
      limits = c(0, 1.05),
      breaks = seq(0, 1, 0.25),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::scale_color_manual(values = c("Precision" = "#E31A1C", "Recall" = "#1F78B4"),
                                name = "Metric") +
    ggplot2::labs(
      title = title,
      x = "VAF Bin",
      y = "Performance",
      caption = paste0("VAF tolerance: ", vaf_analysis_results$vaf_tolerance_pct, "%; ",
                       "Min VAF threshold: ", vaf_analysis_results$min_vaf_threshold)
    ) +
    ggplot2::theme_bw()
  
  if (!is.null(vaf_analysis_results$overall_metrics)) {
    overall <- vaf_analysis_results$overall_metrics
    if (!is.na(overall$precision) && !is.na(overall$sensitivity)) {
      annotation_text <- paste0(
        "Overall: Precision = ", round(overall$precision * 100, 1), "%, ",
        "Recall = ", round(overall$sensitivity * 100, 1), "%"
      )
      p <- p + ggplot2::labs(subtitle = annotation_text)
    }
  }
  
  return(p)
}

plot_precision_recall_ccf <- function(analysis_results, 
                                      title = "Precision and Recall Across CCF Bins",
                                      text_size = 3.5,
                                      point_size = 1.5,
                                      line_size = 1) {
  performance_data <- analysis_results$performance_table
  
  required_cols <- c("CCF_bin", "precision", "sensitivity")
  if (!all(required_cols %in% colnames(performance_data))) {
    stop("Required columns missing: ", 
         paste(setdiff(required_cols, colnames(performance_data)), collapse = ", "))
  }
  
  plot_data <- performance_data %>%
    dplyr::select(CCF_bin, precision, recall = sensitivity) %>%
    dplyr::filter(!is.na(CCF_bin)) %>%
    tidyr::pivot_longer(cols = c(precision, recall),
                        names_to = "metric", values_to = "value") %>%
    dplyr::mutate(
      metric = factor(metric, levels = c("precision", "recall"), 
                      labels = c("Precision", "Recall")),
      value_pct = round(value * 100, 1),
      label_text = paste0(value_pct, "%")
    )
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = CCF_bin, y = value, color = metric, group = metric)) +
    ggplot2::geom_line(size = line_size, alpha = 0.8) +
    ggplot2::geom_point(size = point_size, alpha = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = label_text), vjust = -0.5, size = text_size,
                       show.legend = FALSE, fontface = "bold") +
    ggplot2::scale_y_continuous(
      limits = c(0, 1.05),
      breaks = seq(0, 1, 0.25),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::scale_color_manual(values = c("Precision" = "#E31A1C", "Recall" = "#1F78B4"),
                                name = "Metric") +
    ggplot2::labs(
      title = title,
      x = "CCF Bin",
      y = "Performance",
      caption = paste0("Tolerance: ", analysis_results$tolerance_pct, "%; ",
                       "Min CCF threshold: ", analysis_results$min_ccf_threshold)
    ) +
    ggplot2::theme_bw()
  

  if (!is.null(analysis_results$overall_metrics)) {
    overall <- analysis_results$overall_metrics
    if (!is.na(overall$precision) && !is.na(overall$sensitivity)) {
      annotation_text <- paste0(
        "Overall: Precision = ", round(overall$precision * 100, 1), "%, ",
        "Recall = ", round(overall$sensitivity * 100, 1), "%"
      )
      p <- p + ggplot2::labs(subtitle = annotation_text)
    }
  }
  
  return(p)
}


# Create combined sensitivity and FPR plot
plot_across_vaf <- function(results) {
  plot_data <- results$performance_table %>%
    dplyr::select(VAF_bin, sensitivity_pct, precision_pct) %>%
    tidyr::pivot_longer(
      cols = c(sensitivity_pct, precision_pct),
      names_to = "metric",
      values_to = "percentage"
    )
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = VAF_bin, y = percentage, fill = metric)) +
    ggplot2::geom_col(position = "dodge", alpha = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(!is.na(percentage), paste0(round(percentage, 1), "%"), "NA")),
      position = ggplot2::position_dodge(width = 0.9), vjust = -0.5, size = 3
    ) +
    ggplot2::labs(
      title = "Variant Caller Performance by VAF Range",
      x = "VAF Range (Truth)",
      y = "Percentage (%)",
      fill = "Metric"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::scale_fill_manual(
      values = c("sensitivity_pct" = "steelblue", "precision_pct" = "coral"),
      labels = c("Recall", "Precision")
    ) +
    ggplot2::ylim(0, max(plot_data$percentage, na.rm = TRUE) * 1.1) +
    ggplot2::scale_y_continuous(breaks = c(0, 25, 50, 75, 100))
}


# Create separate FPR plot with improved TN definition
plot_fpr_by_vaf <- function(results) {
  plot_data <- results$performance_table %>%
    dplyr::filter(!is.na(fpr_pct))
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = VAF_bin, y = fpr_pct)) +
    ggplot2::geom_col(fill = "coral", alpha = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(round(fpr_pct, 2), "%")), vjust = -0.5) +
    ggplot2::labs(
      title = "False Positive Rate by VAF Range",
      x = "VAF Range",
      y = "False Positive Rate (%)",
      subtitle = "True Negatives = Ground truth variants below minimum VAF threshold"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::ylim(0, max(plot_data$fpr_pct, na.rm = TRUE) * 1.1)
}


# Create specificity plot (complement of FPR)
plot_specificity_by_vaf <- function(results) {
  plot_data <- results$performance_table %>%
    dplyr::filter(!is.na(specificity_pct))
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = VAF_bin, y = specificity_pct)) +
    ggplot2::geom_col(fill = "lightgreen", alpha = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(round(specificity_pct, 1), "%")), vjust = -0.5) +
    ggplot2::labs(
      title = "Specificity by VAF Range",
      x = "VAF Range",
      y = "Specificity (%)",
      subtitle = "Ability to correctly identify true negatives (low VAF variants)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::ylim(0, 100)
}


# Create VAF accuracy plot
plot_vaf_accuracy <- function(results) {
  accuracy_data <- results$performance_table %>%
    dplyr::select(VAF_bin, mean_abs_error, median_abs_error, n_variants) %>%
    dplyr::filter(!is.na(mean_abs_error) & n_variants > 0) %>%
    tidyr::pivot_longer(
      cols = c(mean_abs_error, median_abs_error),
      names_to = "metric",
      values_to = "error"
    )
  
  ggplot2::ggplot(accuracy_data, ggplot2::aes(x = VAF_bin, y = error, fill = metric)) +
    ggplot2::geom_col(position = "dodge", alpha = 0.7) +
    ggplot2::labs(
      title = "VAF Calling Accuracy by VAF Range (True Positives Only)",
      x = "VAF Range (Truth)",
      y = "Absolute Error",
      fill = "Error Metric"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::scale_fill_manual(
      values = c("mean_abs_error" = "coral", "median_abs_error" = "lightblue"),
      labels = c("Mean Abs Error", "Median Abs Error")
    )
}

plot_multicaller_precision_recall <- function(performance_data, x_name,
                                                  title = "Precision and Recall Across VAF Bins",
                                                  text_size = 3.5,
                                                  point_size = 1.5,
                                                  line_size = 1) {

  required_cols <- c(x_name, "precision", "sensitivity", "caller")
  if (!all(required_cols %in% colnames(performance_data))) {
    stop("Required columns missing from performance table: ", 
         paste(setdiff(required_cols, colnames(performance_data)), collapse = ", "))
  }
  
  plot_data <- performance_data %>%
    dplyr::select(x_name, precision, recall = sensitivity, caller) %>%
    dplyr::filter(!is.na(x_name)) %>%
    tidyr::pivot_longer(cols = c(precision, recall),
                        names_to = "metric",
                        values_to = "value") %>%
    dplyr::mutate(
      metric = factor(metric, levels = c("precision", "recall"), labels = c("Precision", "Recall")),
      value_pct = round(value * 100, 1),
      label_text = paste0(value_pct, "%")
    )
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data[[x_name]], y = value, color = caller, linetype = metric,
                 group = interaction(caller, metric))
  ) +
    ggplot2::geom_line(size = line_size, alpha = 0.8) +
    ggplot2::geom_point(size = point_size, alpha = 0.9) +
    ggplot2::geom_text(
      ggplot2::aes(label = label_text),
      vjust = -0.5, size = text_size,
      show.legend = FALSE, fontface = "bold"
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1.05),
      breaks = seq(0, 1, 0.25),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Precision" = "solid", "Recall" = "dashed"),
      name = "Metric"
    ) +
    ggplot2::scale_color_manual(name = "Caller", values = method_colors) +
    ggplot2::labs(
      title = title,
      x = x_name,
      y = "Performance"
    ) +
    ggplot2::theme_bw()
}

plot_multicaller_rmse_vaf <- function(performance_data, 
                                      title = "RMSE Across VAF Bins",
                                      text_size = 3.5,
                                      point_size = 1.5,
                                      line_size = 1) {
  required_cols <- c("VAF_bin", "rmse", "caller")
  if (!all(required_cols %in% colnames(performance_data))) {
    stop("Required columns missing from performance table: ", 
         paste(setdiff(required_cols, colnames(performance_data)), collapse = ", "))
  }
  
  plot_data <- performance_data %>%
    dplyr::select(VAF_bin, rmse, caller) %>%
    dplyr::filter(!is.na(VAF_bin), !is.na(rmse)) %>%
    dplyr::mutate(
      rmse_rounded = round(rmse, 3),
      label_text = as.character(rmse_rounded)
    )
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = VAF_bin, y = rmse, color = caller, group = caller)) +
    ggplot2::geom_line(size = line_size, alpha = 0.8) +
    ggplot2::geom_point(size = point_size, alpha = 0.9) +
    # Optional: Uncomment if RMSE values should be labeled
    # ggplot2::geom_text(ggplot2::aes(label = label_text), vjust = -0.5,
    #                    size = text_size, show.legend = FALSE, fontface = "bold") +
    ggplot2::scale_y_continuous(limits = c(0, max(plot_data$rmse, na.rm = TRUE) * 1.1)) +
    ggplot2::scale_color_manual(name = "Caller", values = method_colors) +
    ggplot2::labs(
      title = title,
      x = "VAF Bin",
      y = "RMSE"
    ) +
    ggplot2::theme_bw()
}
