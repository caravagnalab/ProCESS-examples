require(patchwork)
require(RColorBrewer)
require(ggExtra)
require(ggVennDiagram)
require(ggalluvial)
require(ggplotify)

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

string_difference <- function(str1, str2) {
  str1_chars <- unlist(strsplit(str1, ""))
  str2_chars <- unlist(strsplit(str2, ""))
  
  for (char in str2_chars) {
    match_idx <- match(char, str1_chars)
    if (!is.na(match_idx)) {
      str1_chars <- str1_chars[-match_idx]
    }
  }
  
  return(paste(str1_chars, collapse = ""))
}

merge_datasets <- function(snp_caller, ground_truth, min_ccf, min_vaf_caller) {
  
  df <- snp_caller %>%
    dplyr::full_join(ground_truth, by = c("mutationID"), suffix = c("_caller", "_truth")) %>%
    dplyr::mutate(
      positive_truth = !is.na(CCF) & CCF > min_ccf,
      positive_call = !is.na(VAF_caller) & VAF_caller > min_vaf_caller,
      true_positive = positive_truth & positive_call,
      false_positive = !positive_truth & positive_call,
      false_negative = positive_truth & !positive_call,
      true_negative = !positive_truth & !positive_call
    ) %>% 
    dplyr::mutate(VAF_truth = ifelse(is.na(VAF_truth), 0, VAF_truth)) %>% 
    dplyr::mutate(VAF_caller = ifelse(is.na(VAF_caller), 0, VAF_caller)) %>% 
    dplyr::mutate(DP_truth = ifelse(is.na(DP_truth), 0, DP_truth)) %>% 
    dplyr::mutate(DP_caller = ifelse(is.na(DP_caller), 0, DP_caller)) %>% 
    dplyr::mutate(CCF = ifelse(is.na(CCF), 0, CCF))
  
  return(df)
}

get_colors = function(df) {
  all_levels <- levels(factor(df$FILTER))
  full_colors <- stats::setNames(RColorBrewer::brewer.pal(8, "Set2")[1:length(all_levels)], all_levels)  
  
  original_names = names(full_colors)
  if ("PASS" %in% names(full_colors)) {
    names(full_colors)[which(names(full_colors)=="PASS")] = original_names[1]
    names(full_colors)[1] = "PASS"
  }
  
  if ("Other" %in% names(full_colors)) {
    names(full_colors)[which(names(full_colors)=="Other")] = original_names[2]
    names(full_colors)[2] = "Other"
  }
  
  full_colors
}

get_sample_info <- function(vcf_path) {
  # Extract caller name from vcf_path
  caller_match <- base::regexpr("strelka|races|freebayes|mutect2|haplotypecaller", vcf_path)
  caller_name <- ifelse(caller_match != -1, base::regmatches(vcf_path, caller_match), NA)
  
  # Extract sample name from vcf_path
  sample_match <- base::regmatches(vcf_path, base::regexpr("SPN[0-9]+_[0-9]+\\.[0-9]+", vcf_path))
  sample_name <- ifelse(length(sample_match) > 0, sample_match, NA)
  
  list(sample_name = sample_name, caller_name = caller_name)
}

#' Generate Comprehensive Variant Calling Performance Report
#'
#' This function creates a detailed report comparing variant calling results
#' between a caller and ground truth data (ProCESS). It generates visualizations
#' and performance metrics for both all mutations and PASS-only mutations.
#'
#' @param seq_res_long Ground truth sequencing results (ProCESS data)
#' @param caller_res Variant caller results 
#' @param sample_info List containing sample metadata (spn, sample_id, caller_name, etc.)
#' @param min_ccf Minimum CCF threshold for considering mutations as true positives
#' @param min_vaf_caller Minimum VAF threshold for considering mutations called by caller
#'
#' @return List containing report_plot, report_metrics, and vaf_comparison
get_report <- function(seq_res_long, caller_res, sample_info, min_ccf, min_vaf_caller) {
  
  ###
  # DATA PREPROCESSING
  ###
  
  # NOTE: Filtering zero VAF mutations is commented out - may be needed for some callers
  # Remove 0 VAF mutations from caller
  #caller_res = caller_res %>% dplyr::filter(VAF != 0)
  
  # Group rare filter categories together for better visualization
  # Filters with frequency < 5% are grouped as "Other"
  caller_res <- caller_res %>% 
    dplyr::group_by(FILTER) %>% 
    dplyr::mutate(f = dplyr::n() / nrow(caller_res)) %>%  # Calculate filter frequency
    dplyr::mutate(FILTER = ifelse(FILTER == "PASS" | f > .05, FILTER, "Other"))
  
  ###
  # DATASET MERGING
  ###
  
  # Create merged dataset for all mutations (with VAF threshold)
  merged_df <- merge_datasets(snp_caller = caller_res, ground_truth = seq_res_long, min_ccf = min_ccf, min_vaf_caller = min_vaf_caller)
  
  # Filter to only PASS mutations for high-quality analysis
  caller_res_pass <- caller_res %>% 
    dplyr::filter(FILTER == "PASS" & !is.na(FILTER))
  
  # Create merged dataset with 0 threshold for PASS mutations (includes all VAFs)
  #merged_df0 <- merge_datasets(caller_res_pass, seq_res_long, min_ccf = 0, min_vaf_caller = 0)
  merged_df0 <- merge_datasets(caller_res_pass, seq_res_long, min_ccf = min_ccf, min_vaf_caller = min_vaf_caller)
  
  ###
  # PLOT GENERATION - FALSE NEGATIVE ANALYSIS
  ###
  
  # Analyze VAF distribution of false negative mutations
  p_false_negative_VAF_dist <- plot_false_negative_vaf_dist(merged_df0)
  
  ###
  # PLOT GENERATION - BASIC DISTRIBUTIONS
  ###
  
  # Prepare data lists for comparison plots
  DP_list = list(seq_res_long$DP, caller_res_pass$DP, caller_res$DP)
  VAF_list = list(seq_res_long$VAF, caller_res_pass$VAF, caller_res$VAF)
  labels = c("ProCESS", sample_info$caller_name, paste0(sample_info$caller_name, " (all)"))
  colors = method_colors[labels]
  
  # Generate distribution comparison plots
  DP_density = plot_density_comparison_multi(DP_list, labels, colors, x_label = "Read Depth")
  DP_ecdf = plot_ecdf_comparison_multi(DP_list, labels, colors, x_label = "Read Depth")
  
  VAF_density = plot_density_comparison_multi(VAF_list, labels, colors, x_label = "VAF")
  VAF_ecdf = plot_ecdf_comparison_multi(VAF_list, labels, colors, x_label = "VAF")
  
  ###
  # UNUSED CODE - INDIVIDUAL COVERAGE AND VAF PLOTS
  ###
  # The following plots were replaced by the multi-comparison plots above
  # Kept for potential future use if individual plots are needed
  
  # races_coverage <- plot_races_coverage(seq_res_long)
  # races_VAF <- plot_races_vaf(seq_res_long)
  # caller_coverage <- plot_caller_coverage(caller_res_pass, sample_info)
  # caller_VAF <- plot_caller_vaf(caller_res_pass, sample_info)
  
  ###
  # PLOT GENERATION - DIFFERENCE ANALYSIS
  ###
  
  # Analyze differences between caller and ground truth
  vaf_differences <- plot_vaf_difference(merged_df)      # VAF difference plots
  cov_differences <- plot_cov_difference(merged_df)      # Coverage difference plots
  
  ###
  # PLOT GENERATION - FILTER AND SCATTER ANALYSIS
  ###
  
  # Get consistent color palette for FILTER categories
  
  # Plot distribution of filter categories
  colors <- get_colors(merged_df)
  p_filter_dist <- plot_filter_distribution(merged_df, colors)
  
  # --- ALL MUTATIONS ANALYSIS ---
  # Generate scatter plots comparing caller vs ground truth for all mutations
  p_scatter_VAF_all <- plot_vaf_scatter_all(merged_df, colors, sample_info)
  p_scatter_DP_all <- plot_dp_scatter_all(merged_df, colors, sample_info)
  
  ###
  # PERFORMANCE METRICS - ALL MUTATIONS
  ###
  
  # Prepare binary classification vectors for performance calculation
  # Convert to 0/1 based on VAF threshold
  y_true <- as.numeric(factor((merged_df$CCF > min_ccf), levels=c(FALSE, TRUE))) - 1
  y_pred <- as.numeric(factor((merged_df$VAF_caller > min_vaf_caller), levels=c(FALSE, TRUE))) - 1
  
  # Generate Venn diagram showing overlap between caller and ground truth
  p_venn_all <- plot_venn_diagram(merged_df, sample_info$caller_name) +
    ggplot2::ggtitle("All called mutations")
  
  # Calculate comprehensive performance metrics (precision, recall, F1, etc.)
  metrics_all <- compute_metrics_from_vectors(y_true, y_pred)
  
  ###
  # PERFORMANCE METRICS - PASS MUTATIONS ONLY
  ###
  
  # Create merged dataset specifically for PASS mutations
  merged_df_pass <- merge_datasets(caller_res_pass, seq_res_long, min_ccf = min_ccf, min_vaf_caller = min_vaf_caller)
  
  # Generate scatter plots for high-quality (PASS) mutations only
  p_scatter_VAF_pass <- plot_vaf_scatter_pass(merged_df_pass, colors, sample_info)
  p_scatter_DP_pass <- plot_dp_scatter_pass(merged_df_pass, colors, sample_info)
  
  # Prepare binary classification vectors for PASS mutations
  y_true_pass <- as.numeric(factor((merged_df_pass$CCF > min_ccf), levels=c(FALSE, TRUE))) - 1
  y_pred_pass <- as.numeric(factor((merged_df_pass$VAF_caller > min_vaf_caller), levels=c(FALSE, TRUE))) - 1
  
  # Generate Venn diagram for PASS mutations only
  p_venn_pass <- plot_venn_diagram(merged_df_pass, sample_info$caller_name) +
    ggplot2::ggtitle("Only PASS mutations")
  
  # Calculate performance metrics for PASS mutations
  metrics_pass <- compute_metrics_from_vectors(y_true_pass, y_pred_pass)
  
  ###
  # METRICS COMBINATION AND ADVANCED ANALYSIS
  ###
  
  # Combine metrics from both all mutations and PASS-only analysis
  metrics <- dplyr::bind_rows(
    metrics_all %>% tidyr::pivot_longer(cols = colnames(metrics_all)) %>% 
      dplyr::mutate(Mutations = "All"),
    metrics_pass %>% tidyr::pivot_longer(cols = colnames(metrics_pass)) %>% 
      dplyr::mutate(Mutations = "Only Pass")  
  )
  
  # Perform detailed VAF performance analysis with tolerance
  # metrics_results = analyze_vaf_performance(seq_res_long, 
  #                                           caller_res, 
  #                                           only_pass = TRUE, 
  #                                           min_vaf_threshold = min_vaf, 
  #                                           vaf_tolerance_pct = 5)  # 5% VAF tolerance
  
  metrics_results = analyze_ccf_performance(seq_res_long, 
                                            caller_res,
                                            min_vaf_caller = min_vaf_caller,
                                            #pi = as.numeric(sample_info$purity), 
                                            only_pass = TRUE,
                                            min_ccf_threshold = min_ccf,
                                            tolerance_pct = 5)  # 5% VAF tolerance
  
  # Generate precision-recall curve
  #precision_recall_plot = plot_precision_recall_vaf(vaf_analysis_results = metrics_results)
  precision_recall_plot = plot_precision_recall_ccf(analysis_results = metrics_results)
  FDR = metrics_results$detection_summary["False Positive"] / (metrics_results$detection_summary["False Positive"] + metrics_results$detection_summary["True Positive"])
  FDR_rate_plot = dplyr::tibble(FDR = FDR) %>% 
    ggplot(mapping = aes(x = "FDR", y = FDR)) +
    geom_col() +
    lims(y = c(0,1)) +
    theme_bw() +
    labs(y = "Value", x = "")
  
  # ccf_bins = c(-Inf, 0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0)
  # ccf_bins = ccf_bins[ccf_bins < (pi * 0.95)]
  # ccf_bins = c(ccf_bins, (pi * 0.95), 1.0) %>% sort()
  # 
  # # Group filter categories
  # caller_res <- caller_res %>% 
  #   dplyr::group_by(FILTER) %>% 
  #   dplyr::mutate(f = dplyr::n() / nrow(caller_res)) %>% 
  #   dplyr::mutate(FILTER = ifelse(FILTER == "PASS" | f > .05, FILTER, "Other"))
  # 
  # if (only_pass) {
  #   caller_res = caller_res %>% dplyr::filter(FILTER == "PASS")
  # }
  # 
  # merged_df = merge_datasets(caller_res, seq_res_long, 0)
  # merged_df = merged_df %>% dplyr::select(mutationID, VAF_caller, VAF_truth, CCF) %>% 
  #   dplyr::mutate(CCF = ifelse(is.na(CCF), 0, CCF))
  # 
  # # Compute CCFs
  # merged_df <- merged_df %>%
  #   dplyr::mutate(
  #     CCF = CCF 
  #     #CCF_truth = VAF_truth / pi * (pi + 2 * (1 - pi)) / 2,
  #     #CCF_caller = VAF_caller / pi * (pi + 2 * (1 - pi)) / 2
  #   ) %>%
  #   dplyr::mutate(
  #     CCF = pmin(CCF, 1)
  #     #CCF_truth = pmin(CCF_truth, 1),
  #     #CCF_caller = pmin(CCF_caller, 1)
  #   )
  # 
  # # Create CCF bins
  # #ccf_labels <- c("0%", "0–5%", "5–10%", "10–20%", "20–30%", "30–50%", "50–100%")
  # 
  # ccf_labels = lapply(2:(length(ccf_bins)-2), function(i) {
  #   paste0(paste0(ccf_bins[i:(i+1)] * 100, collapse = "-"), "%")  
  # }) %>% unlist()
  # ccf_labels = c("0%", ccf_labels, "Clonal")
  # 
  # merged_df$CCF_bin <- cut(merged_df$CCF, 
  #                          breaks = ccf_bins, 
  #                          include.lowest = FALSE,
  #                          right = TRUE,
  #                          labels = ccf_labels)
  # 
  # merged_df %>% 
  #   ggplot(mapping = aes(x = VAF_truth, fill = CCF_bin)) +
  #   geom_histogram(alpha = .5, binwidth = 0.01, position = "identity") +
  #   lims(x = c(0,1))
  
  # Plot metrics comparison between all and PASS mutations
  p_metrics <- plot_metrics_comparison(metrics)
  
  # Show how metrics change across different VAF thresholds
  metrics_over_CCF = plot_metric_over_CCF_threshold(seq_res_long, caller_res, only_pass = TRUE, min_vaf_caller = min_vaf_caller) +
    ggplot2::geom_vline(xintercept = min_ccf, linetype = "dashed")  # Mark current threshold
  
  ###
  # CORRELATION ANALYSIS FOR TRUE POSITIVES
  ###
  
  # Analyze VAF correlation between caller and ground truth for correctly identified variants
  vaf_comparison <- merged_df_pass %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(true_positive) %>%  # Only analyze true positive variants
    dplyr::group_by(chr_caller) %>%   # Group by chromosome
    dplyr::summarise(
      n_variants = dplyr::n(),                             # Count of variants per chromosome
      cor_coeff = safe_cor_test(VAF_caller, VAF_truth),    # Correlation coefficient
      RMSE = safe_rmse(VAF_truth, VAF_caller),             # Root mean square error
      .groups = 'drop'
    )
  
  ###
  # REPORT LAYOUT DESIGN
  ###
  
  # Define grid layout for comprehensive report visualization
  # Each letter represents a plot position in the grid
  # design <- "
  # AABBCC
  # AABBCC
  # DDEEFF
  # DDEEFF
  # GGHHII
  # GGHHII
  # GGHHII
  # LLMMNN
  # LLMMNN
  # LLMMNN
  # OOOPPP
  # OOOPPP
  # "
  design <- "
  AABBCC
  AABBCC
  DDEEFF
  DDEEFF
  GGGHHH
  IIILMM
  IIILMM
  "
  
  # Generate descriptive title and subtitle for the report
  title <- paste0(sample_info$spn, 
                  ", calls by ", sample_info$caller_name)
  subtitle <- paste0("Only ", sample_info$mut_type, " mutations, purity = ", 
                     sample_info$purity, ", coverage = ", sample_info$coverage, "x")
  
  ###
  # UNUSED REPORT LAYOUT - ORIGINAL VERSION
  ###
  # This layout used individual coverage/VAF plots instead of comparison plots
  # Kept for reference in case the original layout is needed
  
  # report_plot <- patchwork::free(races_coverage) + patchwork::free(caller_coverage) + patchwork::free(p_filter_dist) +
  #   patchwork::free(races_VAF) + patchwork::free(caller_VAF) + patchwork::free(p_false_negative_VAF_dist) +
  #   patchwork::free(p_scatter_DP_all) + patchwork::free(p_scatter_VAF_all) + patchwork::free(p_venn_all) +
  #   patchwork::free(p_scatter_DP_pass) + patchwork::free(p_scatter_VAF_pass) + patchwork::free(p_venn_pass) +
  #   patchwork::free(p_metrics) + patchwork::free(metrics_over_VAF) +
  #   patchwork::plot_layout(design = design) +
  #   patchwork::plot_annotation(title, subtitle) & 
  #   ggplot2::theme(text = ggplot2::element_text(size = 12))
  
  ###
  # FINAL REPORT ASSEMBLY
  ###
  
  # Combine all plots into a comprehensive report using patchwork
  # free() function allows each plot to maintain its own scales
  report_plot <- patchwork::free(DP_density) + patchwork::free(DP_ecdf) + patchwork::free(p_scatter_DP_pass) +
    patchwork::free(VAF_density) + patchwork::free(VAF_ecdf) + patchwork::free(p_scatter_VAF_pass) +
    patchwork::free(p_false_negative_VAF_dist) +  patchwork::free(p_venn_pass) +
    patchwork::free(precision_recall_plot) + patchwork::free(FDR_rate_plot) + patchwork::free(metrics_over_CCF) +
    patchwork::plot_layout(design = design) +
    patchwork::plot_annotation(title, subtitle) & 
    ggplot2::theme(text = ggplot2::element_text(size = 12))
  
  ###
  # RETURN RESULTS
  ###
  
  # Return comprehensive results including:
  # - Complete visualization report
  # - Performance metrics table  
  # - VAF correlation analysis by chromosome
  list(
    report_plot = report_plot,                           # Combined visualization report
    report_metrics = metrics_results$performance_table,  # Detailed performance metrics
    vaf_comparison = vaf_comparison                      # Correlation analysis results
  )
}

# Function to perform VAF-stratified performance analysis
analyze_vaf_performance <- function(seq_res_long, caller_res, only_pass, 
                                    vaf_bins = c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0), 
                                    vaf_tolerance_pct = 5, 
                                    min_vaf_threshold = 0.02) {
  
  # Group filter categories for better visualization
  caller_res <- caller_res %>% 
    dplyr::group_by(FILTER) %>% 
    dplyr::mutate(f = dplyr::n() / nrow(caller_res)) %>% 
    dplyr::mutate(FILTER = ifelse(FILTER == "PASS" | f > .05, FILTER, "Other"))
  
  if (only_pass) {
    caller_res = caller_res %>% dplyr::filter(FILTER == "PASS")
  }
  
  merged_df = merge_datasets(caller_res, seq_res_long, 0)
  merged_df = merged_df %>% dplyr::select(mutationID, VAF_caller, VAF_truth, CCF) %>% 
    dplyr::mutate(CCF = ifelse(is.na(CCF), 0, CCF))
  
  # Create VAF bins based on truth VAF
  merged_df$VAF_bin <- cut(merged_df$VAF_truth, 
                           breaks = vaf_bins, 
                           include.lowest = TRUE,
                           right = FALSE,
                           labels = paste0(vaf_bins[-length(vaf_bins)]*100, "-", vaf_bins[-1]*100, "%"))
  
  # Calculate VAF difference (absolute and relative)
  merged_df$vaf_abs_diff <- abs(merged_df$VAF_truth - merged_df$VAF_caller)
  merged_df$vaf_rel_diff_pct <- abs(merged_df$VAF_truth - merged_df$VAF_caller) / merged_df$VAF_truth * 100
  
  # Classify variants based on detection and VAF accuracy
  merged_df$detection_status <- dplyr::case_when(
    # True Positives: Both must have VAF > threshold AND VAF difference within tolerance
    !is.na(merged_df$VAF_truth) & !is.na(merged_df$VAF_caller) & 
      merged_df$VAF_truth > min_vaf_threshold & merged_df$VAF_caller > min_vaf_threshold &
      merged_df$vaf_rel_diff_pct <= vaf_tolerance_pct ~ "True Positive",
    
    # False Negatives: Truth exists but caller missed it OR VAF difference too large
    !is.na(merged_df$VAF_truth) & merged_df$VAF_truth > min_vaf_threshold & 
      (is.na(merged_df$VAF_caller) | merged_df$VAF_caller <= min_vaf_threshold |
         (!is.na(merged_df$VAF_caller) & merged_df$vaf_rel_diff_pct > vaf_tolerance_pct)) ~ "False Negative",
    
    # False Positives: Caller detected but no truth OR truth VAF too low
    !is.na(merged_df$VAF_caller) & merged_df$VAF_caller > min_vaf_threshold &
      (is.na(merged_df$VAF_truth) | merged_df$VAF_truth <= min_vaf_threshold) ~ "False Positive",
    
    # TRUE NEGATIVES: Truth exists but VAF is below threshold AND caller correctly didn't call it (or called below threshold)
    !is.na(merged_df$VAF_truth) & merged_df$VAF_truth <= min_vaf_threshold &
      (is.na(merged_df$VAF_caller) | merged_df$VAF_caller <= min_vaf_threshold) ~ "True Negative",
    
    TRUE ~ "Unknown"
  )
  
  # Create VAF bins for True Negatives based on their truth VAF (even though below threshold)
  merged_df$VAF_bin_tn <- cut(merged_df$VAF_truth, 
                              breaks = vaf_bins, 
                              include.lowest = TRUE,
                              right = FALSE,
                              labels = paste0(vaf_bins[-length(vaf_bins)]*100, "-", vaf_bins[-1]*100, "%"))
  
  # Calculate performance metrics by VAF bin
  performance_by_bin <- merged_df %>%
    dplyr::filter(!is.na(VAF_bin) & VAF_truth > min_vaf_threshold) %>%  # Only consider variants above threshold for main analysis
    dplyr::group_by(VAF_bin) %>%
    dplyr::summarise(
      total_truth = dplyr::n(),
      true_positives = sum(detection_status == "True Positive"),
      false_negatives = sum(detection_status == "False Negative"),
      vaf_discordant = sum(!is.na(VAF_caller) & VAF_caller > min_vaf_threshold & 
                             vaf_rel_diff_pct > vaf_tolerance_pct, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::mutate(
      sensitivity = true_positives / (true_positives + false_negatives),
      sensitivity_pct = round(sensitivity * 100, 1)
    )
  
  # Calculate True Negatives and False Positives by VAF bin
  # For TN: group by the truth VAF bin (even though below threshold)
  tn_by_bin <- merged_df %>%
    dplyr::filter(detection_status == "True Negative") %>%
    dplyr::group_by(VAF_bin_tn) %>%
    dplyr::summarise(true_negatives = dplyr::n(), .groups = 'drop') %>%
    dplyr::rename(VAF_bin = VAF_bin_tn)
  
  # For FP: assign to VAF bins based on their called VAF
  fp_by_bin <- merged_df %>%
    dplyr::filter(detection_status == "False Positive") %>%
    dplyr::mutate(VAF_bin_fp = cut(VAF_caller, 
                                   breaks = vaf_bins, 
                                   include.lowest = TRUE,
                                   right = FALSE,
                                   labels = paste0(vaf_bins[-length(vaf_bins)]*100, "-", vaf_bins[-1]*100, "%"))) %>%
    dplyr::group_by(VAF_bin_fp) %>%
    dplyr::summarise(false_positives = dplyr::n(), .groups = 'drop') %>%
    dplyr::rename(VAF_bin = VAF_bin_fp)
  
  # Merge TN and FP counts with performance metrics
  performance_by_bin <- performance_by_bin %>%
    dplyr::left_join(tn_by_bin, by = "VAF_bin") %>%
    dplyr::left_join(fp_by_bin, by = "VAF_bin") %>%
    dplyr::mutate(
      true_negatives = ifelse(is.na(true_negatives), 0, true_negatives),
      false_positives = ifelse(is.na(false_positives), 0, false_positives),
      # Calculate FPR and Precision per bin
      fpr = false_positives / (false_positives + true_negatives),
      fpr_pct = round(fpr * 100, 3),
      precision = true_positives / (true_positives + false_positives),
      precision_pct = round(precision * 100, 1),
      # Calculate Specificity as well
      specificity = true_negatives / (true_negatives + false_positives),
      specificity_pct = round(specificity * 100, 1)
    )
  
  # Handle cases where denominator is 0 (no TN + FP in a bin)
  performance_by_bin <- performance_by_bin %>%
    dplyr::mutate(
      fpr = ifelse((false_positives + true_negatives) == 0, NA, fpr),
      fpr_pct = ifelse(is.na(fpr), NA, fpr_pct),
      specificity = ifelse((false_positives + true_negatives) == 0, NA, specificity),
      specificity_pct = ifelse(is.na(specificity), NA, specificity_pct),
      precision = ifelse((true_positives + false_positives) == 0, NA, precision),
      precision_pct = ifelse(is.na(precision), NA, precision_pct)
    )
  
  # Calculate overall metrics
  total_tp <- sum(merged_df$detection_status == "True Positive", na.rm = TRUE)
  total_fp <- sum(merged_df$detection_status == "False Positive", na.rm = TRUE)
  total_tn <- sum(merged_df$detection_status == "True Negative", na.rm = TRUE)
  total_fn <- sum(merged_df$detection_status == "False Negative", na.rm = TRUE)
  
  overall_precision <- ifelse((total_tp + total_fp) > 0, total_tp / (total_tp + total_fp), NA)
  overall_sensitivity <- ifelse((total_tp + total_fn) > 0, total_tp / (total_tp + total_fn), NA)
  overall_specificity <- ifelse((total_tn + total_fp) > 0, total_tn / (total_tn + total_fp), NA)
  overall_fpr <- ifelse((total_tn + total_fp) > 0, total_fp / (total_tn + total_fp), NA)
  
  # Calculate VAF accuracy for true positives
  vaf_accuracy_by_bin <- merged_df %>%
    dplyr::filter(detection_status == "True Positive") %>%
    dplyr::group_by(VAF_bin) %>%
    dplyr::summarise(
      n_variants = dplyr::n(),
      vaf_correlation = ifelse(dplyr::n() > 1, stats::cor(VAF_truth, VAF_caller, use = "complete.obs"), NA),
      mean_abs_error = mean(abs(VAF_truth - VAF_caller), na.rm = TRUE),
      median_abs_error = stats::median(abs(VAF_truth - VAF_caller), na.rm = TRUE),
      rmse = sqrt(mean((VAF_truth - VAF_caller)^2, na.rm = TRUE)),
      .groups = 'drop'
    )
  
  # Combine results
  results <- performance_by_bin %>%
    dplyr::left_join(vaf_accuracy_by_bin, by = "VAF_bin")
  
  return(list(
    performance_table = results,
    overall_metrics = list(
      precision = overall_precision,
      sensitivity = overall_sensitivity,
      specificity = overall_specificity,
      fpr = overall_fpr
    ),
    detection_summary = table(merged_df$detection_status),
    confusion_matrix = list(
      TP = total_tp,
      FP = total_fp,
      TN = total_tn,
      FN = total_fn
    ),
    raw_data = merged_df,
    vaf_tolerance_pct = vaf_tolerance_pct,
    min_vaf_threshold = min_vaf_threshold
  ))
}

analyze_ccf_performance <- function(seq_res_long, caller_res, only_pass, 
                                    #pi,
                                    #ccf_bins = c(-Inf, 0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0), 
                                    tolerance_pct = 5, 
                                    min_vaf_caller = 0.0,
                                    min_ccf_threshold = 0.02) {
  
  ccf_bins = c(-Inf, 0, 0.01, 0.02, 0.03, 0.05, 0.1, 0.25, 0.5, .99, 1.0)
  ccf_bins = ccf_bins[ccf_bins < 1.0]
  ccf_bins = c(ccf_bins, 1.0) %>% sort()
  
  # Group filter categories
  caller_res <- caller_res %>% 
    dplyr::group_by(FILTER) %>% 
    dplyr::mutate(f = dplyr::n() / nrow(caller_res)) %>% 
    dplyr::mutate(FILTER = ifelse(FILTER == "PASS" | f > .05, FILTER, "Other"))
  
  if (only_pass) {
    caller_res = caller_res %>% dplyr::filter(FILTER == "PASS")
  }
  
  merged_df = merge_datasets(caller_res, seq_res_long, min_ccf = min_ccf_threshold, min_vaf_caller = min_vaf_caller)
  merged_df = merged_df %>% dplyr::select(mutationID, VAF_caller, VAF_truth, CCF) %>% 
    dplyr::mutate(CCF = ifelse(is.na(CCF), 0, CCF))
  
  ccf_labels = lapply(2:(length(ccf_bins)-2), function(i) {
    paste0(paste0(ccf_bins[i:(i+1)] * 100, collapse = "-"), "%")  
  }) %>% unlist()
  ccf_labels = c("0%", ccf_labels, "Clonal")
  
  merged_df$CCF_bin <- cut(merged_df$CCF, 
                           breaks = ccf_bins, 
                           include.lowest = FALSE,
                           right = TRUE,
                           labels = ccf_labels)
  
  # Calculate differences
  merged_df$VAF_abs_diff <- abs(merged_df$VAF_truth - merged_df$VAF_caller)
  merged_df$VAF_rel_diff_pct <- abs(merged_df$VAF_truth - merged_df$VAF_caller) / merged_df$VAF_truth * 100
  
  # Classify detection status
  merged_df$detection_status <- dplyr::case_when(
    # True Positives: Both must have VAF > threshold AND VAF difference within tolerance
    !is.na(merged_df$VAF_truth) & !is.na(merged_df$VAF_caller) & 
      merged_df$CCF > min_ccf_threshold &
      merged_df$VAF_rel_diff_pct <= tolerance_pct ~ "True Positive",
    
    # False Negatives: Truth exists but caller missed it OR VAF difference too large
    !is.na(merged_df$VAF_truth) & merged_df$CCF > min_ccf_threshold & 
      (is.na(merged_df$VAF_caller) |
         (!is.na(merged_df$VAF_caller) & merged_df$VAF_rel_diff_pct > tolerance_pct)) ~ "False Negative",
    
    # False Positives: Caller detected but no truth OR truth VAF too low
    !is.na(merged_df$VAF_caller) & merged_df$VAF_caller > min_ccf_threshold &
      (is.na(merged_df$VAF_caller) | merged_df$CCF <= min_ccf_threshold) ~ "False Positive",
    
    # TRUE NEGATIVES: Truth exists but VAF is below threshold AND caller correctly didn't call it (or called below threshold)
    !is.na(merged_df$VAF_caller) &
      (is.na(merged_df$VAF_caller) | merged_df$CCF <= min_ccf_threshold) ~ "True Negative",
    
    TRUE ~ "Unknown"
  )
  
  # Bin true negatives
  merged_df$CCF_bin_tn <- cut(merged_df$CCF, 
                              breaks = ccf_bins, 
                              include.lowest = FALSE,
                              right = TRUE,
                              labels = ccf_labels)
  
  # Performance by CCF bin
  performance_by_bin <- merged_df %>%
    dplyr::filter(!is.na(CCF_bin) & CCF > min_ccf_threshold) %>%
    dplyr::group_by(CCF_bin) %>%
    dplyr::summarise(
      total_truth = dplyr::n(),
      true_positives = sum(detection_status == "True Positive"),
      false_negatives = sum(detection_status == "False Negative"),
      false_positive = sum(detection_status == "False Positive"),
      VAF_discordant = sum(!is.na(VAF_caller) & CCF > min_ccf_threshold &
                             VAF_rel_diff_pct > tolerance_pct, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::mutate(
      sensitivity = true_positives / (true_positives + false_negatives),
      sensitivity_pct = round(sensitivity * 100, 1)
    )
  
  tn_by_bin <- merged_df %>%
    dplyr::filter(detection_status == "True Negative") %>%
    dplyr::group_by(CCF_bin_tn) %>%
    dplyr::summarise(true_negatives = dplyr::n(), .groups = 'drop') %>%
    dplyr::rename(CCF_bin = CCF_bin_tn)
  
  fp_by_bin <- merged_df %>%
    dplyr::filter(detection_status == "False Positive") %>%
    dplyr::mutate(CCF_bin_fp = cut(CCF, 
                                   breaks = ccf_bins, 
                                   include.lowest = FALSE,
                                   right = TRUE,
                                   labels = ccf_labels)) %>%
    dplyr::group_by(CCF_bin_fp) %>%
    dplyr::summarise(false_positives = dplyr::n(), .groups = 'drop') %>%
    dplyr::rename(CCF_bin = CCF_bin_fp)
  
  performance_by_bin <- performance_by_bin %>%
    dplyr::left_join(tn_by_bin, by = "CCF_bin") %>%
    dplyr::left_join(fp_by_bin, by = "CCF_bin") %>%
    dplyr::mutate(
      true_negatives = ifelse(is.na(true_negatives), 0, true_negatives),
      false_positives = ifelse(is.na(false_positives), 0, false_positives),
      fpr = false_positives / (false_positives + true_negatives),
      fpr_pct = round(fpr * 100, 3),
      precision = true_positives / (true_positives + false_positives),
      precision_pct = round(precision * 100, 1),
      specificity = true_negatives / (true_negatives + false_positives),
      specificity_pct = round(specificity * 100, 1)
    ) %>%
    dplyr::mutate(
      fpr = ifelse((false_positives + true_negatives) == 0, NA, fpr),
      fpr_pct = ifelse(is.na(fpr), NA, fpr_pct),
      specificity = ifelse((false_positives + true_negatives) == 0, NA, specificity),
      specificity_pct = ifelse(is.na(specificity), NA, specificity_pct),
      precision = ifelse((true_positives + false_positives) == 0, NA, precision),
      precision_pct = ifelse(is.na(precision), NA, precision_pct)
    )
  
  total_tp <- sum(merged_df$detection_status == "True Positive", na.rm = TRUE)
  total_fp <- sum(merged_df$detection_status == "False Positive", na.rm = TRUE)
  total_tn <- sum(merged_df$detection_status == "True Negative", na.rm = TRUE)
  total_fn <- sum(merged_df$detection_status == "False Negative", na.rm = TRUE)
  
  overall_precision <- ifelse((total_tp + total_fp) > 0, total_tp / (total_tp + total_fp), NA)
  overall_sensitivity <- ifelse((total_tp + total_fn) > 0, total_tp / (total_tp + total_fn), NA)
  overall_specificity <- ifelse((total_tn + total_fp) > 0, total_tn / (total_tn + total_fp), NA)
  overall_fpr <- ifelse((total_tn + total_fp) > 0, total_fp / (total_tn + total_fp), NA)
  
  ccf_accuracy_by_bin <- merged_df %>%
    dplyr::filter(detection_status == "True Positive") %>%
    dplyr::group_by(CCF_bin) %>%
    dplyr::summarise(
      n_variants = dplyr::n(),
      vaf_correlation = ifelse(dplyr::n() > 1, stats::cor(VAF_truth, VAF_caller, use = "complete.obs"), NA),
      mean_abs_error = mean(VAF_abs_diff, na.rm = TRUE),
      median_abs_error = stats::median(VAF_abs_diff, na.rm = TRUE),
      rmse = sqrt(mean(VAF_abs_diff^2, na.rm = TRUE)),
      .groups = 'drop'
    )
  
  results <- performance_by_bin %>%
    dplyr::left_join(ccf_accuracy_by_bin, by = "CCF_bin")
  
  return(list(
    performance_table = results,
    overall_metrics = list(
      precision = overall_precision,
      sensitivity = overall_sensitivity,
      specificity = overall_specificity,
      fpr = overall_fpr
    ),
    detection_summary = table(merged_df$detection_status),
    confusion_matrix = list(
      TP = total_tp,
      FP = total_fp,
      TN = total_tn,
      FN = total_fn
    ),
    raw_data = merged_df,
    tolerance_pct = tolerance_pct,
    min_ccf_threshold = min_ccf_threshold
  ))
}

get_multi_caller_report <- function(seq_res_long, caller_res_list, sample_info, min_ccf, min_vaf_caller, only_pass) {
  if (only_pass) {
    names_to_keep = names(caller_res_list)
    caller_res_list = lapply(caller_res_list, function(cr) {
      cr %>% dplyr::filter(FILTER == "PASS")
    })
    names(caller_res_list) = names_to_keep
  }
  
  ###
  # PLOT GENERATION - BASIC DISTRIBUTIONS
  ###
  
  # Prepare data lists for comparison plots
  DP_list = list(seq_res_long$DP, caller_res_list$mutect2$DP, caller_res_list$strelka$DP, caller_res_list$freebayes$DP)
  VAF_list = list(seq_res_long$VAF, caller_res_list$mutect2$VAF, caller_res_list$strelka$VAF, caller_res_list$freebayes$VAF)
  labels = c("ProCESS", "mutect2", "strelka", "freebayes")
  colors <- method_colors[labels]
  
  # Generate distribution comparison plots
  DP_density = plot_density_comparison_multi(DP_list, labels, colors, x_label = "Read Depth")
  DP_ecdf = plot_ecdf_comparison_multi(DP_list, labels, colors, x_label = "Read Depth")
  
  VAF_density = plot_density_comparison_multi(VAF_list, labels, colors, x_label = "VAF")
  VAF_ecdf = plot_ecdf_comparison_multi(VAF_list, labels, colors, x_label = "VAF")
  
  # Generate Venn diagram showing overlap between caller and ground truth
  x = list(
    "ProCESS" = seq_res_long %>% dplyr::filter(CCF > min_ccf) %>% dplyr::pull(mutationID),
    "mutect2" = caller_res_list$mutect2 %>% dplyr::filter(VAF > min_vaf_caller) %>% dplyr::pull(mutationID),
    "strelka" = caller_res_list$strelka %>% dplyr::filter(VAF > min_vaf_caller) %>% dplyr::pull(mutationID),
    "freebayes" = caller_res_list$freebayes %>% dplyr::filter(VAF > min_vaf_caller) %>% dplyr::pull(mutationID)
  )
  
  # Get all unique mutation IDs from all methods
  all_mutations <- unique(c(
    seq_res_long %>% dplyr::filter(CCF > min_ccf) %>% dplyr::pull(mutationID),
    caller_res_list$mutect2 %>% dplyr::filter(VAF > min_vaf_caller) %>% dplyr::pull(mutationID),
    caller_res_list$strelka %>% dplyr::filter(VAF > min_vaf_caller) %>% dplyr::pull(mutationID),
    caller_res_list$freebayes %>% dplyr::filter(VAF > min_vaf_caller) %>% dplyr::pull(mutationID)
  ))
  
  # Create the tibble
  mutation_tibble <- tibble(
    mutationID = all_mutations,
    ProCESS = as.integer(mutationID %in% (seq_res_long %>% dplyr::filter(CCF > min_ccf) %>% dplyr::pull(mutationID))),
    mutect2 = as.integer(mutationID %in% (caller_res_list$mutect2 %>% dplyr::filter(VAF > min_vaf_caller) %>% dplyr::pull(mutationID))),
    strelka = as.integer(mutationID %in% (caller_res_list$strelka %>% dplyr::filter(VAF > min_vaf_caller) %>% dplyr::pull(mutationID))),
    freebayes = as.integer(mutationID %in% (caller_res_list$freebayes %>% dplyr::filter(VAF > min_vaf_caller) %>% dplyr::pull(mutationID)))
  ) %>%
    # Add VAF from seq_res_long
    dplyr::left_join(
      seq_res_long %>% dplyr::select(mutationID, VAF, CCF),
      by = "mutationID"
    ) %>%
    # Replace NA VAF values with 0
    dplyr::mutate(VAF = ifelse(is.na(VAF), 0, VAF), CCF = ifelse(is.na(CCF), 0, CCF))

  #print(head(mutation_tibble))
  
  # Upset plot 
  # CCF binning
  ccf_bins <- c(-Inf, 0, .01, .02, .03, 0.05, 0.1, 0.2, 0.3, 0.5, 0.99, 1.0)
  
  ccf_labels <- sapply(2:(length(ccf_bins) - 2), function(i) {
    paste0(paste0(ccf_bins[i:(i + 1)] * 100, collapse = "–"), "%")
  })
  ccf_labels <- c("0%", ccf_labels, "Clonal")
  
  # Apply bins to your plotting tibble
  mutation_tibble$CCF_bin <- cut(
    mutation_tibble$CCF,
    breaks = ccf_bins,
    include.lowest = FALSE,
    right = TRUE,
    labels = ccf_labels
  )
  
  # Palette that adapts to the number of bins
  base_cols <- c('#f0f9e8','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#08589e')
  pal <- if (length(ccf_labels) <= length(base_cols)) {
    base_cols[seq_len(length(ccf_labels))]
  } else {
    grDevices::colorRampPalette(base_cols)(length(ccf_labels))
  }
  names(pal) <- ccf_labels
  
  # Upset plot
  methods <- names(x)  # as in your code
  upset_plot <- ComplexUpset::upset(
    mutation_tibble,
    methods,
    base_annotations = list(
      "Intersection size" = ComplexUpset::intersection_size(
        counts = FALSE,
        mapping = aes(fill = CCF_bin)
      ) +
        scale_fill_manual(values = pal, drop = FALSE)
    ),
    width_ratio = 0.1
  )
  # ccf_bins = c(-Inf, 0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0)
  # mutation_tibble$CCF = mutation_tibble$VAF * pi^(-1) * (pi * 1 + (1 - pi) * 2) / 2
  # mutation_tibble$CCF[mutation_tibble$CCF > 1] = 1
  # 
  # # Create custom labels
  # ccf_labels = c("0%", "0–5%", "5–10%", "10–20%", "20–30%", "30–50%", "50–100%")
  # 
  # mutation_tibble$CCF_bin <- cut(mutation_tibble$CCF, 
  #                                breaks = ccf_bins, 
  #                                include.lowest = FALSE,
  #                                right = TRUE,
  #                                labels = ccf_labels)
  # methods = names(x)
  # upset_plot <- ComplexUpset::upset(mutation_tibble, methods, base_annotations=list(
  #   'Intersection size'=ComplexUpset::intersection_size(
  #     counts=FALSE,
  #     mapping=aes(fill=CCF_bin)
  #   ) +
  #     scale_fill_manual(values=c('#f0f9e8','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#08589e'))
  # ),
  # width_ratio=0.1) 
  
  # upset_plot = ggVennDiagram(x, force_upset = TRUE)
  # upset_plot = ggplotify::as.ggplot(upset_plot)
  
  ###
  # METRICS COMBINATION AND ADVANCED ANALYSIS
  ###
  
  # Perform detailed VAF performance analysis with tolerance
  # vaf_analysis_results_across_callers = lapply(names(caller_res_list), function(nc) {
  #   caller_res = caller_res_list[[nc]]
  #   metrics_results = analyze_vaf_performance(seq_res_long, 
  #                                             caller_res, 
  #                                             only_pass = TRUE, 
  #                                             min_vaf_threshold = min_vaf, 
  #                                             vaf_tolerance_pct = 5)  # 5% VAF tolerance
  #   metrics_results$performance_table %>% dplyr::mutate(caller = nc)
  # }) %>% do.call("bind_rows", .)
  
  ccf_analysis_results_across_callers = lapply(names(caller_res_list), function(nc) {
    caller_res = caller_res_list[[nc]]
    
    metrics_results = analyze_ccf_performance(
      seq_res_long, caller_res, only_pass = TRUE, tolerance_pct = 5, 
      min_vaf_caller = min_vaf_caller, min_ccf_threshold = min_ccf
    )
    
    metrics_results$performance_table %>% dplyr::mutate(caller = nc)
  }) %>% do.call("bind_rows", .)
  
  # Generate precision-recall curve
  precision_recall_plot = plot_multicaller_precision_recall(performance_data = ccf_analysis_results_across_callers, x_name = "CCF_bin")
  FDR_plot = plot_multicaller_FDR(seq_res_long, caller_res_list, min_vaf_caller, min_ccf)  

  #rmse_plot = plot_multicaller_rmse_vaf(vaf_analysis_results_across_callers)
  
  ###
  # REPORT LAYOUT DESIGN
  ###
  
  # Define grid layout for comprehensive report visualization
  # Each letter represents a plot position in the grid
  
  design <- "
  AAAAABBBBB
  AAAAABBBBB
  CCCCCDDDDD
  CCCCCDDDDD
  EEEEFGGGGG
  EEEEFGGGGG
  EEEEFGGGGG
  "
  
  # Generate descriptive title and subtitle for the report
  title <- paste0(sample_info$spn, ", sample ", sample_info$sample_id)
  subtitle <- paste0(sample_info$mut_type, " mutations, purity = ", 
                     sample_info$purity, ", coverage = ", sample_info$coverage, "x")
  
  # Combine all plots into a comprehensive report using patchwork
  # free() function allows each plot to maintain its own scales
  report_plot <- free(DP_density) + free(DP_ecdf) +
    free(VAF_density) + free(VAF_ecdf) + 
    free(precision_recall_plot) + free(FDR_plot) + free(upset_plot) +
    plot_layout(design = design) +
    plot_annotation(title, subtitle) & 
    theme(text = element_text(size = 12))
  
  report_plot
}
