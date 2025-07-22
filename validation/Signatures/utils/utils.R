### Align ProCESS and tumourevo signatures/exposures data ###

align_sparsesig_res <- function(sparse_list) {
  aligned <- list()
  
  for (spn in names(sparse_list)) {
    aligned[[spn]] <- list()
    
    for (coverage in names(sparse_list[[spn]])) {
      aligned[[spn]][[coverage]] <- list()
      
      for (purity in names(sparse_list[[spn]][[coverage]])) {
        mat <- sparse_list[[spn]][[coverage]][[purity]]
        
        if (is.null(mat) || !is.matrix(mat)) next
        
        # Coerce to numeric matrix
        mat <- apply(mat, 2, as.numeric)
        
        # Assign proper rownames consistently
        n_samples <- nrow(mat)
        sample_ids <- paste0(spn, "_1.", seq_len(n_samples))
        rownames(mat) <- sample_ids
        
        # Normalize rows to proportions
        mat <- t(apply(mat, 1, function(x) if (sum(x) == 0) x else x / sum(x)))
        
        # Store
        aligned[[spn]][[coverage]][[purity]] <- as.matrix(mat)
      }
    }
  }
  
  return(aligned)
}


align_sigprofiler_res <- function(sigprofiler_list) {
  aligned <- list()
  
  for (spn in names(sigprofiler_list)) {
    aligned[[spn]] <- list()
    
    for (coverage in names(sigprofiler_list[[spn]])) {
      aligned[[spn]][[coverage]] <- list()
      
      for (purity in names(sigprofiler_list[[spn]][[coverage]])) {
        df <- sigprofiler_list[[spn]][[coverage]][[purity]][["SigProfiler_COSMIC_exposure"]]
        
        if (is.null(df)) next
        
        df <- as.data.frame(df)
        if (!"Samples" %in% colnames(df)) next
        
        df <- df %>%
          tibble::column_to_rownames("Samples") %>%
          mutate(across(everything(), as.numeric))
        
        df_norm <- t(apply(df, 1, function(x) if (sum(x) == 0) x else x / sum(x)))
        if (is.vector(df_norm)) {
          df_norm <- matrix(df_norm, nrow = 1, dimnames = list(rownames(df), colnames(df)))
        }
        
        df <- as.data.frame(df_norm)
        
        # Clean sample names: extract "SPN01_1.1" from long names
        pattern <- "^.*?_(SPN\\d+_\\d+\\.\\d+)$"
        rn <- rownames(df)
        rn_new <- ifelse(grepl(pattern, rn), sub(pattern, "\\1", rn), rn)
        rownames(df) <- rn_new
        
        aligned[[spn]][[coverage]][[purity]] <- df
      }
    }
  }
  
  return(aligned)
}


### Reshape data ###

reshape_exposures_long <- function(exposures_mat, spn, coverage = NA, purity = NA, method_name) {
  df <- as.data.frame(exposures_mat)
  df$Sample_ID <- rownames(exposures_mat)
  
  df_long <- df %>%
    pivot_longer(-Sample_ID, names_to = "Signature", values_to = "Exposure") %>%
    mutate(SPN = spn,
           Coverage = coverage,
           Purity = purity,
           Method = method_name)
  
  return(df_long)
}

extract_ground_truth_long <- function(ground_truth_list) {
  out <- list()

  for (spn in names(ground_truth_list)) {
    for (coverage in names(ground_truth_list[[spn]])) {
      for (purity in names(ground_truth_list[[spn]][[coverage]])) {
        exposures_mat <- ground_truth_list[[spn]][[coverage]][[purity]]

        if (is.null(exposures_mat)) next

        long_df <- reshape_exposures_long(exposures_mat, spn, coverage, purity, "ProCESS")
        out[[paste(spn, coverage, purity, sep = "_")]] <- long_df
      }
    }
  }

  do.call(rbind, out)
}


extract_sparsesig_long <- function(sparsesig_aligned) {
  out <- list()

  for (spn in names(sparsesig_aligned)) {
    for (coverage in names(sparsesig_aligned[[spn]])) {
      for (purity in names(sparsesig_aligned[[spn]][[coverage]])) {
        df <- sparsesig_aligned[[spn]][[coverage]][[purity]]

        if (is.null(df)) next

        # df is samples x signatures matrix/data.frame, samples are rownames
        long_df <- reshape_exposures_long(df, spn, coverage, purity, "SparseSignatures")
        out[[paste(spn, coverage, purity, sep = "_")]] <- long_df
      }
    }
  }

  do.call(rbind, out)
}


extract_sigprofiler_long <- function(sigprof_aligned) {
  out <- list()

  for (spn in names(sigprof_aligned)) {
    for (coverage in names(sigprof_aligned[[spn]])) {
      for (purity in names(sigprof_aligned[[spn]][[coverage]])) {
        df <- sigprof_aligned[[spn]][[coverage]][[purity]]

        if (is.null(df)) next

        long_df <- reshape_exposures_long(df, spn, coverage, purity, "SigProfiler")
        out[[paste(spn, coverage, purity, sep = "_")]] <- long_df
      }
    }
  }

  do.call(rbind, out)
}


prepare_sankey_data_sbs <- function(ground_truth_list, sparsesig_aligned, sigprof_aligned) {
  gt_long <- extract_ground_truth_long(ground_truth_list)
  sparsesig_long <- extract_sparsesig_long(sparsesig_aligned)
  sigprof_long <- extract_sigprofiler_long(sigprof_aligned)

  combined <- bind_rows(gt_long, sparsesig_long, sigprof_long)

  # Filter to keep only exposures > 0
  combined <- combined %>% filter(Exposure > 0)

  return(combined)
}


prepare_sankey_data_id <- function(ground_truth_list, sigprof_aligned) {
  gt_long <- extract_ground_truth_long(ground_truth_list)
  sigprof_long <- extract_sigprofiler_long(sigprof_aligned)

  combined <- bind_rows(gt_long, sigprof_long)

  # Filter to keep only exposures > 0
  combined <- combined %>% filter(Exposure > 0)

  return(combined)
}



