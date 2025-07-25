pkgs <- c("ProCESS", "tidyverse", "ggplot2", 
          "caret", "ggtext", "reshape2", 
          "lsa", "Metrics", "ggalluvial", 
          "patchwork")
sapply(pkgs, require, character.only = TRUE)

source("../../getters/process_getters.R")
source("../../getters/tumourevo_getters.R")
source("utils/utils_getters.R")
source("utils/utils_validation.R")
source("utils/utils.R")
source("utils/utils_plots.R")
source("utils/utils_sparsesig.R")

### Get ProCESS exposure data ###
base_path <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
coverage_list <- c(50)
purity_list <- c(0.9)
dataset <- "SCOUT"
spn <- c("SPN01", "SPN02", "SPN03")
context <- 'SBS'
if (context == 'SBS'){
  context_sig = 'SBS96'
}
vcf_caller <- "mutect2"
cna_caller <- "ascat"

# Load GT data
ground_truth_nested <- list()
for (spn_id in spn) {
  print(spn_id)
  
  if (is.null(ground_truth_nested[[spn_id]])) {
    ground_truth_nested[[spn_id]] <- list()
  }
  
  for (cov in coverage_list) {
    print(cov)
    coverage_key <- paste0("coverage_", cov)
    
    if (is.null(ground_truth_nested[[spn_id]][[coverage_key]])) {
      ground_truth_nested[[spn_id]][[coverage_key]] <- list()
    }
    
    for (pur in purity_list) {
      print(pur)
      
      purity_key <- paste0("purity_", pur)
      gt_exposure <-  get_process_exposures(spn = spn_id, 
                                            coverage = cov,
                                            purity = pur)
      ground_truth_nested[[spn_id]][[coverage_key]][[purity_key]] <- gt_exposure[[context]]  %>%
        tibble::column_to_rownames("Sample_ID") %>% 
        as.matrix()
    }
  }
}


### Load tumourevo signature data ###
tumourevo_signature_res <- list()
for (spn_id in spn) {
  tumourevo_signature_res[[spn_id]] <- list()
  
  for (cov in coverage_list) {
    tumourevo_signature_res[[spn_id]][[paste0("coverage_", cov)]] <- list()
    
    for (pur in purity_list) {
      message("Processing ", spn_id, " | Coverage: ", cov, " | Purity: ", pur)
      
      result <- tryCatch({
        # Get SparseSignatures paths
        sparsesig <- get_tumourevo_signatures(
          spn = spn_id,
          coverage = cov,
          purity = pur,
          vcf_caller = vcf_caller,
          cna_caller = cna_caller,
          tool = "SparseSignatures"
        )
        
        # Get SigProfiler paths
        sigprofiler <- get_tumourevo_signatures(
          spn = spn_id,
          coverage = cov,
          purity = pur,
          vcf_caller = vcf_caller,
          cna_caller = cna_caller,
          tool = "SigProfiler",
          context = context_sig
        )
        
        # Combine and load paths
        paths <- c(
          sparsesig$nmf_Lasso_out,
          sparsesig$cv_means_mse,
          sparsesig$best_params_config,
          sparsesig$mut_counts,
          sigprofiler$COSMIC_exposure,
          sigprofiler$COSMIC_signatures,
          sigprofiler$denovo_exposure,
          sigprofiler$denovo_signatures
        )
        
        data <- load_signature_data(paths)
        
        names(data) <- c(
          "SparseSig_nmf_Lasso_out",
          "SparseSig_cv_means_mse",
          "SparseSig_best_params_config",
          "SparseSig_mut_counts",
          "SigProfiler_COSMIC_exposure",
          "SigProfiler_COSMIC_signatures",
          "SigProfiler_denovo_exposure",
          "SigProfiler_denovo_signatures"
        )
        
        data
      }, error = function(e) {
        message("Failed for ", spn_id, " | cov: ", cov, " | pur: ", pur)
        message("Reason: ", e$message)
        NULL
      })
      
      # Store result if successful
      if (!is.null(result)) {
        tumourevo_signature_res[[spn_id]][[paste0("coverage_", cov)]][[paste0("purity_", pur)]] <- result
      }
    }
  }
}


### Map De novo SparseSignatures results to COSMIC using cosine similarity ###
# Iterate across combinations of coverage and purity for a given SPN
cosmic_path <- 'COSMIC_ref/COSMIC_v3.4_SBS_GRCh38.txt'
mapped_res <- list()
for (spn_id in spn) {
  mapped_res[[spn_id]] <- list()
  
  for (cov in coverage_list) {
    c <- paste0("coverage_", cov)
    mapped_res[[spn_id]][[c]] <- list()
    
    for (pur in purity_list) {
      p <- paste0("purity_", pur)
      sig_out <- tumourevo_signature_res[[spn_id]][[c]][[p]][["SparseSig_nmf_Lasso_out"]]
      
      if (is.null(sig_out)) {
        warning(paste("Missing sig_out for", spn_id, cov, pur))
        mapped_res[[spn_id]][[c]][[p]] <- NULL
        next
      }
      
      remapped <- tryCatch({
        map_sparsesig_to_cosmic(
          sparsesig_out = sig_out,
          cosmic_path = cosmic_path,
          threshold = 0.5
        )
      }, error = function(e) {
        message(paste("Error in", spn_id, cov, pur, ":", e$message))
        return(NULL)
      })
      
      mapped_res[[spn_id]][[c]][[p]] <- remapped
    }
  }
}

### Summarize Signatures across combinations ###
sparsig_cosmic <- mapped_res
sparsesig_aligned <- align_sparsesig_res(sparsig_cosmic)
sigprof_aligned <- align_sigprofiler_res(tumourevo_signature_res)

### Compare exposure of estimated and true signatures  ###
sankey_df <- prepare_sankey_data_sbs(ground_truth_nested, sparsesig_aligned, sigprof_aligned)

<<<<<<< HEAD
plots <- list()
for (spn_id in spn){
  for (cov in coverage_list){
    for (pur in purity_list){
      sankey_plots <- generate_sankey_sbs(sankey_df = sankey_df, spn_id = spn_id, cov = cov, pur = pur)
      plots[[paste0('coverage_',cov, "_purity_", pur)]][[spn_id]] <- sankey_plots
=======
sankey_df <- prepare_sankey_data(ground_truth_nested, sparsesig_aligned, sigprof_aligned)

# Substitution: SPN##_X.1 -> SPN##_1.X
sankey_df$Sample_ID <- gsub("^(SPN\\d+)_([0-9]+)\\.1$", "\\1_1.\\2", sankey_df$Sample_ID)

saveRDS(sankey_df, file = "validation/sankey_signatures_sbs.rds")

sankey_df <- sankey_df %>%
  dplyr::mutate(
    Coverage = as.numeric(gsub("coverage_", "", Coverage)),
    Purity = as.numeric(gsub("purity_", "", Purity))
  ) %>% 
  dplyr::filter(Coverage == 50, Purity == 0.9)


spns <- c("SPN01")
sankey_plots <- lapply(spns, function(spn) generate_sankey_sbs(spn, cov = 50, pur = 0.9))

wrapped_sankey <- wrap_plots(sankey_plots, ncol = 1)
wrapped_sankey

saveRDS(wrapped_sankey, file = "validation/sankey_SPN01_50x_0.9p_sbs.rds")
ggsave("validation/sankey_SPN01_50x_0.9p_sbs.pdf", plot = wrapped_sankey, width = 12, height = 5, units = "in")



### Exposure validation ###

# Align exposure data

aligned_exposures <- list(
  SparseSignatures = list(),
  SigProfiler = list()
)

spn_list <- names(ground_truth_nested)

for (spn in spn_list) {
  coverage_list <- names(ground_truth_nested[[spn]])

  for (coverage in coverage_list) {
    purity_list <- names(ground_truth_nested[[spn]][[coverage]])

    for (purity in purity_list) {
      # Access exposure matrices
      gt <- ground_truth_nested[[spn]][[coverage]][[purity]]
      sparse <- sparsesig_aligned[[spn]][[coverage]][[purity]]
      sigprof <- sigprof_aligned[[spn]][[coverage]][[purity]]

      # Find shared signatures
      shared_sparse <- intersect(colnames(gt), colnames(sparse))
      shared_sigprof <- intersect(colnames(gt), colnames(sigprof))

      # Subset to shared signatures
      gt_sparse <- gt[, shared_sparse, drop = FALSE]
      sparse <- sparse[, shared_sparse, drop = FALSE]

      gt_sigprof <- gt[, shared_sigprof, drop = FALSE]
      sigprof <- sigprof[, shared_sigprof, drop = FALSE]

      # Create unique key for storage
      key <- paste(spn, coverage, purity, sep = "_")

      # Store aligned pairs
      aligned_exposures$SparseSignatures[[key]] <- list(gt = gt_sparse, tool = sparse)
      aligned_exposures$SigProfiler[[key]]     <- list(gt = gt_sigprof, tool = sigprof)
>>>>>>> 4d0c19519efbff18dbcd6c28e31b0b55d23837e2
    }
  }
}

plots_all <- list()
for (cov in coverage_list){
  for (pur in purity_list){
    plots_all[[paste0('coverage_',cov, "_purity_", pur)]] <- wrap_plots(plots[[paste0('coverage_',cov, "_purity_", pur)]], nrow = length(spn)) 
  }
}

# Signature and Exposure Validation
metrics_sample <- tibble()
metrics_spn <- tibble()
cosine_mse <- tibble()
for (spn_id in spn) {
  
  for (cov in coverage_list) {
    c <- paste0("coverage_", cov)
    
    for (pur in purity_list) {
      p <- paste0("purity_", pur)
      inf_sig_sparsesig <- sparsesig_aligned[[spn_id]][[c]][[p]] %>% as.data.frame()
      inf_sig_sigprofiler <- sigprof_aligned[[spn_id]][[c]][[p]] 
      sim_sig <- ground_truth_nested[[spn_id]][[c]][[p]]  %>% as.data.frame()
      
      sample_metrics_sparsesig <- per_sample_metrics(inferred_df = inf_sig_sparsesig, simulated_df = sim_sig) %>% 
        mutate(caller = 'SparseSignatures', spn = spn_id, coverage = cov, purity = pur)
      sample_metrics_sigprofiler <- per_sample_metrics(inferred_df = inf_sig_sigprofiler, simulated_df = sim_sig) %>% 
        mutate(caller = 'SigProfiler', spn = spn_id, coverage = cov, purity = pur)
      summary_sparsesig <- summary_stats(sample_metrics_sparsesig) %>% mutate(caller = 'SparseSignatures', spn = spn_id, coverage = cov, purity = pur)
      summary_sigprofiler <- summary_stats(sample_metrics_sigprofiler) %>% mutate(caller = 'SigProfiler', spn = spn_id, coverage = cov, purity = pur)
      
      cosine_mse_sigprofiler <- compute_cosine_mse(inferred = inf_sig_sigprofiler, simulated = sim_sig) %>% 
        mutate(caller = 'SigProfiler', spn = spn_id, coverage = cov, purity = pur)
      cosine_mse_sparsesig <- compute_cosine_mse(inferred = inf_sig_sparsesig, simulated = sim_sig) %>% 
        mutate(caller = 'SparseSignatures', spn = spn_id, coverage = cov, purity = pur)
      
      metrics_sample <- bind_rows(metrics_sample, sample_metrics_sparsesig, sample_metrics_sigprofiler)
      metrics_spn <- bind_rows(metrics_spn, summary_sparsesig, summary_sigprofiler)
      
      cosine_mse <- bind_rows(cosine_mse, cosine_mse_sigprofiler, cosine_mse_sparsesig)
    }
  }
}

col_spn <- c('SPN01' = 'steelblue', 'SPN02' ='seagreen', 'SPN03' ='goldenrod', 
             'SPN04' ='coral', 'SPN06' ='palevioletred', 'SPN07' ='indianred3')

### Generate final plots ###
cosine_mse_plot <- cosine_mse %>% ggplot() +
  geom_boxplot(aes(x = caller, y = cosine, col = spn)) +
  scale_color_manual(values = col_spn) +
  facet_grid(coverage~purity) + 
  ggtitle(paste0('Simulated vs inferred exposures per sample')) + 
  theme_bw() + 

cosine_mse %>% ggplot() +
  geom_boxplot(aes(x = caller, y = mse, col = spn)) +
  scale_color_manual(values = col_spn) +
  facet_grid(coverage~purity) + 
  theme_bw() + 
  plot_layout(guides = 'collect')

metric_plot <- metrics_sample %>% 
  ggplot() +
  geom_col(aes(x = caller, y = recall, fill = spn),position=position_dodge()) +
  scale_fill_manual(values = col_spn) +
  facet_grid(coverage~purity) + 
  ggtitle(paste0('Simulated vs inferred SBS signature')) +
  theme_bw() + 
  
  metrics_sample %>% 
  ggplot() +
  geom_col(aes(x = caller, y = precision, fill = spn),position=position_dodge()) +
  scale_fill_manual(values = col_spn) +
  facet_grid(coverage~purity) + 
  theme_bw() + 
  plot_layout(guides = 'collect')

wrap_plots(metric_plot, cosine_mse_plot, nrow = 2) 

# # Validate Mutation Counts
# inf_sparse_sig <- tumourevo_signature_res$SPN03$coverage_50$purity_0.9$SparseSig_mut_counts
# 
# inf_sig_profiler <- tumourevo_signature_res$SPN03$coverage_50$purity_0.9$SigProfiler_COSMIC_signatures
# exp_sig_profiler <- tumourevo_signature_res$SPN03$coverage_50$purity_0.9$SigProfiler_COSMIC_exposure
# 
# # recostruct signature profiler for sigprofiler ####################
# signature_matrix <- inf_sig_profiler
# rownames(signature_matrix) <- signature_matrix$MutationType
# signature_matrix$MutationType <- NULL  # drop the column for matrix math
# exposures <- exp_sig_profiler
# rownames(exposures) <- exposures$Samples
# exposures$Samples <- NULL
# signature_matrix <- signature_matrix[, colnames(exposures)]
# reconstructed_counts <- as.matrix(signature_matrix) %*% t(as.matrix(exposures))
# reconstructed_df <- data.frame(
#   MutationType = rownames(reconstructed_counts),
#   reconstructed_counts,
#   row.names = NULL
# )
# reconstructed_df[,-1] <- round(reconstructed_df[,-1])
# head(reconstructed_df)
# ######################################################################
# 
# sim_sig <- read.table('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN03/tumourevo/50x_0.9p_mutect2_ascat/signature_deconvolution/SigProfiler/SCOUT/results/SBS96/SBS96/Samples.txt',header = T)
# head(sim_sig)
# 
# normalize_profile <- function(df) {
#   mat <- as.matrix(df[,-1])
#   mat <- apply(mat, 2, function(x) x / sum(x))
#   df_norm <- data.frame(MutationType = df$MutationType, mat)
#   return(df_norm)
# }
# 
# reconstructed_rel <- normalize_profile(reconstructed_df)
# sim_rel <- normalize_profile(sim_sig)
# 
# cosine_similarity <- function(x, y) {
#   sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
# }
# 
# mse <- function(x, y) {
#   mean((x - y)^2)
# }
# 
# reconstructed_rel <- reconstructed_rel[order(reconstructed_rel$MutationType), ]
# sim_rel <- sim_rel[order(sim_rel$MutationType), ]
# 
# recon_mat <- as.matrix(reconstructed_rel[,-1])
# sim_mat   <- as.matrix(sim_rel[,-1])
# 
# n_samples <- min(ncol(recon_mat), ncol(sim_mat))
# results <- data.frame(
#   Sample_Reconstructed = colnames(recon_mat)[1:n_samples],
#   Sample_Simulated     = colnames(sim_mat)[1:n_samples],
#   Cosine_Similarity     = NA,
#   MSE                   = NA
# )
# 
# for (i in 1:n_samples) {
#   x <- recon_mat[, i]
#   y <- sim_mat[, i]
#   
#   results$Cosine_Similarity[i] <- cosine_similarity(x, y)
#   results$MSE[i] <- mse(x, y)
# }
# 
# print(results)
# 
# 
