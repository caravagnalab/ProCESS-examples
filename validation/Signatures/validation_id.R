setwd(setwd("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/")

pkgs <- c("ProCESS", "tidyverse", "ggplot2", "caret", "ggtext", "reshape2", "lsa", "Metrics", 
"ggalluvial", "patchwork")
sapply(pkgs, require, character.only = TRUE)

#setwd("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/")

source("getters/process_getters.R")
source("getters/tumourevo_getters.R")
source("validation/signatures/utils/utils_getters.R")
source("validation/signatures/utils/utils_validation.R")
source("validation/signatures/utils/utils.R")
source("validation/signatures/utils/utils_plots.R")


##### Use getters to get a pairs of RACES and tumourevo data #####

### Get ProCESS exposure data ###

base_path <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
spn <- c("SPN01", "SPN03", "SPN04")
context <- "ID"

process_exposures_list <- list()

for (spn_id in spn) {
  exposure_id <- tryCatch({
    get_exposures_by_context(spn = spn_id, base_path = base_path, context = context)
  }, error = function(e) {
    warning("Failed for ", spn_id, ": ", e$message)
    return(NULL)
  })

  process_exposures_list[[spn_id]] <- exposure_id
}

### Load tumourevo signature data ###

base_path <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
spn_list <- c("SPN01", "SPN03", "SPN04")
coverage_list <- c(50)
purity_list <- c(0.9)
dataset <- "SCOUT"
context <- 'ID83'
vcf_caller <- "mutect2"
cna_caller <- "ascat"

tumourevo_signature_res <- list()

# Iterate over all combinations
for (spn in spn_list) {
  tumourevo_signature_res[[spn]] <- list()
  
  for (cov in coverage_list) {
    tumourevo_signature_res[[spn]][[paste0("coverage_", cov)]] <- list()
    
    for (pur in purity_list) {
      message("Processing ", spn, " | Coverage: ", cov, " | Purity: ", pur)
      
      result <- tryCatch({
        # Get SigProfiler paths
        sigprofiler <- get_tumourevo_signatures(
          spn = spn,
          coverage = cov,
          purity = pur,
          vcf_caller = vcf_caller,
          cna_caller = cna_caller,
          tool = "SigProfiler",
          context = context,  # ID83!
          base_path = base_path
        )
        
        # Combine and load paths
        paths <- c(
          sigprofiler$COSMIC_exposure,
          sigprofiler$COSMIC_signatures,
          sigprofiler$denovo_exposure,
          sigprofiler$denovo_signatures
        )
        
        data <- load_signature_data(paths)
        
        names(data) <- c(
          "SigProfiler_COSMIC_exposure",
          "SigProfiler_COSMIC_signatures",
          "SigProfiler_denovo_exposure",
          "SigProfiler_denovo_signatures"
        )
        
        data
      }, error = function(e) {
        message("Failed for ", spn, " | cov: ", cov, " | pur: ", pur)
        message("Reason: ", e$message)
        NULL
      })
      
      # Store result if successful
      if (!is.null(result)) {
        tumourevo_signature_res[[spn]][[paste0("coverage_", cov)]][[paste0("purity_", pur)]] <- result
      }
    }
  }
}

### Summarize Signatures across combinations ###

spn_list <- c("SPN01", "SPN03", "SPN04")
coverage <- c(50)
purity_values <- c(0.9)

ground_truth <- list()

for (spn in spn_list) {
  if (spn %in% names(process_exposures_list)) {
    gt <- process_exposures_list[[spn]] %>%
      tibble::column_to_rownames("Sample_ID") %>%
      as.matrix()
    
    ground_truth[[spn]] <- gt
  } else {
    warning(paste("Missing ground truth for", spn))
  }
}

sigprof_aligned <- align_sigprofiler_res(tumourevo_signature_res)

ground_truth_nested <- list()

for (spn in names(ground_truth)) {
  if (is.null(ground_truth_nested[[spn]])) {
    ground_truth_nested[[spn]] <- list()
  }
  
  for (cov in coverage) {
    coverage_key <- paste0("coverage_", cov)
    
    if (is.null(ground_truth_nested[[spn]][[coverage_key]])) {
      ground_truth_nested[[spn]][[coverage_key]] <- list()
    }
    
    for (purity in purity_values) {
      purity_key <- paste0("purity_", purity)
      ground_truth_nested[[spn]][[coverage_key]][[purity_key]] <- ground_truth[[spn]]
    }
  }
}


### Compare exposures of estimated and true signatures  ###

sankey_df <- prepare_sankey_data_id(ground_truth_nested, sigprof_aligned)

saveRDS(sankey_df, file = "validation/sankey_signatures_id.rds")

# Plot sankey for each SPN

sankey_df <- sankey_df %>%
  dplyr::mutate(
    Coverage = as.numeric(gsub("coverage_", "", Coverage)),
    Purity = as.numeric(gsub("purity_", "", Purity))
  ) %>% 
  dplyr::filter(Coverage == 50, Purity == 0.9, SPN == "SPN01")


spns <- c("SPN01")
sankey_plots <- lapply(spns, function(spn) generate_sankey_id(spn, cov = 50, pur = 0.9))

wrapped_sankey <- wrap_plots(sankey_plots, ncol = 1)
wrapped_sankey

saveRDS(wrapped_sankey, file = "validation/sankey_SPN01_50x_0.9p_id.rds")
ggsave("validation/sankey_SPN01_50x_0.9p_id.pdf", plot = wrapped_sankey, width = 12, height = 5, units = "in")

