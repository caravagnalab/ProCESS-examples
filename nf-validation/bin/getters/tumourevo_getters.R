################################################################################
# UTILS
################################################################################
# return corresponding matched directory based on spn, coverage, purity, vcf_caller and cna_caller
dir_getter <- function(
    spn, 
    coverage, 
    purity, 
    vcf_caller, 
    cna_caller, 
    base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
) {
  
  vcf_caller_list <- c("mutect2", "strelka")
  cna_caller_list <- c("ascat")
  
  if (!is.numeric(coverage) || coverage %% 1 != 0) {
    stop("Error: 'coverage' must be integer") # coverage : integer
  }
  if (!is.numeric(purity) || (purity > 1 || purity < 0)) {
    stop("Error: 'purity' must be a number between 0 and 1") # purity : float [0, 1]
  }
  if (!(vcf_caller %in% vcf_caller_list)) {
    stop("Error: 'vcf_caller' is not valid")
  }
  if (!(cna_caller %in% cna_caller_list)) {
    stop("Error: 'typecna_caller' is not valid")
  }
  
  # EXTEND PATH
  MAIN_PATH <- file.path(base_path, spn, "tumourevo")
  
  ptr <- paste(paste0(coverage, "x"), paste0(purity, "p"), vcf_caller, cna_caller, sep = "_")
  all_dirs <- list.dirs(path = MAIN_PATH, full.names = T, recursive = FALSE)
  matched_dir <- all_dirs[grepl(ptr, basename(all_dirs))]
  
  return(matched_dir)
}



# given a directory path as input, returns the named list of files residing in that directory
get_named_file_list <- function(dir_path) {
  files <- list.files(dir_path, full.names = TRUE)
  files <- files[!file.info(files)$isdir]
  file_names <- basename(files)
  names(files) <- tools::file_path_sans_ext(file_names)
  return(as.list(files))
}



mobster_named_list <- function(input_list) {
  output_list <- list()
  
  for (i in seq_along(input_list)) {
    file_path <- input_list[[i]]
    file_name <- basename(file_path)
    
    # Extract extension
    file_ext <- tools::file_ext(file_name)
    
    # Get base name without extension
    name_without_ext <- tools::file_path_sans_ext(file_name)
    
    # Strip everything before "mobsterh" or "REPORT"
    simplified_key <- sub(".*?(mobsterh.*|REPORT.*)", "\\1", name_without_ext)
    
    # Make key unique if already present
    #if (simplified_key %in% names(output_list)) {
    simplified_key <- paste0(simplified_key, "_", file_ext)
    #}
    
    output_list[[simplified_key]] <- file_path
  }
  
  return(output_list)
}


ctree_named_list <- function(input_list) {
  output_list <- list()
  
  for (i in seq_along(input_list)) {
    file_path <- input_list[[i]]
    file_name <- basename(file_path)
    
    # Extract extension
    file_ext <- tools::file_ext(file_name)
    
    # Get base name without extension
    name_without_ext <- tools::file_path_sans_ext(file_name)
    
    # Strip everything before "mobsterh" or "REPORT"
    simplified_key <- sub(".*?(ctree*|REPORT.*)", "\\1", name_without_ext)
    
    # Make key unique if already present
    #if (simplified_key %in% names(output_list)) {
    simplified_key <- paste0(simplified_key, "_", file_ext)
    #}
    
    output_list[[simplified_key]] <- file_path
  }
  
  return(output_list)
}



pyclone_input_list <- function(input_list) {
  output_list <- list()
  
  for (i in seq_along(input_list)) {
    file_path <- input_list[[i]]
    file_name <- basename(file_path)
    
    # Remove "SCOUT_SPNxx_" prefix (xx can be any number)
    name_without_prefix <- sub("^SCOUT_SPN\\d+_", "", file_name)
    
    # Replace dots in extensions with underscores
    name_cleaned <- gsub("\\.", "_", name_without_prefix)
    
    output_list[[name_cleaned]] <- file_path
  }
  
  return(output_list)
}


viber_input_list <- function(input_list) {
  output_list <- list()
  
  for (i in seq_along(input_list)) {
    file_path <- input_list[[i]]
    file_name <- basename(file_path)
    
    # Remove "SCOUT_SPNxx_" prefix (xx can be any number)
    name_without_prefix <- sub("^SCOUT_SPN\\d+_", "", file_name)
    
    # Replace dots in extensions with underscores
    name_cleaned <- gsub("\\.", "_", name_without_prefix)
    
    output_list[[name_cleaned]] <- file_path
  }
  
  return(output_list)
}

qc_input_list <- function(input_list) {
  output_list <- list()
  
  for (i in seq_along(input_list)) {
    file_path <- input_list[[i]]
    file_name <- basename(file_path)
    
    # Remove "SCOUT_SPNxx_" prefix (xx can be any number)
    name_without_prefix <- sub("^SCOUT_SPN\\d+_SPN\\d+_SPN\\d+_\\d+\\.\\d+_", "", file_name)
    
    # Replace dots in extensions with underscores
    name_cleaned <- gsub("\\.", "_", name_without_prefix)
    
    output_list[[name_cleaned]] <- file_path
  }
  
  return(output_list)
}


################################################################################

################################################################################
# DRIVERS
################################################################################
get_tumourevo_driver <- function(
    spn, 
    coverage, 
    purity, 
    vcf_caller, 
    cna_caller, 
    sample, 
    base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
) {
  
  MAIN_PATH <- dir_getter(
    spn, 
    coverage, 
    purity, 
    vcf_caller, 
    cna_caller, 
    base_path
  )
  
  MAIN_PATH <- file.path(MAIN_PATH, "driver_annotation/annotate_driver/SCOUT", spn)
  
  output <- list.files(MAIN_PATH, pattern = paste0("\\", sample, "_driver.rds$"), full.names = T, recursive=T)
  
  return(output)
}



################################################################################
# SUBCLONAL
################################################################################
get_tumourevo_subclonal <- function(
    spn, 
    coverage, 
    purity, 
    tool, 
    vcf_caller, 
    cna_caller, 
    sample, 
    base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
) {
  
  # quality control
  tool_list <- c("mobster", "pyclonevi", "ctree", "viber")
  if (!(tool %in% tool_list)) {
    stop("ERROR: wrong tool name!")
  }
  
  MAIN_PATH <- dir_getter(
    spn, 
    coverage, 
    purity, 
    vcf_caller, 
    cna_caller, 
    base_path
  )
  
  MAIN_PATH <- file.path(MAIN_PATH, "subclonal_deconvolution", tool, "SCOUT", spn)
  
  if (tool=="ctree") {
    if (sample != spn) {
      MAIN_PATH <- file.path(MAIN_PATH, paste0(spn, "_", sample))
    }
    output <- get_named_file_list(MAIN_PATH)
    output <- ctree_named_list(output)
    
  } else if (tool=="mobster") {
    MAIN_PATH <- file.path(MAIN_PATH, paste0(spn, "_", sample))
    output <- get_named_file_list(MAIN_PATH)
    output <- mobster_named_list(output)
  } else if (tool=="pyclonevi") {
    output <- get_named_file_list(MAIN_PATH)
    output <- pyclone_input_list(output)
  } else if (tool=="viber") {
    output <- get_named_file_list(MAIN_PATH)
    output <- viber_input_list(output)
  }
  return(output)
}


################################################################################
# QC
################################################################################
get_tumourevo_qc <- function(
    spn, 
    coverage, 
    purity, 
    tool, 
    vcf_caller, 
    cna_caller, 
    sample, 
    base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
) {
  
  # quality control
  tool_list <- c("CNAqc", "join_CNAqc", "tinc")
  if (!(tool %in% tool_list)) {
    stop("ERROR: wrong tool name!")
  }
  
  MAIN_PATH <- dir_getter(
    spn, 
    coverage, 
    purity, 
    vcf_caller, 
    cna_caller, 
    base_path
  )
  
  MAIN_PATH <- file.path(MAIN_PATH, "QC", tool, "SCOUT", spn)
  
  if (tool != "join_CNAqc") {
    all_entries <- list.dirs(MAIN_PATH, full.names = FALSE, recursive = FALSE)
    matching_dir <- all_entries[grepl(sample, all_entries)]
    MAIN_PATH <- file.path(MAIN_PATH, matching_dir)
    output <- get_named_file_list(MAIN_PATH)
    output <- qc_input_list(output)
  } else {
    output <- get_named_file_list(MAIN_PATH)
    output <- viber_input_list(output)
  }

  
  return(output)
  
}




################################################################################
# SIGNATURES
################################################################################
get_tumourevo_signatures <- function(
    spn, 
    coverage, 
    purity, 
    tool, 
    context, 
    vcf_caller, 
    cna_caller, 
    base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
) {
  
  # quality control
  tool_list <- c("SigProfiler", "SparseSignatures")
  if (!(tool %in% tool_list)) {
    stop("ERROR: wrong tool name!")
  }
  
  MAIN_PATH <- dir_getter(
    spn, 
    coverage, 
    purity, 
    vcf_caller, 
    cna_caller, 
    base_path
  )
  
  MAIN_PATH <- file.path(MAIN_PATH, "signature_deconvolution")
  
  
  
  if (tool == "SigProfiler") {
    
    MAIN_PATH <- file.path(MAIN_PATH, tool, "SCOUT", "results", context, context, "Suggested_Solution")
    
    COSMIC <- file.path(MAIN_PATH, paste0("COSMIC_", context, "_Decomposed_Solution"))
    COSMIC_exposure <- paste0(COSMIC, "/Activities/COSMIC_", context, "_Activities.txt")
    COSMIC_signatures <- paste0(COSMIC, "/Signatures/COSMIC_", context, "_Signatures.txt")
    
    denovo <- file.path(MAIN_PATH, paste0(context, "_De-Novo_Solution"))
    denovo_exposure <- paste0(denovo, "/Activities/", context, "_De-Novo_Activities_refit.txt")
    denovo_signatures <- paste0(denovo, "/Signatures/", context, "_De-Novo_Signatures.txt")
    
    return(list(
      COSMIC_exposure=COSMIC_exposure, 
      COSMIC_signatures=COSMIC_signatures, 
      denovo_exposure=denovo_exposure, 
      denovo_signatures=denovo_signatures
    ))
    
  } else if (tool == "SparseSignatures") {
    
    MAIN_PATH <- file.path(MAIN_PATH, tool, "SCOUT")
    output <- list() 
    output[['nmf_Lasso_out']] <- file.path(MAIN_PATH, "SCOUT_nmf_Lasso_out.rds")
    output[['cv_means_mse']] <- file.path(MAIN_PATH, "SCOUT_cv_means_mse.rds")
    output[['best_params_config']] <- file.path(MAIN_PATH, "SCOUT_best_params_config.rds")
    output[['mut_counts']] <- file.path(MAIN_PATH, "SCOUT_mut_counts.rds")
    
    return(output)
  } else {
    stop("ERROR: wrong tool name!")
  }
  
}

