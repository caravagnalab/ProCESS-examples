library(ProCESS)
# github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples'
source(file.path(github_path, "getters/process_getters.R"))
source(file.path(github_path, "getters/tumourevo_getters.R"))
source(file.path(github_path, "validation/Subclonal_deconvolution/utils_tables.R"))

my_get_tumourevo_subclonal <- function(
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
  MAIN_PATH <- sub("tumourevo", "tumourevo_old", MAIN_PATH)
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

