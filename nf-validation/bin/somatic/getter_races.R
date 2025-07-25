

#rm(list = ls()) # clears objects from the workspace
library(dplyr)
library(ProCESS)


#===============================================================================
# AUXILIARY FUNCTIONS
#===============================================================================

# quality control the input parameters
#-------------------------------------------------------------------------------
qc_args <- function(name, coverage, purity, timepoint, sample_id, type) {
  if (!is.numeric(name) || name %% 1 != 0) {
    stop("Error: 'name' must be integer") # name : integer
  }
  if (!is.numeric(coverage) || coverage %% 1 != 0) {
    stop("Error: 'coverage' must be integer") # coverage : integer
  }
  if (!is.numeric(purity) || (purity > 1 || purity < 0)) {
    stop("Error: 'purity' must be a number between 0 and 1") # purity : float [0, 1]
  }
  if (!is.null(timepoint) && (!is.numeric(timepoint) || timepoint %% 1 != 0)) {
    stop("Error: 'timepoint' must be integer") # timepoint : integer
  }
  if (!is.null(sample_id) && (!is.numeric(sample_id) || sample_id %% 1 != 0)) {
    stop("Error: 'sample_id' must be integer") # sample_id : integer
  }
  if (!(type %in% c("snv", "cna", "phylo"))) {
    stop("Error: 'type' must be one of 'SNV', 'CNA' and 'phylo'")
  }
  return(TRUE)
}


# generate the pattern for string matching for each required scenario
#-------------------------------------------------------------------------------

gen_ptr <- function(MAIN_PATH, SPN_ID, coverage, purity, timepoint, sample_id, type) {
  
  if (type == "cna") {
    path <- paste0(MAIN_PATH, "SPN0", SPN_ID, "/races/cna_data/")
    timepoint_ph <- ifelse(test = is.null(timepoint), yes = "\\d+", no = timepoint)
    sample_id_ph <- ifelse(test = is.null(sample_id), yes = "\\d+", no = sample_id)
    ptr <- paste0("^SPN0", SPN_ID, "_", timepoint_ph, "\\.", sample_id_ph, "_cna\\.rds$")
    return( c(path, ptr) )
  }
  
  else if (type == "snv") {
    path <- paste0(MAIN_PATH, "SPN0", SPN_ID, "/races/purity_", purity, "/")
    coverage_ph <- ifelse(test = is.null(coverage), yes = "\\d+", no = coverage)
    #purity_ph <- ifelse(test = is.null(purity), yes = "\\d+", no = purity)
    ptr <- paste0("seq_results_muts_merged_coverage_", coverage_ph, "x.rds")
    return( c(path, ptr) )
  }
  
  else if (type == "phylo") {
    path <- paste0(MAIN_PATH, "SPN0", SPN_ID, "/races/")
    ptr <- paste0("phylo_forest.sff")
    return( c(path, ptr) )
  }
  
  else {
    stop("Error in 'gen_ptr' function")
  }
}



# return list of absolute file paths as a result of string matching
#-------------------------------------------------------------------------------
search_files <- function(directory, ptr) {
  files <- list.files(path = directory, pattern = ptr, full.names = TRUE, recursive = TRUE)
  return(files)
}


# return the file object given the file path
#-------------------------------------------------------------------------------
load_multiple_rds <- function(file_paths) {
  ext <- tolower( tail(strsplit(basename(file_paths), "\\.")[[1]], 1) )
  if (ext != "rds") {
    stop("Invalid file type!")
  }
  data_list <- lapply(file_paths, readRDS)
  names(data_list) <- basename(file_paths)  # Use file names as list names
  return(data_list)
}



#===============================================================================
# MAIN GETTER FUNCTION
#===============================================================================
# input:
# PATH      : (string) - abosolute path
# SPN_ID    : (integer) - SPN ID
# coverage  : (integer)
# purity    : (float)
# timepoint : (integer)
# sample_id : (integer)
# type      : ("cna", "snv", "phylo")
# output:
# list of file object
getter_races <- function(MAIN_PATH = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/", 
                         SPN_ID = 1, 
                         coverage = 100, 
                         purity = 0.6, 
                         timepoint = NULL, 
                         sample_id = NULL, 
                         type = "cna") {
  
  type <- tolower(type)
  
  tryCatch(
    {
      qc_args(SPN_ID, coverage, purity, timepoint, sample_id, type)
      
      #files_list <- search_files(directory=PATH, ptr)
      
      if ( type %in% c("snv", "cna") ) {
        x <- gen_ptr(MAIN_PATH, SPN_ID, coverage, purity, timepoint, sample_id, type)
        path <- x[1]
        ptr <- x[2]
        files_list <- search_files(path, ptr)
        if (length(files_list) != 0) {
          return(load_multiple_rds(files_list))
        }
      }
      else if (type=="phylo") {
        x <- gen_ptr(MAIN_PATH, SPN_ID, coverage, purity, timepoint, sample_id, type)
        #f <- paste0(PATH, "SPN0", SPN_ID, "/races/phylo_forest.sff")
        return(rRACES::load_phylogenetic_forest(paste0(x[1], x[2])))
      }
      
      # if nothing found
      else {
        message("NOT FOUND!")
        return(NULL)
      }
      
    },
    error = function(cond) {
      message(conditionMessage(cond))
    }, 
    warning = function(cond) {
      message(conditionMessage(cond))
    }, 
    finally = {
      #message("JOB DONE!")
    }
  )
}