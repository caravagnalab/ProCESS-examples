#!/usr/bin/env Rscript
get_sample_names <- function(spn, base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"){
  search_path <- file.path(base_path, spn, "process", spn)
  name_files <- dir(search_path, pattern = "\\.rff$")
  sample_names <- gsub('.{4}$', '', name_files)
  return(sample_names)
}

get_process_cna <- function(spn, sample=NULL, base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"){
    if (!is.character(spn)){
        stop("Error: spn must be character")
    }
    if (!is.character(base_path)){
       stop("Error: base_path must be character") 
    }
    # add error handling for sample
    if (is.character(sample)) {
        search_path <- file.path(base_path, spn, "process/cna_data", paste0(sample, "_cna.rds"))
        if (!file.exists(search_path)){
          stop("Error: Could not find cna file")
        } else {
          return(search_path)
        }
    } else if (is.null(sample)) {
      samples <- get_sample_names(spn)
      file_list <- list()
      for (i in 1:length(samples)){
        search_path <- file.path(base_path, spn, "process/cna_data", paste0(samples[i], "_cna.rds"))
        if (!file.exists(search_path)){
          stop("Error: Could not find cna file")
        } else {
          file_list[[samples[i]]] <- search_path
        }
      }
      return(file_list)
    }
}

get_process_gender <- function(spn, base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/") {
  search_path <- file.path(base_path, spn, "process", "subject_gender.txt")
  if (file.exists(search_path)) {
    con <- file(search_path,"r")
    first_line <- readLines(con,n=1)
    close(con)
    return(first_line)
    # check that this is a character and not a chr vector
  } else {
    stop("Error: subject_gender file not found.")
  }
}

get_phylo_forest <- function(spn, base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"){
  search_path <- file.path(base_path, spn, "process", "phylo_forest.sff")
  if (file.exists(search_path)) {
    return(search_path)
  } else {
    stop("Error: phylo_forest file not found.")
  }
}

get_sample_forest <- function(spn, base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"){
  search_path <- file.path(base_path, spn, "process", "sample_forest.sff")
  if (file.exists(search_path)) {
    return(search_path)
  } else {
    stop("Error: sample_forest file not found.")
  }
}

get_mutations <- function(spn = NA, 
                          base_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/",
                          coverage = NA,
                          purity = NA,
                          type = NA ){
  if (is.na(spn)){
    stop('spn name is required')
  }
  
  if (is.na(type)){
    stop('type is required: normal or tumour')
  }
  
  if (type == 'tumour'){
    if (is.na(coverage)){
      stop('coverage is required')
    }
    if (is.na(purity)){
        stop('purity is required')
    }
    path <- file.path(base_path, spn, paste0('sequencing/tumour/purity_', purity), paste0('data/mutations/seq_results_muts_merged_coverage_',coverage, 'x.rds'))
  } else{
    path <- file.path(base_path, spn, 'sequencing/normal/purity_1/data/mutations/seq_results_muts_merged_coverage_30x.rds')
  }
  if (!file.exists(path)){
    stop('The required file do not exist, check purity and coverage values')
  }
  return(path)
}
