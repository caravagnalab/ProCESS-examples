
process_seq_results <- function(gt_path, spn, purity, coverage, chromosome, base_path, outdir) {
  # Construct the output folder path
  combination = paste0(coverage, "x_", purity, "p")
  #combination <- paste0("purity_", purity, "_coverage_", coverage, "x")
  folder_path <- file.path(outdir, spn, combination, "process")
  dir.create(folder_path, recursive = TRUE, showWarnings = TRUE)
  
  # Load sequencing results
  message("Reading sequencing results...")
  seq_res <- readRDS(gt_path)
  
  # Filter out germinal mutations
  message("Filtering out germinal mutations...")
  seq_res <- seq_res %>% dplyr::filter(classes!="germinal")
  
  # Extract sample names from column headers
  samples <- str_replace(colnames(seq_res)[grepl(".VAF", colnames(seq_res))], ".VAF", "")
  
  # Iterate through each sample and process mutations
  for (sample in samples) {
    message(paste0("Parsing sample ", sample, "..."))
    sample_path <- file.path(folder_path, sample)
    dir.create(sample_path, recursive = TRUE, showWarnings = TRUE)
    
    for (mutation in c("SNV", "INDEL")) {
      message(paste0("Parsing ", mutation, " mutations..."))
      mut_path <- file.path(sample_path, mutation)
      dir.create(mut_path, recursive = TRUE, showWarnings = TRUE)
      
      message(paste0("Parsing chr", chromosome, "..."))
      process_sample_mutation_chromosome(sample, mutation, chromosome, seq_res, mut_path)
    }
  }
}

process_sample_mutation_chromosome <- function(sample, mutation, chromosome, seq_res, mut_path) {
  # Filter sequencing results for the current chromosome
  sample_data <- seq_res %>% dplyr::filter(chr == chromosome)
  
  # Further filter based on mutation type
  if (mutation == "SNV") {
    sample_data <- sample_data %>% dplyr::filter(alt %in% c("A", "C", "T", "G"))
  } else if (mutation == "INDEL") {
    sample_data <- sample_data %>% dplyr::filter(!(alt %in% c("A", "C", "T", "G")))
  } else {
    stop("Mutation type not recognized")
  }
  
  # Convert sequencing data to long format if necessary
  if ("sample_name" %in% colnames(sample_data)) {
    seq_res_long <- sample_data
  } else {
    seq_res_long <- ProCESS::seq_to_long(sample_data)
  }
  
  # Filter and annotate mutation data
  seq_res_long <- seq_res_long %>%
    dplyr::filter(sample_name == sample) %>%
    dplyr::filter(NV != 0) %>% 
    dplyr::mutate(mutationID = paste0("chr", chr, ":", from, ":", ref, ":", alt))
  
  # Save the processed mutation data
  file_name <- file.path(mut_path, paste0("chr", chromosome, ".rds"))
  saveRDS(seq_res_long, file_name)
}

# 
# get_seq_res = function(gt_path, sample_id, chromosome, mut_type) {
#   # Read ground truth
#   message("reading seq res")
#   seq_res = readRDS(gt_path)
#   
#   message("keeping desired chromosome")
#   seq_res = seq_res %>% 
#     dplyr::filter(classes!="germinal") %>% 
#     dplyr::filter(chr == str_replace(chromosome, "chr", ""))
#   
#   if (mut_type == "SNV") {
#     seq_res = seq_res %>% dplyr::filter(alt %in% c("A", "C", "T", "G"))
#   } else if (mut_type == "INDEL") {
#     seq_res = seq_res %>% dplyr::filter(!(alt %in% c("A", "C", "T", "G")))
#   } else {
#     stop("mut_type parameter not recognized")
#   }
#   
#   if ("sample_name" %in% colnames(seq_res)){
#     seq_res_long <- seq_res
#   } else {
#     seq_res_long <- seq_to_long(seq_res)
#   }
#   
#   message("keeping SNVs and desired sample")
#   seq_res_long %>%
#     dplyr::filter(sample_name==sample_id) %>%
#     dplyr::filter(NV!=0) %>% 
#     dplyr::mutate(mutationID=paste0("chr",chr,":",from,":",ref,":",alt))
# }
