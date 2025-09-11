library(dplyr)

process_seq_results <- function(spn, purity, coverage, chromosome, outdir,rds_path,sample_id) {
  # Construct the output folder path
  combination = paste0(coverage, "x_", purity, "p")

  # Load sequencing results
  message("Reading sequencing results...")
  seq_res <- readRDS(rds_path)
  
  # Filter out germinal mutations
  message("Filtering out germinal mutations...")
  seq_res <- seq_res %>% dplyr::filter(classes!="germinal")
  message(paste0("Parsing sample ", sample_id, "..."))

  mut_snv = process_sample_mutation_chromosome(sample_id,'SNV', chromosome, seq_res)
  mut_indel = process_sample_mutation_chromosome(sample_id,'INDEL', chromosome, seq_res)
  
  # Add CCF status
  sample_forest = ProCESS::load_sample_forest(get_sample_forest(spn))
  phylo_forest = ProCESS::load_phylogenetic_forest(get_phylo_forest(spn))
  
  sample_cell_ids = sample_forest$get_nodes() %>% 
    dplyr::filter(sample == sample_id) %>% 
    dplyr::pull(cell_id) %>% 
    unique()
  
  mut_df = phylo_forest$get_sampled_cell_mutations()
  mut_df = mut_df %>% dplyr::mutate(mutationID = paste0("chr", chr, ":", chr_pos, ":", ref ,":", alt))
  
  mut_snv$CCF = lapply(1:nrow(mut_snv), function(i) {
    cells_w_mutation = mut_df %>% dplyr::filter(mutationID == mut_snv[i,]$mutationID) %>% 
      dplyr::pull(cell_id)
    n = (cells_w_mutation %in% sample_cell_ids) %>% sum()
    Ntot = length(sample_cell_ids)
    n / Ntot
  }) %>% unlist()
  
  mut_indel$CCF = lapply(1:nrow(mut_indel), function(i) {
    cells_w_mutation = mut_df %>% dplyr::filter(mutationID == mut_indel[i,]$mutationID) %>% 
      dplyr::pull(cell_id)
    n = (cells_w_mutation %in% sample_cell_ids) %>% sum()
    Ntot = length(sample_cell_ids)
    n / Ntot
  }) %>% unlist()
  
  mut = list('SNV' = mut_snv, 'INDEL' = mut_indel)
  saveRDS(object = mut, file = paste0("chr", chromosome, ".rds"))
}

process_sample_mutation_chromosome <- function(sample, mutation, chromosome, seq_res) {
  # Filter sequencing results for the current chromosome
  sample_data <- seq_res %>% dplyr::filter(chr == chromosome)
  
  # Further filter based on mutation type
  if (mutation == "SNV") {
    sample_data <- sample_data %>% dplyr::filter(alt %in% c("A", "C", "T", "G"), ref %in% c("A", "C", "T", "G") )
  } else if (mutation == "INDEL") {
    sample_data <- sample_data %>% dplyr::filter(!(alt %in% c("A", "C", "T", "G")) | !(ref %in% c("A", "C", "T", "G")))
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
  
  return(seq_res_long)
}