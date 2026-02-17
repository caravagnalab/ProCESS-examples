tool = "process_univariate_w_private"
source("/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/generate_table_main.R")
source('/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/get_mrca_branches.R')
get_cell_id = function(mutation_object) {
  tryCatch(
    expr = { phylo_forest$get_first_occurrences(mutation_object)[[1]] },
    error = function(e) return(NA)
  )
}
# source("~/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/generate_table_main.R")
# spn = 'SPN01'
# coverage=100
# purity=0.3

subforest_path = file.path("/orfeo/cephfs/scratch/cdslab/shared/SCOUT", spn, "process")
cli::cli_text("Generating {tool} table for {spn} and simulation {simulation_id}")

mut_process = get_mutations(spn=spn, type="tumour", coverage=coverage, purity=purity)
seq_results = readRDS(mut_process)

mut_process = seq_results %>%
  mutate(mutation_id=paste0(spn, ":", chr, ":", chr_pos,  ":", alt),
         is_driver_process=classes=="driver")

## Sticks #####

sample_forest = load_sample_forest(get_sample_forest(spn))

sample_names = (sample_forest$get_samples_info())$name

final_table_complete = data.frame()
# Other samples ####
# for(j in length(sample_names)-1){
  # i = j+1

# if(!(spn %in% c('SPN03', 'SPN04'))){
phylo_forest = load_phylogenetic_forest(get_phylo_forest(spn=spn))
mrca_branches_complete = get_mrca_branches(phylo_forest)
# }
sample_names = sort(sample_names)

i = 1
for(i in 1:length(sample_names)){

  sample = sample_names[[i]]
  print(sample)
  test_sample_forest = sample_forest$get_subforest_for(sample) #"SPN04_1.1"
  
  if(spn %in% c('SPN03', 'SPN04')){
    test_forest = load_phylogenetic_forest(paste0(subforest_path, "/subforest_", spn, "_", stringr::str_remove(sample, ".*_"), '.sff'))
    }else{
    # phylo_forest = load_phylogenetic_forest(get_phylo_forest(spn=spn))
    test_forest = phylo_forest$get_subforest_for(sample)
  }
  # plot_forest(test_sample_forest) %>% annotate_forest(test_sample_forest)
  # mut_process = mut_process_all %>% filter(sample_id == sample)
  
  sticks <- test_forest$get_sticks()
  relevant_branches = get_relevant_branches(test_forest)
  # plot_sticks(test_sample_forest, relevant_branches)
  
  if (is.null(relevant_branches)) {
    relevant_branches=data.frame(
      'cell_id' = 0,
      'birth_time' = 0,
      'mutant' = 'Clone 1',
      'label' = 'Truncal'
    )
  }
  # ----------- New branches ----------- #
  # mrca_branches = get_mrca_branches(test_forest)
  
  mrca_branches = keep(mrca_branches_complete, ~ .x$sample == sample)
  
  if(length(mrca_branches)>0){
    for(i in 1:length(mrca_branches)){
      tmp = mrca_branches[[i]]
      ids = tmp$ids
      target_label = paste0(tmp$mutant, '_', tmp$sample)
      print(target_label)
      for (id in ids) {
        if (id %in% relevant_branches$cell_id) { # if it is already present in relevant_branches, replace the label with the new one
          print(paste0(id, ' Already present'))
          relevant_branches <- relevant_branches %>%
            mutate(label = if_else(cell_id == id, target_label, label))
        } else {
          print('Not already present')
          relevant_branches <- relevant_branches %>%
            add_row(
              cell_id    = id,
              ancestor   = NA_integer_,
              mutant     = NA_character_,
              epistate   = NA_character_,
              sample     = NA_character_,
              birth_time = NA_real_,
              label      = target_label
            )
        }
      }
    }
  }
  
  # plot_sticks(test_sample_forest, relevant_branches)
  # ----------- New branches ----------- #
  
  mut_process_with_clusterid = mut_process %>% 
    filter(.data[[paste0(sample, ".VAF")]] > 0) %>%
    filter(classes != "germinal") %>%
    rowwise() %>%
    mutate(cell_id=get_cell_id(Mutation(chr, chr_pos, ref, alt))) %>%
    ungroup() %>% 
    left_join(relevant_branches) %>% 
    ungroup() %>%
    select(cell_id, mutation_id, causes, is_driver_process, label, contains(".VAF")) %>%
    pivot_longer(
      cols=ends_with(".VAF"),
      names_to="sample_id",
      names_pattern="(.*)\\.VAF", # remove matching text "VAF" from the start of each variable name
      values_to="vaf_process" # this is the VAF!
    ) %>%
    rename(cluster_id_process=label)
  
  print(paste0('unique(mut_process_with_clusterid$cluster_id_process) ', unique(mut_process_with_clusterid$cluster_id_process)))
  
  # all sampled cells
  sample_cell_ids = test_sample_forest$get_nodes() %>%
    filter(!is.na(sample)) %>% 
    group_by(sample) %>%
    summarise(sample_cell_ids=list(unique(cell_id)),
              Ntot=length(unique(cell_id)),
              .groups="drop") %>% 
    rename(sample_id=sample)
  
  print('sample_cell_ids')
  
  # all sequenced mutations and cell ids
  mut_cells = test_forest$get_sampled_cell_mutations() %>%
    mutate(mutation_id=paste0(spn, ":", chr, ":", chr_pos,  ":", alt)) %>% 
    filter(mutation_id %in% mut_process_with_clusterid$mutation_id) %>% 
    group_by(mutation_id) %>%
    summarise(cells_w_mutation=list(unique(cell_id)),
              .groups="drop")
  
  print('mut_cells')
  
  # keep only mutations with non-zero frequency
  final_table = mut_process_with_clusterid %>%
    filter(vaf_process > 0) %>%
    left_join(mut_cells) %>% 
    left_join(sample_cell_ids) %>%
    rowwise() %>% 
    mutate(n_cells=length(intersect(cells_w_mutation, sample_cell_ids)),
           ccf_process=n_cells / Ntot) %>% 
    ungroup() %>% 
    group_by(cluster_id_process, sample_id) %>%
    # group_by(cluster_id_tool) %>%
    mutate(ccf_process=mean(ccf_process, na.rm=TRUE)) %>%
    ungroup()
  
  print('final_table')
  print(unique(final_table$sample_id))
  
  final_table = final_table %>% 
    # mutate(is_clonal_process=ifelse(cluster_id_process==get_clonal_cluster_process(final_table), TRUE, FALSE)) %>% 
    # select(-cell_id, -cells_w_mutation, -sample_cell_ids, -Ntot, -n_cells) %>% 
    select(-cells_w_mutation, -sample_cell_ids, -Ntot, -n_cells) %>% 
    mutate(n_clones_process=length(unique(cluster_id_process))) %>% 
    group_by(cluster_id_process) %>% 
    mutate(n_mutations_process=length(unique(mutation_id))) %>% 
    ungroup() %>% 
    mutate(#sample_id=paste0(spn, "_", sample_id),
           driver_label_process=ifelse(is_driver_process, mutation_id, NA)) %>% 
    
    mutate(patient_id=spn, coverage=coverage, purity=purity) %>% 
    select(patient_id, sample_id, coverage, purity, mutation_id, dplyr::everything())%>% 
    filter(sample_id==sample) %>% # metti?
    mutate(sample_id=paste0(spn, "_", sample_id))
      
    print('final_table2')
    print(unique(final_table$sample_id))
    final_table_complete = bind_rows(final_table_complete, final_table)
    
}
tool = 'process_univariate_w_private'
out_path = get_table_path(save_path, tool, spn, simulation_id)
cli::cli_text("Saving {tool} table for {spn} and simulation {simulation_id} in {out_path}")

saveRDS(final_table_complete, out_path)

