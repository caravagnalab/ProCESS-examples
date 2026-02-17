tool = "process"
source("/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/generate_table_main.R")
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

# tool = "process"
out_path = get_table_path(save_path, tool, spn, simulation_id)

cli::cli_text("Generating {tool} table for {spn} and simulation {simulation_id}")

mut_process = get_mutations(spn=spn, type="tumour", coverage=coverage, purity=purity)
seq_results = readRDS(mut_process)

mut_process = seq_results %>%
  mutate(mutation_id=paste0(spn, ":", chr, ":", chr_pos,  ":", alt),
         is_driver_process=classes=="driver")

## Sticks #####

sample_forest = load_sample_forest(get_sample_forest(spn))
phylo_forest = load_phylogenetic_forest(get_phylo_forest(spn=spn))

sticks <- phylo_forest$get_sticks()
relevant_branches = get_relevant_branches(phylo_forest)

mut_process_with_clusterid = mut_process %>% 
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

# all sampled cells
sample_cell_ids = sample_forest$get_nodes() %>%
  filter(!is.na(sample)) %>% 
  group_by(sample) %>%
  summarise(sample_cell_ids=list(unique(cell_id)),
            Ntot=length(unique(cell_id)),
            .groups="drop") %>% 
  rename(sample_id=sample)


# all sequenced mutations and cell ids
mut_cells = phylo_forest$get_sampled_cell_mutations() %>%
  mutate(mutation_id=paste0(spn, ":", chr, ":", chr_pos,  ":", alt)) %>% 
  filter(mutation_id %in% mut_process_with_clusterid$mutation_id) %>% 
  group_by(mutation_id) %>%
  summarise(cells_w_mutation=list(unique(cell_id)),
            .groups="drop")

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

final_table = final_table %>% 
  # mutate(is_clonal_process=ifelse(cluster_id_process==get_clonal_cluster_process(final_table), TRUE, FALSE)) %>% 
  select(-cell_id, -cells_w_mutation, -sample_cell_ids, -Ntot, -n_cells) %>% 
  mutate(n_clones_process=length(unique(cluster_id_process))) %>% 
  group_by(cluster_id_process) %>% 
  mutate(n_mutations_process=length(unique(mutation_id))) %>% 
  ungroup() %>% 
  mutate(sample_id=paste0(spn, "_", sample_id),
         driver_label_process=ifelse(is_driver_process, mutation_id, NA)) %>% 
  
  mutate(patient_id=spn, coverage=coverage, purity=purity) %>% 
  select(patient_id, sample_id, coverage, purity, mutation_id, dplyr::everything())


cli::cli_text("Saving {tool} table for {spn} and simulation {simulation_id} in {out_path}")

saveRDS(final_table, out_path)

