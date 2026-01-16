get_mrca_branches = function(phylo_forest){
  
  mutants = phylo_forest$get_nodes() %>% filter(!is.na(sample)) %>% group_by(sample,mutant) %>% 
    summarize(n = length(cell_id)) %>% ungroup()
  
  ids = phylo_forest$get_nodes() %>% filter(!is.na(sample)) %>% as_tibble() %>% 
    dplyr::select(mutant,sample,cell_id)
  
  # muts = phylo_forest$get_sampled_cell_mutations() %>% as_tibble() %>% 
  #   filter(class != "pre-neoplastic", type != 'indel')
  
  muts = phylo_forest$get_sampled_cell_mutations() %>% as_tibble() %>% 
    filter(class != "pre-neoplastic")
  
  tab = lapply(1:nrow(mutants),function(i){
    
    cell = ids %>% filter(sample == mutants$sample[i],mutant == mutants$mutant[i]) %>% pull(cell_id)
    
    muts %>% filter(cell_id %in% cell) %>% mutate(id = paste0(chr_pos,":",allele,":",ref,":",alt)) %>%
      group_by(id,chr,chr_pos,alt,ref) %>% summarize(n = length(cell_id)) %>% filter(n == mutants$n[i]) %>%
      mutate(mutant = mutants$mutant[i],sample = mutants$sample[i])
    
  }) %>% bind_rows() %>% ungroup()
  
  trunc = tab %>% dplyr::select(id,mutant,sample) %>% unique() %>% group_by(id) %>% 
    summarize(n = length(paste0(sample,"_",mutant))) %>% filter(n == nrow(mutants)) %>% pull(id)
  
  tab = tab %>% filter(!id %in% trunc) %>% rowwise() %>%
    mutate(cell_origin=phylo_forest$get_first_occurrences(Mutation(chr, chr_pos, ref, alt))[[1]]) %>%
    ungroup()
  
  # tab = tab %>% filter(!id %in% trunc) %>% rowwise() %>%
  #   mutate(cell_origin = phylo_forest$get_first_occurrences(SNV(
  #     chr, chr_pos, alt, ref
  #   ))[[1]]) %>%
  #   ungroup()

  
  tab = full_join(tab,phylo_forest$get_nodes() %>% group_by(mutant) %>% summarize(o = min(cell_id)),
                  by = "mutant") %>% filter(cell_origin  > o, o > 0) 
  
  mutants = tab %>% dplyr::select(mutant,sample) %>% unique()
  
  branches = lapply(1:nrow(mutants),function(j){
    
    
    list(ids = tab %>% filter(mutant == mutants$mutant[j],sample == mutants$sample[j]) %>% 
           pull(cell_origin) %>%
           sort() %>% unique(),mutant = mutants$mutant[j],sample = mutants$sample[j])
    
  })
  
  for(j in 1:nrow(mutants)){
    names(branches)[j] = paste0(mutants$mutant,"-",mutants$sample[j])
  }
  
  
  return(branches)
  
}