library(ggplot2)
library(tidyverse)
library(ProCESS)

tool = "process_univariate_w_private"

github_path = "~/GitHub/ProCESS-examples/"
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(main_path, "validation_subclonal/")
save_path_shared = file.path(main_path, "validation_subclonal/")

setwd(main_path)
source(file.path(github_path, "validation/Subclonal_deconvolution/utils_tables.R"))
source(file.path(github_path, "validation/Subclonal_deconvolution/generate_table_main.R"))
source(file.path(github_path, "validation/Subclonal_deconvolution/get_mrca_branches.R"))


coverage_list = c(50, 100, 150)
purity_list = c(0.3, 0.6, 0.9)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = paste("SPN", 3:7, sep="0")
# spn_list = paste("SPN", 1:7, sep="0")

combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller=cna_caller_list)


for (spn in spn_list) {
  
  print(spn)

  subforest_path = file.path("/orfeo/cephfs/scratch/cdslab/shared/SCOUT", spn, "process")
  cli::cli_text("Subforest path: {subforest_path}")
  
  sample_forest = load_sample_forest(get_sample_forest(spn))
  sample_names = (sample_forest$get_samples_info())$name %>% sort()
  phylo_forest = load_phylogenetic_forest(get_phylo_forest(spn=spn))
  mrca_branches_complete = get_mrca_branches(phylo_forest)
  
  mutation_rates_table = lapply(sample_names, function(sample) {
    print(sample)
    test_sample_forest = sample_forest$get_subforest_for(sample)
    
    if(spn %in% c('SPN03', 'SPN04')) {
      test_forest = load_phylogenetic_forest(paste0(subforest_path, "/subforest_", spn, "_", stringr::str_remove(sample, ".*_"), '.sff'))
    } else {
      test_forest = phylo_forest$get_subforest_for(sample)
    }
    
    mutation_rates = test_forest$get_species_info()
    
    sticks = test_forest$get_sticks()
    relevant_branches = get_relevant_branches(test_forest)
    
    #######
    
    if (is.null(relevant_branches)) {
      relevant_branches = data.frame(
        "cell_id" = 0,
        "birth_time" = 0,
        "mutant" = "Clone 1",
        "label" = "Truncal"
      )
    }
    
    mrca_branches = keep(mrca_branches_complete, ~ .x$sample == sample)
    
    if (length(mrca_branches) > 0) {
      for (i in 1:length(mrca_branches)) {
        tmp = mrca_branches[[i]]
        ids = tmp$ids
        target_label = paste0(tmp$mutant, "_", tmp$sample)
        print(target_label)
        for (id in ids) {
          if (id %in% relevant_branches$cell_id) { # if it is already present in relevant_branches, replace the label with the new one
            print(paste0(id, " already present"))
            relevant_branches = relevant_branches %>%
              mutate(label = if_else(cell_id == id, target_label, label))
          } else {
            print("Not already present")
            relevant_branches = relevant_branches %>%
              add_row(
                cell_id = id,
                ancestor = NA_integer_,
                mutant = NA_character_,
                epistate = NA_character_,
                sample = NA_character_,
                birth_time = NA_real_,
                label = target_label
              )
          }
        }
      }
    }
    
    print(mutation_rates)
    
    label_names = relevant_branches %>% select(mutant, label) %>% unique() %>% filter(mutant != label)
    mutation_rates %>% inner_join(label_names) %>% as_tibble() %>%
      select(-epistate, -mutant) %>%
      rename(cluster_id_process=label) %>% 
      mutate(patient_id=spn, sample_id=sample) %>% 
      select(patient_id, sample_id, cluster_id_process, everything())
    # mutate(patient_id=spn, sample_id=sample, coverage=coverage, purity=purity) %>% 
    # select(patient_id, sample_id, coverage, purity, cluster_id_process, everything())
  }) %>% bind_rows() %>% 
    rename_with(.cols=ends_with("rate"), .fn=function(x) paste0(x, "_process")) %>% 
    mutate(sample_id=paste0(spn, "_", sample_id))
  
  for(i in 1:nrow(combs)){
    coverage = combs[i, "coverage"]
    purity = combs[i, "purity"]
    vcf_caller = combs[i, "vcf_caller"]
    cna_caller = combs[i, "cna_caller"]
    
    simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
    
    old_table_path = get_table_path(save_path_shared, tool, spn, simulation_id)
    old_table = readRDS(old_table_path) %>% as_tibble() %>% 
      select(-ends_with("_rate_process"))
    
    new_table = old_table %>% left_join(mutation_rates_table %>% unique())
    
    old_nrow = nrow(old_table)
    new_nrow = nrow(new_table)
    stopifnot(old_nrow == new_nrow)
    
    # Save updated table #####
    cli::cli_text("Save `new_table` in {old_table_path}")
    saveRDS(new_table, old_table_path)
  }
}


