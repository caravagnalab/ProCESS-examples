tool = "pyclonevi"

# source("~/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/generate_table_main.R")
source("/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/generate_table_main.R")

# tool = "pyclonevi"
out_path = get_table_path(save_path, tool, spn, simulation_id)

cli::cli_text("Generating {tool} table for {spn} and simulation {simulation_id}")

path_p = get_tumourevo_subclonal(spn,
                                 coverage=coverage,
                                 purity=purity, 
                                 vcf_caller=vcf_caller,
                                 cna_caller=cna_caller,
                                 sample=spn,
                                 tool=tool)

input = read.delim(path_p$pyclone_input_tsv, sep="\t", row.names=1) %>% as_tibble()
# for each mutation_id there is a row for each sample in which it is present (i.e. 4 because even if vaf = 0 the row is present)

best_fit = read.delim(path_p$best_fit_txt) %>% as_tibble() # like the input, with also cluster_id

# I think I need to take input and add cluster_id and cellular_prevalence which are inside pyclone
# Columns that I need from input: 
# patient_id, sample_id, mutation_id, driver_label, is_driver

# Columns that I need from pyclone: cluster_id, cellular_prevalence
# What I need to add: tool = pyclone, purity = 0.6, coverage = 50

final_table = input %>%
  left_join(best_fit, by=c("sample_id", "mutation_id")) %>%
  
  mutate(vaf_tool=alt_counts / (ref_counts + alt_counts)) %>% 
  
  select(patient_id, sample_id, mutation_id, driver_label, is_driver, cluster_id, cellular_prevalence, vaf_tool) %>%
  
  rename(ccf_tool=cellular_prevalence, cluster_id_tool=cluster_id, 
         is_driver_tool=is_driver, driver_label_tool=driver_label) %>%
  mutate(is_driver_tool=ifelse(is_driver_tool=="False", FALSE, TRUE),
         driver_label_tool=replace(driver_label_tool, !is_driver_tool, NA)) %>% 
  group_by(cluster_id_tool) %>% 
  mutate(n_mutations_tool=length(unique(mutation_id))) %>% ungroup() %>% 
  filter(!is.na(cluster_id_tool)) %>% # not sure
  mutate(purity=purity,
         coverage=coverage,
         tool=tool,
         n_clones_tool=n_distinct(cluster_id_tool),
         is_tail_tool=FALSE) %>% 
  
  group_by(cluster_id_tool, sample_id) %>%
  mutate(ccf_tool=mean(ccf_tool, na.rm = TRUE)) %>%
  ungroup()

# add clonal cluster
final_table = final_table %>% 
  mutate(is_clonal_tool=ifelse(cluster_id_tool==get_clonal_cluster_tool(final_table), TRUE, FALSE),
         is_subclonal_tool=!is_tail_tool & !is_clonal_tool) %>% 
  
  select(patient_id, sample_id, coverage, purity, tool, mutation_id, dplyr::everything())



cli::cli_text("Saving {tool} table for {spn} and simulation {simulation_id} in {out_path}")

saveRDS(final_table, out_path)
