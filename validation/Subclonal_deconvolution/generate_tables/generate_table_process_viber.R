tool = "process_viber"

# source("~/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/generate_table_main.R")
source("/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/generate_table_main.R")

# tool = "viber"
out_path = get_table_path(save_path, tool, spn, simulation_id)

cli::cli_text("Generating {tool} table for {spn} and simulation {simulation_id}")

path_v = paste0("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/test_subclonal/v4_process/", 
                spn,"_", coverage, "x_", purity, "_viber_best_fit.rds")

obj = readRDS(path_v)

final_table = obj$data %>% 
  mutate(cluster_id_tool=obj$labels$cluster.Binomial) %>% 
  mutate(chr=sub("^chr", "", chr),
         mutation_id=paste(spn, chr, from, alt, sep=":"),
         patient_id=spn,
         coverage=coverage,
         purity=purity,
         tool=tool) %>% 
  mutate(driver_label_tool=gene, is_driver_tool=driver) %>%
  select(patient_id, coverage, purity, tool, mutation_id, cluster_id_tool, starts_with("VAF"), driver_label_tool, is_driver_tool) %>% 
  
  # mutate(driver_label_tool=NA, is_driver_tool=FALSE) %>%
  
  pivot_longer(cols=starts_with("VAF"), names_to="sample_id", values_to="ccf_tool", names_prefix="VAF.") %>%
  
  add_count(cluster_id_tool, name="n_mutations_tool") %>%
  mutate(vaf_tool=ccf_tool,
         is_tail_tool=cluster_id_tool=="Tail",
         n_clones_tool = n_distinct(cluster_id_tool)) %>% 
  
  # group_by(cluster_id_tool, sample_id) %>%
  group_by(cluster_id_tool) %>%
  mutate(ccf_tool=mean(ccf_tool, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(ccf_tool=ccf_tool*2 / purity)

cli::cli_text("Saving {tool} table for {spn} and simulation {simulation_id} in {out_path}")

saveRDS(final_table, out_path)
