tool = "viber"

source("~/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/generate_table_main.R")

# tool = "viber"
out_path = get_table_path(save_path, tool, spn, simulation_id)

cli::cli_text("Generating {tool} table for {spn} and simulation {simulation_id}")

path_v = get_tumourevo_subclonal(spn,
                                 coverage=coverage,
                                 purity=purity, 
                                 vcf_caller=vcf_caller,
                                 cna_caller=cna_caller,
                                 sample=spn,
                                 tool=tool)

obj = readRDS(path_v$viber_best_st_fit_rds)

final_table = obj$data %>% 
  mutate(cluster_id_tool=obj$labels$cluster.Binomial) %>% 
  mutate(chr=sub("^chr", "", chr),
         mutation_id=paste(spn, chr, from, alt, sep=":"),
         patient_id=spn,
         coverage=coverage,
         purity=purity,
         tool=tool) %>% 
  select(patient_id, coverage, purity, tool, mutation_id, cluster_id_tool, starts_with("VAF")) %>% 
  
  mutate(driver_label_tool=NA, is_driver_tool=FALSE) %>%
  
  pivot_longer(cols=starts_with("VAF"), names_to="sample_id", values_to="ccf_tool", names_prefix="VAF.") %>%
  
  add_count(cluster_id_tool, name="n_mutations_tool") %>%
  mutate(vaf_tool=ccf_tool,
         is_tail_tool=cluster_id_tool=="Tail",
         n_clones_tool = n_distinct(cluster_id_tool)) %>% 
  
  group_by(cluster_id_tool, sample_id) %>%
  mutate(ccf_tool=mean(ccf_tool, na.rm=TRUE)) %>%
  ungroup() %>% 
  
  mutate(ccf_tool=ccf_tool*2 / purity)

final_table = final_table %>% 
  mutate(is_clonal_tool=ifelse(cluster_id_tool==get_clonal_cluster_tool(final_table), TRUE, FALSE),
         is_subclonal_tool=!is_tail_tool & !is_clonal_tool)



cli::cli_text("Saving {tool} table for {spn} and simulation {simulation_id} in {out_path}")

saveRDS(final_table, out_path)
