tool = "mobster"

# source("~/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/generate_table_main.R")
source("/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/generate_table_main.R")

# tool = "mobster"
out_path = get_table_path(save_path, tool, spn, simulation_id)

cli::cli_text("Generating {tool} table for {spn} and simulation {simulation_id}")

path_m = file.path(spn, "tumourevo", simulation_id, "subclonal_deconvolution", tool, "SCOUT", spn)
samples = list.dirs(path_m, recursive=F, full.names=F)  # get_sample_names(spn, base_path=main_path)

print(samples)

if (!dir.exists(path_m)) cli::cli_abort("Mobster directory does not exist")

final_table = lapply(samples, function(sample_name) {
  
  print(file.path(path_m,
                  paste0(sample_name),
                  paste0("SCOUT_", spn, "_", sample_name, "_mobsterh_st_best_fit.rds")))
  
  obj = readRDS(file.path(path_m,
                          paste0(sample_name),
                          paste0("SCOUT_", spn, "_", sample_name, "_mobsterh_st_best_fit.rds")))
  
  final_table_s = obj$data %>% 
    mutate(chr=sub("^chr", "", chr)) %>% 
    mutate(mutation_id=paste0(spn, ":", chr, ":", from,  ":", alt),
           patient_id=spn,
           coverage=coverage,
           purity=purity,
           tool=tool,
           ccf_tool=VAF) %>% 
    select(patient_id, sample_id, coverage, purity, tool, mutation_id, driver_label, is_driver, cluster, VAF, ccf_tool) %>% 
    
    rename(vaf_tool=VAF, cluster_id_tool=cluster,
           driver_label_tool=driver_label,
           is_driver_tool=is_driver) %>%
    add_count(cluster_id_tool, name="n_mutations_tool") %>% 
    mutate(is_tail_tool=cluster_id_tool=="Tail",
           n_clones_tool=n_distinct(cluster_id_tool)) %>% 

    group_by(cluster_id_tool) %>%
    mutate(ccf_tool=mean(ccf_tool, na.rm=TRUE)) %>%
    ungroup() %>% 
    
    mutate(ccf_tool=ccf_tool*2 / purity)
  
  final_table_s %>% 
    mutate(is_clonal_tool=ifelse(cluster_id_tool==get_clonal_cluster_tool(final_table_s), TRUE, FALSE),
           is_subclonal_tool=!is_tail_tool & !is_clonal_tool)
}) %>% bind_rows()



cli::cli_text("Saving {tool} table for {spn} and simulation {simulation_id} in {out_path}")

if (dir.exists(path_m)) saveRDS(final_table, out_path)
