get_clonal_cluster_tool = function(final_table) {
  theta_long = final_table %>%
    group_by(sample_id, cluster_id_tool) %>%
    summarize(mean_ccf=mean(ccf_tool, na.rm=TRUE), .groups="drop")
  
  theta = theta_long %>%
    pivot_wider(
      names_from=cluster_id_tool,
      values_from=mean_ccf
    ) %>%
    select(-sample_id)
  
  max_colnames = apply(theta, 1, function(row) {
    names(row)[which(row == max(row))]
  }) # Extract all clusters which have max ccf for each sample (because in one sample there can be more than one cluster with ccf == 1)
  
  names(which.max(table(unlist(max_colnames)))) # extract the cluster which appear more frequently (i.e. possibly in all the samples)
}


get_clonal_cluster_process = function(final_table) {
  theta_long = final_table %>%
    group_by(sample_id, cluster_id_process) %>%
    summarize(mean_ccf=mean(ccf_process, na.rm=TRUE), .groups="drop")
  
  theta = theta_long %>%
    pivot_wider(
      names_from=cluster_id_process,
      values_from=mean_ccf
    ) %>%
    select(-sample_id)
  
  max_colnames = apply(theta, 1, function(row) {
    names(row)[which(row == max(row))]
  }) # Extract all clusters which have max ccf for each sample (because in one sample there can be more than one cluster with ccf == 1)
  
  names(which.max(table(unlist(max_colnames)))) # extract the cluster which appear more frequently (i.e. possibly in all the samples)
}


get_cell_id = function(mutation_object) {
  tryCatch(
    expr = { phylo_forest$get_first_occurrences(mutation_object)[[1]] },
    error = function(e) return(NA)
  )
}



get_table = function(path, tool, env) {
  filename = get_table_path(path, tool, env$spn, env$simulation_id)
  
  if (file.exists(filename)) {
    cli::cli_text("Loading existing table {file.path(path, filename)}")
    if (nrow(readRDS(filename)) == 0) {
      unlink(filename)
      return(NULL)
    }
    readRDS(filename)
  } else {
    cli::cli_text("Generating table {file.path(path, filename)}")
    tryCatch(expr={
      sys.source(file.path(path, paste0("generate_table_", tool, ".R")), envir=env)
      readRDS(filename)
    },
    error=function(e) {
      cli::cli_alert_warning("Not able to generate the table due to the following error\n{e}")
      return(NULL)
    })
  }
}


get_table_path = function(path, tool, spn, simulation_id) {
  file.path(path, "tables", paste0("table_", tool, "_", spn, "_", simulation_id, ".rds"))
}


select_mutation_specific_columns = function(table_joined) {
  table_joined %>% 
    select(all_of(get_id_colnames()), mutation_id, contains("driver"),
           contains("cluster_id"), contains("is_tail"), contains("is_clonal"),
           contains("is_subclonal"), "causes", contains("vaf")) %>% 
    unique()
}


select_cluster_specific_columns = function(table_joined) {
  table_joined %>% 
    select(all_of(get_id_colnames()), contains("cluster_id"),
           contains("n_clones"), contains("n_mutations"),
           contains("ccf")) %>% 
    unique()
}


get_id_colnames = function() {
  return(
    c("patient_id", "sample_id", "coverage", "purity", "tool")
  )
}


