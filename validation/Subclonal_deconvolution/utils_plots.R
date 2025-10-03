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
  file.path(save_path, "tables", paste0("table_", tool, "_", spn, "_", simulation_id, ".rds"))
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

