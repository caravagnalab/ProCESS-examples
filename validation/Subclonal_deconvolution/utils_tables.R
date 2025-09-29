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
