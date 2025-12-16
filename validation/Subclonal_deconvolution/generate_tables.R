library(ggplot2)
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
cat(paste("\nArguments:", paste(args, collapse=", "), "\n"))

join_tables = as.logical(args[1])
i = as.integer(args[2])
print(i)

# github_path = "~/GitHub/ProCESS-examples/"
github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")
# save_path_shared = file.path(main_path, "subclonal/")

setwd(main_path)
source(file.path(save_path, "utils_plots.R"))

# coverage_list = c(50, 100, 150, 200)
# purity_list = c(0.3, 0.6, 0.9)
# vcf_caller_list = c("mutect2", "strelka", "freebayes")
# cna_caller_list = c("ascat", "sequenza", "battenberg")
# spn_list = paste("SPN", 3:7, sep="0")

coverage_list = c(50,100, 150)
purity_list = c(0.3, 0.6, 0.9)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = paste("SPN", 1:7, sep="0")

combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller=cna_caller_list,
                    spn=spn_list)

# if (!is.na(i)) {
for(i in 1:nrow(combs)){
  # missing_files = readRDS(file.path(save_path, "missing_files.rds"))
  # if (!i %in% missing_files$i) cli::cli_abort("File already present")
  
  coverage = combs[i, "coverage"]
  purity = combs[i, "purity"]
  vcf_caller = combs[i, "vcf_caller"]
  cna_caller = combs[i, "cna_caller"]
  spn = combs[i, "spn"]
  
  simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
  
  env = list2env(list(spn=spn, coverage=coverage, purity=purity,
                      vcf_caller=vcf_caller, cna_caller=cna_caller,
                      simulation_id=simulation_id),
                 parent=globalenv())
  
  # table_pyclonevi = get_table(save_path, "pyclonevi", env)
  
  table_mobster = get_table(save_path, "mobster", env)
  # 
  # table_viber = get_table(save_path, "viber", env)
  # 
  # table_viber_h = get_table(save_path, "viber_heuristics", env)
  # 
  # table_process = get_table(save_path, "process", env)
}


if (join_tables) {
  full_table = data.frame()
  for(i in 1:nrow(combs)) {
    coverage = combs[i, "coverage"]
    purity = combs[i, "purity"]
    vcf_caller = combs[i, "vcf_caller"]
    cna_caller = combs[i, "cna_caller"]
    spn = combs[i, "spn"]
    simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)

    env = list2env(list(spn=spn, coverage=coverage, purity=purity,
                        vcf_caller=vcf_caller, cna_caller=cna_caller,
                        simulation_id=simulation_id),
                   parent=globalenv())

    table_pyclonevi = get_table(save_path, "pyclonevi", env)
    table_mobster = get_table(save_path, "mobster", env)
    table_viber = get_table(save_path, "viber", env)
    table_process = get_table(save_path, "process", env)

    if (all(sapply(list(table_pyclonevi, table_mobster, table_viber), is.null)))
      next
    
    if (is.null(table_process))
      next

    table_joined_tools = rbind(table_pyclonevi, table_mobster, table_viber)

    table_joined = table_joined_tools %>%
      full_join(table_process) %>%
      mutate(cluster_id_process=ifelse(is.na(cluster_id_process), "FP", cluster_id_process),
             tool=ifelse(is.na(tool), "FN", tool)) %>%
      filter(tool!="FN") %>%
      mutate(cluster_id_tool=as.character(cluster_id_tool))

    full_table = bind_rows(full_table, table_joined) %>% 
      mutate(vcf_caller=vcf_caller, cna_caller=cna_caller)
    rm(table_joined)
  }

  saveRDS(full_table, file.path(save_path, "full_table_stats.rds"))
}

