library(dplyr)
library(tidyr)
setwd('/u/cdslab/erivar00/scratch/GitHub/ProCESS-examples/')
source("getters/process_getters.R")
source("getters/tumourevo_getters.R")

main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
setwd(main_path)

tools = c("mobster", "pyclonevi", "viber")
getwd()

coverage = 50
purity = 0.6
vcf_caller = 'mutect2'
cna_caller = "ascat"
spn = "SPN04"

# Process ####
mut_process = get_mutations(spn = spn, type = 'tumour', coverage = coverage, purity = purity)
mut_process = readRDS(mut_process)

mut_process = mut_process %>% 
  mutate(mutation_id = paste0(spn, ":", chr, ":", chr_pos,  ":", alt))

mut_process$is_driver <- mut_process$classes == "driver"

mut_process <- mut_process %>%
  mutate(causes = if_else(classes == "driver", causes, NA))

View(mut_process)

mut_process_new <- mut_process %>%
  ungroup() %>% 
  select(mutation_id, causes, is_driver, contains(".VAF")) %>%
  pivot_longer(
    cols = ends_with(".VAF"),
    names_to = "sample_id",
    names_pattern = "(.*)\\.VAF", # remove matching text "VAF" from the start of each variable name
    values_to = "ccf_process" # this is the VAF!
  ) %>%
  mutate(sample_id = paste0(spn, "_", sample_id))
View(mut_process_new)

# library(ProCESS)
# forest_path = get_phylo_forest(spn = spn)
# phylo_forest = load_phylogenetic_forest(forest_path)
# 
# # forest_cna = phylo_forest$get_sampled_cell_CNAs()
# # saveRDS(forest_cna, paste0("/u/cdslab/erivar00/scratch/GitHub/subclonal_validation_data/forest_cna_", spn, ".Rds"))
# forest_cna = readRDS(paste0("/u/cdslab/erivar00/scratch/GitHub/subclonal_validation_data/forest_cna_", spn, ".Rds"))
# View(forest_cna)
# 
# # forest_cells = phylo_forest$get_sampled_cell_mutations()
# forest_cells = readRDS(paste0("/u/cdslab/erivar00/scratch/GitHub/subclonal_validation_data/forest_cells_mut_", spn, ".Rds"))
# View(forest_cells)
# # saveRDS(forest_cells, paste0("/u/cdslab/erivar00/scratch/GitHub/subclonal_validation_data/forest_cells_mut_", spn, ".Rds"))


# PyClone ####
tool = "pyclonevi"
path_p = paste0(spn, "/tumourevo/", coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller, "/subclonal_deconvolution/", tool, "/SCOUT/", spn, "/")
input = read.delim(paste0(path_p,"SCOUT_", spn, "_pyclone_input.tsv"), sep="\t")
# View(input) # for each mutation_id there is a row for each sample in which it is present (i.e. 4 because even if vaf = 0 the row is present)

# cluster = read.csv(paste0(path_p,"SCOUT_", spn, "_cluster_table.csv"), sep="\t")
# View(cluster) # cluster per mutation_id, so there is a single row for each mutation (even if it is present in more samples)

best_fit = read.delim(paste0(path_p,"SCOUT_", spn, "_best_fit.txt"))
# View(best_fit) # like the input, with also cluster_id.

# I think I need to take input and add cluster_id and cellular_prevalence which are inside pyclone
# Columns that I need from input: 
# patient_id, sample_id, mutation_id, driver_label, is_driver

# Columns that I need from pyclone: cluster_id, cellular_prevalence
# What I need to add: tool = pyclone, purity = 0.6, coverage = 50

final_table = input %>%
  left_join(best_fit, by = c("sample_id", "mutation_id")) %>%
  # select(patient_id, sample_id, mutation_id, driver_label, is_driver, cluster_id, cellular_prevalence) %>%
  select(patient_id, sample_id, mutation_id, driver_label, cluster_id, cellular_prevalence) %>%
  rename(ccf_tool = cellular_prevalence, cluster_id_tool = cluster_id) %>%
  add_count(cluster_id_tool, name = "n_mutations_tool") %>% 
  filter(!is.na(cluster_id_tool)) %>% # not sure
  mutate(
    purity = purity,
    coverage = coverage,
    tool = tool,
    n_clones_tool = n_distinct(cluster_id_tool)
  )
View(final_table)

final_table %>%
  distinct(sample_id, mutation_id) %>%
  nrow()
nrow(final_table)

### Join pyclonevi and process ####
pyclone_join <- inner_join(mut_process_new, final_table, by = c("mutation_id", "sample_id"))
View(pyclone_join)
nrow(mut_process_new)
nrow(final_table)
nrow(pyclone_join)
pyclone_join %>% count(causes)

# Viber ####

tool = "viber"
path_v = paste0(spn, "/tumourevo/", coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller, "/subclonal_deconvolution/", tool, "/SCOUT/", spn, "/")
viber = readRDS(paste0(path_v, "SCOUT_", spn, "_viber_best_st_fit.rds"))

viber_fit <- bind_cols(viber$data, cluster = viber$labels$cluster.Binomial)
View(viber_fit)

viber_fit = viber_fit %>% 
  mutate(chr = sub("^chr", "", chr)) %>%
  mutate(mutation_id = paste0(spn, ":", chr, ":", from,  ":", alt))

viber_fit_long <- viber_fit %>%
  select(mutation_id, cluster, matches("^VAF\\.")) %>%
  pivot_longer(
    cols = starts_with("VAF."),
    names_to = "sample_id",
    names_prefix = "VAF.", # remove matching text "VAF" from the start of each variable name
    values_to = "ccf_tool" # this is the VAF!
  ) %>%
  rename(cluster_id_tool = cluster)  %>%
  add_count(cluster_id_tool, name = "n_mutations_tool") %>%
  mutate(
    purity = purity,
    coverage = coverage,
    patient_id = spn,
    tool = tool,
    n_clones_tool = n_distinct(cluster_id_tool)
  )

viber_fit_long %>%
  distinct(sample_id, mutation_id) %>%
  nrow()
nrow(viber_fit_long)

### Join viber and process ####
viber_join <- inner_join(mut_process_new, viber_fit_long, by = c("mutation_id", "sample_id"))
View(viber_join)
nrow(mut_process_new)
nrow(viber_fit_long)
nrow(viber_join)

# Mobster ####
tool = "mobster"
path_m = paste0(spn, "/tumourevo/", coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller, "/subclonal_deconvolution/", tool, "/SCOUT/", spn, "/")

subdirs <- list.dirs(path_m, full.names = FALSE, recursive = FALSE)
sample_names <- sub("^.*?_.*?_(.*)$", "\\1", subdirs)

mobster_results <- list()
for (sample_name in sample_names) {
  mobster = readRDS(paste0(path_m, spn,  "_", spn, "_", sample_name, "/",  "SCOUT_",  spn, "_", spn, "_", spn, "_", sample_name, "_mobsterh_st_best_fit.rds"))
  
  mobster_fit = mobster$data
  
  mobster_fit = mobster_fit %>% 
    mutate(chr = sub("^chr", "", chr)) %>%
    mutate(mutation_id = paste0(spn, ":", chr, ":", from,  ":", alt))
  
  final_mobster = mobster_fit %>%
    select(sample_id, mutation_id, driver_label, is_driver, cluster, VAF) %>%
    rename(ccf_tool = VAF, cluster_id_tool = cluster) %>%
    add_count(cluster_id_tool, name = "n_mutations_tool") %>% 
    mutate(
      purity = purity,
      coverage = coverage,
      patient_id = spn,
      tool = tool,
      n_clones_tool = n_distinct(cluster_id_tool)
    )
  mobster_results[[sample_name]] <- final_mobster
}

final_table <- bind_rows(mobster_results)
View(final_table)

final_table %>%
  distinct(sample_id, mutation_id) %>%
  nrow()
nrow(final_table)

### Join viber and process ####
mobster_join <- inner_join(mut_process_new, final_table, by = c("mutation_id", "sample_id"))
View(mobster_join)
nrow(mut_process_new)
nrow(final_table)
nrow(mobster_join)

