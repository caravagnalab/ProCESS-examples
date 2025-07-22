library(dplyr)
library(tidyr)
library(mobster)
library(ggplot2)
setwd('/u/cdslab/erivar00/scratch/GitHub/ProCESS-examples/')
source("getters/process_getters.R")
source("getters/tumourevo_getters.R")

main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
setwd(main_path)

tools = c("mobster", "pyclonevi", "viber")
getwd()

# Use getters here!!
coverage = 100
purity = 0.9
vcf_caller = 'mutect2'
cna_caller = "sequenza"
spn = "SPN03"

samples <- get_sample_names(SPN)
# Process ####
mut_process = get_mutations(spn = spn, type = 'tumour', coverage = coverage, purity = purity)
mut_process = readRDS(mut_process)

mut_process = mut_process %>% 
  mutate(mutation_id = paste0(spn, ":", chr, ":", chr_pos,  ":", alt))

mut_process$is_driver <- mut_process$classes == "driver"

mut_process <- mut_process %>%
  mutate(causes = replace(causes, classes != 'driver', NA))

View(mut_process)

mut_process_new <- mut_process %>%
  ungroup() %>% 
  select(mutation_id, causes, is_driver, contains(".VAF")) %>%
  pivot_longer(
    cols = ends_with(".VAF"),
    names_to = "sample_id",
    names_pattern = "(.*)\\.VAF", # remove matching text "VAF" from the start of each variable name
    values_to = "vaf_process" # this is the VAF!
  ) %>%
  mutate(sample_id = paste0(spn, "_", sample_id))
# saveRDS(mut_process_new, paste0("/u/cdslab/erivar00/scratch/GitHub/subclonal_validation_data/mut_process_new_", spn, ".Rds"))

spn = 'SPN04'
# Use getters here!!
mut_process_new = readRDS(paste0("/u/cdslab/erivar00/scratch/GitHub/subclonal_validation_data/mut_process_new_", spn, ".Rds"))
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

final_table <- final_table %>%
  group_by(cluster_id_tool, sample_id) %>%
  mutate(ccf_tool = mean(ccf_tool, na.rm = TRUE)) %>%
  ungroup()
View(final_table)

final_table %>%
  distinct(cluster_id_tool, sample_id, ccf_tool)

final_table %>%
  distinct(sample_id, mutation_id) %>%
  nrow()
nrow(final_table)


# Test if all the ccf for all the samples are present
final_table %>%
  distinct(cluster_id_tool, sample_id, ccf_tool) %>% nrow()

#### Extract clonal cluster ####
theta_long <- final_table %>%
  group_by(sample_id, cluster_id_tool) %>%
  summarize(mean_ccf = mean(ccf_tool, na.rm = TRUE), .groups = "drop")
# 
theta <- theta_long %>%
  pivot_wider(
    names_from = cluster_id_tool,
    values_from = mean_ccf
  ) %>%
  select(-sample_id)

max_colnames <- apply(theta, 1, function(row) {
  names(row)[which(row == max(row))]
}) # Extract all clusters which have max ccf for each sample (because in one sample there can be more than one cluster with ccf == 1)

clonal_cluster = names(which.max(table(unlist(max_colnames)))) # extract the cluster which appear more frequently (i.e. possibly in all the samples)

# theta <- as.data.frame(theta)
# clonal_cluster = apply(theta, 1, which.max) # 1 indicates row so it finds the maximum column for each row, i.e. the top cluster for each sample
# clonal_cluster = colnames(theta)[clonal_cluster] # takes the column name of theta for the indices found above, i.e. the top cluster for each sample
# clonal_cluster = which.max(table(clonal_cluster)) %>% names # table(clonal_cluster) creates a table of frequencies of top cluster in the samples

final_table$is.clonal = FALSE

final_table <- final_table %>%
  mutate(is.clonal = ifelse(cluster_id_tool == clonal_cluster, TRUE, FALSE))

### Join pyclonevi and process ####
pyclone_join <- inner_join(mut_process_new, final_table, by = c("mutation_id", "sample_id"))
View(pyclone_join)
nrow(mut_process_new)
nrow(final_table)
nrow(pyclone_join)
pyclone_join %>% count(causes)


### Plot fit pyclonevi ####
df = input %>%
  left_join(best_fit, by = c("sample_id", "mutation_id"))
View(df)


df <- df %>%
  mutate(vaf = alt_counts / (alt_counts + ref_counts))
View(df)
vaf_wide <- df %>%
  select(mutation_id, sample_id, cluster_id, vaf) %>%
  pivot_wider(names_from = sample_id, values_from = vaf)
View(vaf_wide)
samples = unique(df$sample_id)

ggplot(vaf_wide, aes(x = .data[[samples[[1]]]], y = .data[[samples[[4]]]], color = as.factor(cluster_id))) +
  geom_point(alpha = 0.8, size = 2) +
  labs(x = samples[[1]], y = samples[[3]], color = "Cluster") +
  theme_minimal()

pairs <- combn(samples, 2, simplify = FALSE)
library(purrr)
plots <- map(pairs, function(pair) {
  ggplot(vaf_wide, aes(x = .data[[pair[1]]], y = .data[[pair[2]]], color = as.factor(cluster_id))) +
    geom_point(alpha = 0.7, size = 1.8) +
    labs(
      x = pair[1],
      y = pair[2],
      title = paste("VAF:", pair[1], "vs", pair[2]),
      color = "Cluster"
    ) +
    theme_minimal()
})

library(patchwork)
wrap_plots(plots, ncol = 2)

# Viber ####
tool = "viber"
path_v = paste0(spn, "/tumourevo/", coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller, "/subclonal_deconvolution/", tool, "/SCOUT/", spn, "/")
viber = readRDS(paste0(path_v, "SCOUT_", spn, "_viber_best_st_fit.rds"))

viber_fit <- bind_cols(viber$data, cluster = viber$labels$cluster.Binomial)
View(viber_fit)
# saveRDS(viber_fit, paste0("/u/cdslab/erivar00/scratch/GitHub/subclonal_validation_data/viber_fit.Rds"))

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

viber_fit_long <- viber_fit_long %>%
  group_by(cluster_id_tool, sample_id) %>%
  mutate(ccf_tool = mean(ccf_tool, na.rm = TRUE)) %>%
  ungroup()

View(viber_fit_long)

viber_fit_long$ccf_tool = (viber_fit_long$ccf_tool*(purity+2))/2

viber_fit_long %>%
  distinct(sample_id, mutation_id) %>%
  nrow()
nrow(viber_fit_long)

#### Extract clonal cluster ####

theta = viber$theta_k
clonal_cluster = apply(theta, 1, which.max) # 1 indicates row so it finds the maximum column for each row, i.e. the top cluster for each sample
clonal_cluster = colnames(theta)[clonal_cluster] # takes the column name of theta for the indices found above, i.e. the top cluster for each sample
clonal_cluster = which.max(table(clonal_cluster)) %>% names # table(clonal_cluster) creates a table of frequencies of top cluster in the samples

viber_fit_long = viber_fit_long %>%
  mutate(is.clonal = ifelse(cluster_id_tool == clonal_cluster, TRUE, FALSE))


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
  sample_name = sample_names[1]
  mobster = readRDS(paste0(path_m, spn,  "_", spn, "_", sample_name, "/",  "SCOUT_",  spn, "_", spn, "_", spn, "_", sample_name, "_mobsterh_st_best_fit.rds"))
  
  #
  # x = mobster
  # 
  # # mobster:::is_mobster_fit(x)
  # 
  # x$Clusters
  # cluster_table = x$Clusters %>%
  #   dplyr::filter(cluster != 'Tail', type == 'Mean') %>%
  #   dplyr::select(cluster, fit.value) %>%
  #   rename(ccf = fit.value)
  # cluster_table$nMuts = x$N.k[cluster_table$cluster]
  # 
  # # Clonality status - maximum fit is the clonal
  # cluster_table$is.clonal = FALSE
  # cluster_table$is.clonal[which.max(cluster_table$ccf)] = TRUE
  
  
  # all_fit = readRDS(paste0(path_m, spn,  "_", spn, "_", sample_name, "/",  "SCOUT_",  spn, "_", spn, "_", spn, "_", sample_name, "_mobsterh_st_fit.rds"))
  # 
  # evolutionary_parameters(all_fit)
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
  
  final_mobster <- final_mobster %>%
    group_by(cluster_id_tool) %>%
    mutate(ccf_tool = mean(ccf_tool, na.rm = TRUE)) %>%
    ungroup()
  
  View(final_mobster)
  
  final_mobster$ccf_tool = (final_mobster$ccf_tool*(purity+2))/2
  
  mobster_results[[sample_name]] <- final_mobster
}

final_table <- bind_rows(mobster_results)
View(final_table)

final_table %>%
  distinct(sample_id, mutation_id) %>%
  nrow()
nrow(final_table)

### Join mobster and process ####
mobster_join <- inner_join(mut_process_new, final_table, by = c("mutation_id", "sample_id"))
View(mobster_join)
nrow(mut_process_new)
nrow(final_table)
nrow(mobster_join)
