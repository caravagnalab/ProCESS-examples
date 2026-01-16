library(ggplot2)
library(tidyverse)

# table_path = "/orfeo/cephfs/scratch/cdslab/ebusca00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/full_table_stats.rds"

spn = 'SPN04'
coverage=100
purity=0.9
# vcf_caller = "mutect2"
# cna_caller = "ascat"
simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

tool = 'process'
table_process = readRDS(get_table_path(save_path, tool, spn, simulation_id))

tool = 'viber'
table_viber = readRDS(get_table_path(save_path, tool, spn, simulation_id))

### Plots process ####
plot_scatter = function(table_wide, s1,s2){
  ggplot(table_wide, aes(x = eval(parse(text = s1)),
                                y = eval(parse(text = s2)),
                                color = cluster_id_process)) +
    geom_point() +
    theme_minimal() +
    labs(
      x = s1,
      y = s2
    )

}
# Read table
sample_names <- final_table$sample_id %>% unique()

table_wide = final_table %>%
  select(patient_id, sample_id, mutation_id, cluster_id_process, vaf_process) %>%
  pivot_wider(values_from="vaf_process", names_from="sample_id")

table_wide <- table_wide %>%
  filter(!is.na(cluster_id_process))
table_wide[is.na(table_wide)] = 0.0

sample_names = as.character(sample_names)
cm = combn(sample_names, 2)

plots <- apply(
  cm,
  2,
  function(w) plot_scatter(table_wide, s1 = w[1], s2 = w[2]) # w is the sample name
)

library(ggpubr)

ggarrange(
  plotlist = plots,
  ncol = 3,
  nrow = 2)


### Map mutations on tree #### 

sample_forest = load_sample_forest(get_sample_forest(spn))
phylo_forest = load_phylogenetic_forest(get_phylo_forest(spn=spn))

viber_table = readRDS("/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/tables/table_viber_SPN04_100x_0.9p_mutect2_ascat.rds")
# process_table = readRDS("/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/tables/table_process_SPN04_100x_0.9p_mutect2_ascat.rds")

mut_process = get_mutations(spn=spn, type="tumour", coverage=coverage, purity=purity)
seq_results = readRDS(mut_process)

process_table = seq_results %>%
  mutate(mutation_id=paste0(spn, ":", chr, ":", chr_pos,  ":", alt),
         is_driver_process=classes=="driver")

table_wide = viber_table %>%
  select(patient_id, sample_id, mutation_id, cluster_id_tool, vaf_tool) %>%
  pivot_wider(values_from="vaf_tool", names_from="sample_id")

table_wide <- table_wide %>%
  filter(!is.na(cluster_id_tool))

table_wide[is.na(table_wide)] = 0.0

s1 = "SPN04_SPN04_1.1"
s2 = "SPN04_SPN04_2.1"
ggplot(table_wide, aes(x = eval(parse(text = s1)),
                       y = eval(parse(text = s2)),
                       color = cluster_id_tool)) +
  geom_point() +
  theme_minimal() +
  labs(
    x = s1,
    y = s2
  ) +
  scale_color_manual(values=color_palette)


# Here I need to do a join between the process table of this simulation and the viber table
join_table = viber_table %>%
  inner_join(process_table) %>% 
  select(chr, chr_pos, ref, alt, mutation_id, cluster_id_tool, vaf_tool)

mutations_with_cell = join_table %>% 
  rowwise() %>%
  mutate(cell_id=phylo_forest$get_first_occurrences(Mutation(
    chr, chr_pos, ref, alt
  ))[[1]]) %>%
  ungroup()

cells_labels = mutations_with_cell %>% 
  select(mutation_id, cell_id, cluster_id_tool) %>% 
  group_by(cell_id) %>% 
  summarise(label_list=list(cluster_id_tool)) %>% 
  rowwise() %>% 
  mutate(label=names(sort(table(label_list[[1]]), decreasing=TRUE))[1]) %>% 
  ungroup() %>% 
  select(-label_list)

final_labels = sample_forest$get_nodes() %>% as_tibble() %>% 
  left_join(cells_labels)

color_palette = c(
  "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00","#a65628",
  "#FFD700",  "#999999", "#000000", "#f781bf", # First 10 colors (Set1)  
  "#46f0f0", "#f032e6", "#bcf60c", "#fabed4", "#008080", "#e6beff",  
  "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",  
  "#000075", "#808080", "#d3a6f3", "#ff9cdd", "#73d7b0"  
) %>% setNames(str_sort(unique(viber_table$cluster_id_tool), numeric=T))

pl_sticks = plot_sticks(sample_forest, labels=final_labels, cls = color_palette) %>%
  annotate_forest(sample_forest, samples=TRUE, drivers=TRUE)
pl_sticks
