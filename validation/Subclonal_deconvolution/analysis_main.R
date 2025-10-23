library(ggplot2)
library(tidyverse)
library(ProCESS)
library(ggpubr)
library(patchwork)
library(randnet)
library(scales)

spn = 'SPN03'
coverage=100
purity=0.9
vcf_caller = "mutect2"
cna_caller = "ascat"
simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

source(file.path(save_path, "utils_plots_final.R"))
source(file.path(save_path, "generate_table_main.R"))

drivers_table = readRDS(file.path(main_path,"/drivers/all_drivers.rds"))

drivers_table = drivers_table %>% filter(SPN==spn)  %>%
  mutate(mutation_id=paste0(spn, ":", chr, ":", start,  ":", alt)) %>% 
  select(mutation_id, code)

# Get process table
mut_process = get_mutations(spn=spn, type="tumour", coverage=coverage, purity=purity)
table_process = readRDS(get_table_path(save_path, 'process', spn, simulation_id)) # process table in folder tables/

# Join process table with drivers
table_process = table_process %>% left_join(drivers_table, by = 'mutation_id') %>% 
  select(-driver_label_process) %>% 
  mutate(driver_label_process=code)

# SPN03:9:136496197: deve diventare SPN03:9:136496196:C
if(spn=='SPN03'){
  table_process = table_process %>% mutate(mutation_id = if_else(
    mutation_id=="SPN03:9:136496197:",
    'SPN03:9:136496196:C',
    mutation_id
  ))
}

tool = 'viber'
table_tool = readRDS(get_table_path(save_path, tool, spn, simulation_id))

if(tool =='viber_heuristics'){
  table_tool = table_tool %>%
    mutate(cluster_id_tool = replace_na(cluster_id_tool, 'NA'))
}

join_table_tool = table_tool %>% left_join(table_process) # keep all mut in tool
join_table_process = join_table_tool %>% filter(!is.na(cluster_id_process)) # only mutations present in both

nmi_complete = randnet::NMI(as.factor(join_table_process$cluster_id_tool),
                    as.factor(join_table_process$cluster_id_process))

ari_complete = aricode::ARI(as.factor(join_table_process$cluster_id_tool), 
            as.factor(join_table_process$cluster_id_process))

### Find cluster/driver in tool
driver_clusters_tool = join_table_tool %>%
  distinct(cluster_id_tool, is_driver_tool) %>% 
  filter(is_driver_tool == TRUE) %>% 
  pull(cluster_id_tool)

final_table = join_table_tool %>%
  mutate(cluster_id_tool_interpreted = if_else(
    !(cluster_id_tool %in% driver_clusters_tool),
    'Subclonal',
    as.character(cluster_id_tool)
  ))


### Plot tool with new labels ####
color_palette_tool = c(
  "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00","#a65628",
  "#FFD700", "#000000", "#f781bf", # First 10 colors (Set1)  
  "#46f0f0", "#f032e6", "#bcf60c", "#fabed4", "#008080", "#e6beff",  
  "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",  
  "#000075", "#808080", "#d3a6f3", "#ff9cdd", "#73d7b0"  ) %>% setNames(str_sort(unique(final_table$cluster_id_tool_interpreted), numeric=T))

color_palette_tool["Subclonal"] = "#cccccc"

sample_names = sort(unique(table_process$sample_id))

# Scatterplot tool
plot_tool = plot_scatter_tool(final_table, color_palette_tool, sample_names, type ='interpreted')
plot_tool

# Scatterplot process
color_palette_process = hue_pal()(length(unique(table_process$cluster_id_process))) %>% 
  setNames(str_sort(unique(table_process$cluster_id_process), numeric=T))
color_palette_process["Subclonal"] = "#cccccc"

scatter_process = plot_scatter_process(join_table_process, sample_names, color_palette_process, driver=T)
scatter_process

final_table_interpreted = join_table_process %>%
  mutate(cluster_id_tool_interpreted = if_else(
    !(cluster_id_tool %in% driver_clusters_tool),
    'Subclonal',
    as.character(cluster_id_tool)
  ))
nmi_interpreted = randnet::NMI(as.factor(final_table_interpreted$cluster_id_tool_interpreted),
                             as.factor(final_table_interpreted$cluster_id_process))

ari_interpreted = aricode::ARI(as.factor(final_table_interpreted$cluster_id_tool_interpreted),
                            as.factor(final_table_interpreted$cluster_id_process))

# La tabella da salvare è final_table_interpreted, dove devo salvare le colonne:
# patient_id, sample_id, coverage, purity, tool, mutation_id, driver_label_tool,
# is_driver_tool, cluster_id_tool, vaf_tool, is_driver_process, cluster_id_process,
# vaf_process, driver_label_process, cluster_id_tool_interpreted
table_to_save = final_table_interpreted %>% select(patient_id, sample_id,coverage, purity, 
                                                   tool, mutation_id, driver_label_tool,
                                                   is_driver_tool, cluster_id_tool, vaf_tool, 
                                                   is_driver_process, cluster_id_process,
                                                   vaf_process, driver_label_process, 
                                                   cluster_id_tool_interpreted)

saveRDS(table_to_save, file.path(main_path, "subclonal/metrics_tables", paste0(tool, "_", spn, "_", simulation_id, ".rds")))

if(spn=='SPN03'){
  width=40
  height=50
  design="aaaaaa\nbbbbbb"
}else{
  width=20
  height=8
  design="aaabbb"
}

patch_t = patchwork::wrap_plots(
  plot_tool +labs(title=paste0(tool, " clusters")),
  scatter_process + labs(title="Process clusters"),
  design=design)&
  patchwork::plot_annotation(tag_levels="a", 
                             title = paste0(tool, "_", spn, "_", simulation_id, " - interpreted"),
                             subtitle = paste0("NMI = ", nmi_interpreted, "\nARI = ", ari_interpreted)) &
  theme(plot.tag=element_text(size=12, face="bold"),
        plot.title = element_text(size=12, face="bold", hjust=0.5))

ggsave(get_plots_path(save_path, tool, spn, simulation_id, plot_name="interpreted"), plot = patch_t,
       device=png, width=width, height=height, units="cm")

ggsave(get_plots_path_shared(main_path, tool, spn, simulation_id, plot_name="interpreted"), plot = patch_t,
       device=png, width=width, height=height, units="cm")

### Write NMI values ####

file_path = paste0(save_path, "/results/metrics.csv")
columns = c("spn", "coverage", "purity", "vcf_caller", "cna_caller", "tool", "nmi_general", "nmi_interpreted", "ari_general","ari_interpreted")

if (file.exists(file_path)) {
  metrics = read.csv(file_path, stringsAsFactors = FALSE)
  
  # Ensure all columns exist
  missing_cols <- setdiff(columns, names(metrics))
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      metrics[[col]] <- NA
    }
  }
  
} else {
  metrics = data.frame(matrix(ncol = length(columns), nrow = 0))
  colnames(metrics) = columns
}

df = data.frame(
  spn = spn,
  coverage = coverage,
  purity = purity,
  vcf_caller = vcf_caller,
  cna_caller = cna_caller,
  tool = tool,
  nmi_general = nmi_complete,
  nmi_interpreted = nmi_interpreted, 
  ari_general = ari_complete, 
  ari_interpreted = ari_interpreted
)

key_cols = c("spn", "coverage", "purity", "vcf_caller", "cna_caller","tool")
# Check if a row with the same combination exists
idx = which(
  metrics$spn == spn &
    metrics$coverage == coverage &
    metrics$purity == purity &
    metrics$vcf_caller == vcf_caller &
    metrics$cna_caller == cna_caller&
    metrics$tool == tool
)

# If exists, overwrite it; otherwise append
if (length(idx) > 0) {
  metrics[idx, ] = df
} else {
  metrics = rbind(metrics, df)
}

metrics$nmi_general <- round(metrics$nmi_general, 4)
metrics$nmi_interpreted <- round(metrics$nmi_interpreted, 4)
metrics$ari_general <- round(metrics$ari_general, 4)
metrics$ari_interpreted <- round(metrics$ari_interpreted, 4)

write.csv(metrics, file_path, row.names = FALSE)
write.csv(metrics, file.path(main_path,"subclonal/metrics_tables/metrics.csv"), row.names= FALSE)
          