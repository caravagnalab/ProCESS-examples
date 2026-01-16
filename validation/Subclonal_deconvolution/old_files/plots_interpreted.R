library(ggplot2)
library(tidyverse)
library(ProCESS)
library(ggpubr)
library(patchwork)
library(randnet)
library(scales)
library(ggrepel)

# coverage_list = c(50, 100, 150, 200)
# purity_list = c(0.3, 0.6, 0.9)
# vcf_caller_list = c("mutect2", "strelka", "freebayes")
# cna_caller_list = c("ascat", "sequenza", "battenberg")
# spn_list = paste("SPN", 3:7, sep="0")

spn = 'SPN01'
purity=0.6
vcf_caller = "mutect2"
cna_caller = "ascat"
coverage = 150


spns = c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN06', 'SPN07')
spns = c('SPN03', 'SPN04')

coverage_list = c(100)
coverage_list = c(50,150)
purity_list = c(0.3, 0.6, 0.9)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN06', 'SPN07')
tool = 'viber'

combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller=cna_caller_list,
                    spn=spn_list)

# for(spn in spns){
for(i in 1:nrow(combs)){
  
  # if(spn=='SPN02'){
  #   coverage=100
  # }else{
  #   coverage=100
  # }
  coverage = combs[i, "coverage"]
  purity = combs[i, "purity"]
  vcf_caller = combs[i, "vcf_caller"]
  cna_caller = combs[i, "cna_caller"]
  spn = combs[i, "spn"]
  
  simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
  
simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

source(file.path(save_path, "utils_plots_final.R"))
source(file.path(save_path, "generate_table_main.R"))

# drivers_table = readRDS(file.path(main_path,"/drivers/all_drivers.rds"))

# True process drivers
true_drivers_table = readRDS(file.path(main_path,"drivers", spn, "process_drivers.rds"))

true_drivers_table = true_drivers_table %>% filter(SPN==spn)  %>%
  mutate(mutation_id=paste0(spn, ":", chr, ":", start,  ":", alt)) %>% 
  select(mutation_id, code)

# Get process table
mut_process = get_mutations(spn=spn, type="tumour", coverage=coverage, purity=purity)
table_process = readRDS(get_table_path(save_path, 'process_new', spn, simulation_id)) # process table in folder tables/

# Join process table with drivers
  # now in process_table we have a column "code" with the drivers gene names
table_process = table_process %>% left_join(true_drivers_table, by = 'mutation_id') %>% 
  select(-driver_label_process) %>% 
  mutate(driver_label_process=code)

tool = 'viber'
# SPN03:9:136496197: deve diventare SPN03:9:136496196:C
if(spn=='SPN03' & tool != 'process_viber'){
  table_process = table_process %>% mutate(mutation_id = if_else(
    mutation_id=="SPN03:9:136496197:",
    'SPN03:9:136496196:C',
    mutation_id
  ))
}

table_tool = readRDS(get_table_path(save_path, tool, spn, simulation_id))

if(tool =='viber_heuristics'){
  table_tool = table_tool %>%
    mutate(cluster_id_tool = replace_na(cluster_id_tool, 'NA'))
}
if(tool =='process_viber'){
  table_tool = table_tool %>%
    mutate(sample_id = paste0(spn,"_",sample_id))
}

join_table_tool = table_tool %>% left_join(table_process) # keep all mut in tool
join_table_process = join_table_tool %>% filter(!is.na(cluster_id_process)) # only mutations present in both

nmi_complete = randnet::NMI(as.factor(join_table_process$cluster_id_tool),
                    as.factor(join_table_process$cluster_id_process))

ari_complete = aricode::ARI(as.factor(join_table_process$cluster_id_tool), 
            as.factor(join_table_process$cluster_id_process))

### Find cluster/driver in tool and add column cluster_id_tool_interpreted
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

### Find cluster/driver in process and add column cluster_id_tool_interpreted_driver
driver_clusters_process = join_table_tool %>%
  distinct(cluster_id_tool, is_driver_process) %>% 
  filter(is_driver_process == TRUE) %>% 
  pull(cluster_id_tool)

final_table = final_table %>%
  mutate(cluster_id_tool_interpreted_driver = if_else(
    !(cluster_id_tool %in% driver_clusters_process),
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
# color_palette = color_palette_tool
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
  )) %>%
  mutate(cluster_id_tool_interpreted_driver = if_else(
    !(cluster_id_tool %in% driver_clusters_process),
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
# table_to_save = final_table_interpreted %>% select(patient_id, sample_id,coverage, purity, 
#                                                    tool, mutation_id, driver_label_tool,
#                                                    is_driver_tool, cluster_id_tool, vaf_tool, 
#                                                    is_driver_process, cluster_id_process,
#                                                    vaf_process, driver_label_process, 
#                                                    cluster_id_tool_interpreted)

table_to_save = final_table_interpreted 
saveRDS(table_to_save, file.path(main_path, "subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds")))
saveRDS(table_to_save, file.path(save_path, "tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds")))
# }

if(spn=='SPN03'){
  width=40
  height=50
  design="aaaaaa\nbbbbbb"
}else if(spn=='SPN01'){
  width=40
  height=25
  design="aaaaaa\nbbbbbb"
}else if(spn=='SPN06' | spn=='SPN07'){
  width=40
  height=70
  design="aaaaaa\nbbbbbb"
}else if(spn=='SPN02' | spn=='SPN04'){
  width=20
  height=12
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

if(startsWith(tool, 'process_')){
  # CAN BE CHANGED
  plot_name = 'interpreted_drivers/interpreted'
}else{
  plot_name = 'interpreted/interpreted'
}
ggsave(get_plots_path(save_path, tool, spn, simulation_id, plot_name=plot_name), plot = patch_t,
       device="png", width=width, height=height, units="cm")

ggsave(get_plots_path_shared(main_path, tool, spn, simulation_id, plot_name=plot_name), plot = patch_t,
       device="png", width=width, height=height, units="cm")

}

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
          