library(ggplot2)
library(tidyverse)
library(ProCESS)
library(ggpubr)
library(patchwork)
library(randnet)
library(scales)

spn = 'SPN04'
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

# ------ #
# Do not run again if already loaded for the same combination of parameters
sample_forest = load_sample_forest(get_sample_forest(spn))
phylo_forest = load_phylogenetic_forest(get_phylo_forest(spn=spn))

mut_process = get_mutations(spn=spn, type="tumour", coverage=coverage, purity=purity)
seq_results = readRDS(mut_process)

process_seq = seq_results %>%
  mutate(mutation_id=paste0(spn, ":", chr, ":", chr_pos,  ":", alt))
# ------ #

### Plot process labels ####
table_process = readRDS(get_table_path(save_path, 'process', spn, simulation_id)) # process table in folder tables/

tool = 'viber_heuristics'
table_tool = readRDS(get_table_path(save_path, tool, spn, simulation_id))
if(tool =='viber_heuristics'){
  table_tool = table_tool %>%
    mutate(cluster_id_tool = replace_na(cluster_id_tool, 'NA'))
}

join_table_tool = table_tool %>% left_join(table_process) # keep all mut in tool
# join_table_process = table_process %>% left_join(table_tool) # keep all mut in process
join_table_process = join_table_tool %>% filter(!is.na(cluster_id_process)) # only mutations present in both


color_palette_tool = c(
  "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00","#a65628",
  "#FFD700",  "#999999", "#000000", "#f781bf", # First 10 colors (Set1)  
  "#46f0f0", "#f032e6", "#bcf60c", "#fabed4", "#008080", "#e6beff",  
  "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",  
  "#000075", "#808080", "#d3a6f3", "#ff9cdd", "#73d7b0"  ) %>% setNames(str_sort(unique(table_tool$cluster_id_tool), numeric=T))

if(tool =='viber_heuristics'){
  color_palette_tool["NA"] <- "#cccccc"
}
color_palette_process = hue_pal()(length(unique(table_process$cluster_id_process))) %>% 
  setNames(str_sort(unique(table_process$cluster_id_process), numeric=T))

sample_names = sort(unique(table_process$sample_id))
scatter_process = plot_scatter_process(join_table_process, sample_names, color_palette_process)
scatter_process


### Plot mutations on tree #### 

plot_tool = plot_scatter_tool(join_table_tool, color_palette_tool, sample_names)
plot_tool

plot_sticks = plot_mutations_on_tree(join_table_tool, 
                                     process_seq, 
                                     sample_forest,
                                     phylo_forest,
                                     color_palette_tool, 
                                     color_palette_process)
plot_sticks_tool = plot_sticks[[1]]
plot_sticks_process = plot_sticks[[2]]

# new_table = table_tool %>% full_join(table_process)
join_table_tool2 = join_table_tool %>% filter(!is.na(cluster_id_process))

nmi <- randnet::NMI(as.factor(join_table_tool2$cluster_id_tool),
                    as.factor(join_table_tool2$cluster_id_process))

ari = aricode::ARI(as.factor(join_table_tool2$cluster_id_tool), 
            as.factor(join_table_tool2$cluster_id_process))

if(spn=='SPN03'){
  width=40
  height=50
  design="aaaaaa\ncccccc"
}else{
  width=40
  height=50
  design="aa####\ncccccc\ncccccc"
}

patch_t = patchwork::wrap_plots(
  plot_tool +labs(title=paste0(tool, " clusters")),
  plot_sticks_tool + labs(title=paste0(tool, " tree")),
  design=design)&
  patchwork::plot_annotation(tag_levels="a", 
                             title = paste0(tool, "_", spn, "_", simulation_id),
                             subtitle = paste0("NMI = ", nmi, "\nARI = ", ari)) &
  theme(plot.tag=element_text(size=12, face="bold"),
        plot.title = element_text(size=12, face="bold", hjust=0.5))

ggsave(get_plots_path(save_path, tool, spn, simulation_id, plot_name="final"), plot = patch_t,
       device=png, width=width, height=height, units="cm")

patch_pr = patchwork::wrap_plots(
  scatter_process + theme(legend.position="bottom") + labs(title="Process clusters"), 
  plot_sticks_process + labs(title="Process tree"),
  design=design)&
  patchwork::plot_annotation(tag_levels="a", title = paste0('Process', "_", spn, "_", simulation_id)) &
  theme(plot.tag=element_text(size=12, face="bold"),
        plot.title = element_text(size=12, face="bold", hjust=0.5))

ggsave(get_plots_path(save_path, 'process', spn, simulation_id, plot_name=paste0("final_",tool)), plot = patch_pr,
       device=png, width=width, height=height, units="cm")

