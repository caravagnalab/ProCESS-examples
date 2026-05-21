# Plot multivariate of:
# tool blind
# tool interpreted
# process
library(ggplot2)
library(tidyverse)
library(ProCESS)
library(ggpubr)
library(patchwork)
library(randnet)
library(scales)
library(ggrepel)
library(RColorBrewer)
library(png)
library(grid)

tool = 'mobster'
coverage = 150
purity = 0.9
vcf_caller = "mutect2"
cna_caller = 'ascat'

spn = 'SPN04'

source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/report/plotting/utils.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")
source(file.path(save_path, "utils_plots_final.R"))
source(file.path(save_path, "generate_table_main.R"))

# for(spn in spns){
subforest_path = file.path("/orfeo/cephfs/scratch/cdslab/shared/SCOUT", spn, "process")
sample_forest_full = load_sample_forest(get_sample_forest(spn))
phylo_forest_full = load_phylogenetic_forest(get_phylo_forest(spn=spn))

if(spn == 'SPN05'){
  purity = 0.3
}

if(spn == 'SPN07'){
  coverage = 150
}

mut_process = get_mutations(spn=spn, type="tumour", coverage=coverage, purity=purity)

# Get interpreted table
simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)

# final_table = readRDS(file.path(main_path, "validation_subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds"))) # process table in folder tables/

final_table = tryCatch(
  readRDS(file.path(main_path, "validation_subclonal_new/tables_interpreted", paste0(tool, "_univariate_", spn, "_", simulation_id, ".rds"))),
  error = function(e) {
    message("Skipping simulation_id: ", simulation_id,
            " (", e$message, ")")
    return(NULL)
  }
)

if (is.null(final_table)) {
  next
}

mut_process = get_mutations(spn=spn, type="tumour", coverage=coverage, purity=purity)

seq_results = readRDS(mut_process)

process_seq = seq_results %>%
  mutate(mutation_id=paste0(spn, ":", chr, ":", chr_pos,  ":", alt))

sample_names = sort(unique(final_table$sample_id))
sample = "SPN04_SPN04_2.1"

sample = sub("^([A-Z0-9]+)_\\1_", "\\1_", sample)
sample_forest = sample_forest_full$get_subforest_for(sample) #"SPN04_1.1"

if(spn %in% c('SPN03', 'SPN04')){
  phylo_forest = load_phylogenetic_forest(paste0(subforest_path, "/subforest_", spn, "_", stringr::str_remove(sample, ".*_"), '.sff'))
}else{
  phylo_forest = phylo_forest_full$get_subforest_for(sample)
}



colors_cluster = c('indianred', 
                   'steelblue', 
                   'palevioletred',                   
                   'goldenrod', 
                   'darkorange3', 
                   'palevioletred', 
                   'mediumpurple', 
                   'cornsilk4', 
                   'olivedrab3', 
                   'steelblue4', 
                   'indianred4',
                   'aquamarine3',
                   'saddlebrown',
                   'deeppink2',
                   'cornflowerblue',
                   'black')
names(colors_cluster) = paste0('C',0:15)
tail = c('Tail' = 'gray70')
colors_cluster = c(colors_cluster, tail)

color_palette_process_uni = RColorBrewer::brewer.pal(n = 8, name = "Set1") 
names(color_palette_process_uni) <- c('Clonal', paste0('Clone ', 1:7))
color_palette_process_uni['Subclonal'] = 'gray70'


plot_sticks = plot_mutations_on_subtree(final_table,
                                        process_seq,
                                        sample_forest,
                                        phylo_forest,
                                        colors_cluster,
                                        color_palette_process_uni)

plot_sticks_tool = plot_sticks[[1]] + my_ggplot_theme() + theme(legend.position = 'none')
plot_sticks_process = plot_sticks[[2]] + my_ggplot_theme() + theme(legend.position = 'none')

patch_t = patchwork::wrap_plots(
  plot_sticks_tool ,
  plot_sticks_process ,
  design='b\na')

ggsave(plot = patch_t, filename = 'phylo_spn04_bias.png',
       width = 4.5, height = 5.5, units = 'in',
       dpi = 600)
