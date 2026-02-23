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

args = commandArgs(trailingOnly = TRUE)
cat(paste("\nArguments:", paste(args, collapse=", "), "\n"))

j = as.integer(args[1])
print(j)

coverage_list = c(50,100, 150)
purity_list = c(0.3, 0.6, 0.9)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spns = c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN06', 'SPN07')

tool = 'viber'

combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller=cna_caller_list)

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")
source(file.path(save_path, "utils_plots_final.R"))
source(file.path(save_path, "generate_table_main.R"))

# spn = 'SPN02'
# coverage = 50
# purity = 0.3
# 
# vcf_caller = "mutect2"
# cna_caller = "ascat"

# for(spn in spns){
if (!is.na(j)) {
  spn = spns[j]
  sample_forest = load_sample_forest(get_sample_forest(spn))
  phylo_forest = load_phylogenetic_forest(get_phylo_forest(spn=spn))
  
  for(i in 1:nrow(combs)){
    
    coverage = combs[i, "coverage"]
    purity = combs[i, "purity"]
    vcf_caller = combs[i, "vcf_caller"]
    cna_caller = combs[i, "cna_caller"]
    
    mut_process = get_mutations(spn=spn, type="tumour", coverage=coverage, purity=purity)
  
    # Get interpreted table
    simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
    final_table = readRDS(file.path(main_path, "validation_subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds"))) # process table in folder tables/
    
    # process_seq = final_table %>%
    #   select(patient_id, sample_id, mutation_id, cluster_id_process, vaf_process,cluster_id_tool, cluster_id_tool_interpreted) %>%
    #   pivot_wider(
    #     id_cols    = mutation_id,
    #     names_from = "sample_id",
    #     values_from = "vaf_process"
    #   )
    
    mut_process = get_mutations(spn=spn, type="tumour", coverage=coverage, purity=purity)
    
    seq_results = readRDS(mut_process)
    
    process_seq = seq_results %>%
      mutate(mutation_id=paste0(spn, ":", chr, ":", chr_pos,  ":", alt))
    
    sample_names = sort(unique(final_table$sample_id))
    
    ### Plot tool with new labels ####
    color_palette_tool = c(
      "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00","#a65628",
      "#FFD700", "#000000", "#f781bf", # First 10 colors (Set1)  
      "#46f0f0", "#f032e6", "#bcf60c", "#fabed4", "#008080", "#e6beff",  
      "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",  
      "#000075", "#808080", "#d3a6f3", "#ff9cdd", "#73d7b0"  ) %>% setNames(str_sort(unique(final_table$cluster_id_tool), numeric=T))
    
    color_palette_tool["Subclonal"] = "#cccccc"
    
    # Scatterplot process
    set3 <- brewer.pal(12, "Set3")
    
    extra_colors <- c(
      "#8B5A2B",  # brown
      "#000000",  # black
      "#7F0000"   # dark red
    )
    
    set3_extended <- c(set3, extra_colors)
    set3_extended[2] = 'goldenrod'
    
    color_palette_process = set3_extended[1:length(unique(final_table$cluster_id_process))] %>% 
      setNames(str_sort(unique(final_table$cluster_id_process), numeric=T))
    
    color_palette_process["Subclonal"] = "#cccccc"
    
  
    plot_sticks = plot_mutations_on_tree(final_table,
                                         process_seq,
                                         sample_forest,
                                         phylo_forest,
                                         color_palette_tool,
                                         color_palette_process)
    plot_sticks_tool = plot_sticks[[1]]
    plot_sticks_process = plot_sticks[[2]]
    
    
    patch_t = patchwork::wrap_plots(
      plot_sticks_tool +labs(title=paste0(tool, " tree")),
      plot_sticks_process +labs(title='Process tree'),
      design='ab')&
      patchwork::plot_annotation(tag_levels="a", 
                                 title = paste0(tool, "_", spn, "_", simulation_id)) &
      theme(plot.tag=element_text(size=12, face="bold"),
            plot.title = element_text(size=12, face="bold", hjust=0.5))
    
    if(spn=='SPN01' | spn=='SPN06'| spn=='SPN07'| spn=='SPN02'){
      width = 40
    }else{
      width = 30
    }
    
    height = 30
    plot_name = 'trees/process'
    
    ggsave(get_plots_path(save_path, tool, spn, simulation_id, plot_name=plot_name), plot = patch_t,
           device="png", width=width, height=height, units="cm")
    
    # ggsave(get_plots_path_shared(main_path, tool, spn, simulation_id, plot_name=plot_name), plot = patch_t,
    #        device="png", width=width, height=height, units="cm")
    
  }
  # Save here a unique pdf file with all the previous plots (i.e. at the end I want one pdf for each spn)
  png_dir <- dirname(
    get_plots_path(save_path, tool, spn, simulation_id, plot_name = plot_name)
  )
  
  pdf_dir <- file.path(save_path, "plots/trees/")
  
  pdf_file = file.path(pdf_dir, paste0(tool, "_", spn, ".pdf"))
  
  
  pdf(pdf_file, width = width, height = height)
  
  
  png_files <- list.files(
    png_dir,
    pattern = paste0("^process_", tool, "_", spn, ".*\\.png$"),
    full.names = TRUE
  )
  
  png_files <- sort(png_files)
  
  
  for (f in png_files) {
    img <- png::readPNG(f)
    grid::grid.newpage()
    grid::grid.raster(img)
  }
  
  dev.off()
  
  message("Saved PDF: ", pdf_file)
}
