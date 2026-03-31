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

tool = 'mobster'

# coverage_list = c(50,100, 150)
# purity_list = c(0.3, 0.6, 0.9)
# vcf_caller_list = c("mutect2")
# cna_caller_list = c("ascat")
coverage = 100
purity = 0.9
vcf_caller = "mutect2"
cna_caller = 'ascat'

spns = c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN05','SPN06', 'SPN07')

# combs = expand.grid(coverage=coverage_list,
#                     purity=purity_list,
#                     vcf_caller=vcf_caller_list,
#                     cna_caller=cna_caller_list)

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")
source(file.path(save_path, "utils_plots_final.R"))
source(file.path(save_path, "generate_table_main.R"))

# for(spn in spns){
if (!is.na(j)) {
  spn = spns[j]
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
    
    for(sample in sample_names){
      sample = sub("^([A-Z0-9]+)_\\1_", "\\1_", sample)
      sample_forest = sample_forest_full$get_subforest_for(sample) #"SPN04_1.1"
      
      if(spn %in% c('SPN03', 'SPN04')){
        phylo_forest = load_phylogenetic_forest(paste0(subforest_path, "/subforest_", spn, "_", stringr::str_remove(sample, ".*_"), '.sff'))
      }else{
        phylo_forest = phylo_forest_full$get_subforest_for(sample)
      }
      
      color_palette_process_uni = RColorBrewer::brewer.pal(n = 8, name = "Set1") 
      names(color_palette_process_uni) <- c('Clonal', paste0('Clone ', 1:7))
      color_palette_process_uni['Subclonal'] = 'gray70'
      
      # color_palette_tool = c(
      #   C1 = "red",
      #   C2 = "blue",
      #   C3 = "green",
      #   C4 = "orange",
      #   Other = "grey",
      #   Tail = "gainsboro"
      # )
      
      
      # Extract Set1 palette
      set1 <- brewer.pal(9, "Set1")
      
      color_palette_tool <- c(
        C1 = set1[1],  # red
        C2 = set1[2],  # blue
        C3 = set1[3],  # green
        C4 = set1[5],  # orange
        Other = "grey",
        Tail = "gainsboro"
      )
      
      plot_sticks = plot_mutations_on_subtree(final_table,
                                             process_seq,
                                             sample_forest,
                                             phylo_forest,
                                             color_palette_tool,
                                             color_palette_process_uni)
     
      plot_sticks_tool = plot_sticks[[1]]
      plot_sticks_process = plot_sticks[[2]]
      
      patch_t = patchwork::wrap_plots(
        plot_sticks_tool +labs(title='Mobster tree'),
        plot_sticks_process +labs(title='Process tree'),
        design='ab')&
        patchwork::plot_annotation(tag_levels="a", 
                                   title = paste0(tool, "_", sample)) &
        theme(plot.tag=element_text(size=12, face="bold"),
              plot.title = element_text(size=12, face="bold", hjust=0.5))
   
      width = 25
      height = 15
      plot_name = 'subtrees/png/process'
      
      ggsave(get_plots_path_process(save_path, tool, sample, simulation_id, plot_name=plot_name), plot = patch_t,
             device="png", width=width, height=height, units="cm")
    # }
    
    # ggsave(get_plots_path_process_shared(main_path, tool, spn, simulation_id, plot_name=plot_name), plot = patch_t,
    #        device="png", width=width, height=height, units="cm")
    
    }
  
  width = 25
  height = 15
  plot_name = 'subtrees/process'
  
  # Save here a unique pdf file with all the previous plots (i.e. at the end I want one pdf for each spn)
  # png_dir <- dirname(
  #   get_plots_path_process(save_path, tool, spn, simulation_id, plot_name = plot_name)
  # )
  png_dir = paste0(save_path, 'plots/subtrees/png/')
  
  pdf_dir <- file.path(save_path, "plots/subtrees/")
  
  pdf_file = file.path(pdf_dir, paste0(spn, ".pdf"))
  
  pdf(pdf_file, width = width, height = height)
  
  
  png_files <- list.files(
    png_dir,
    pattern = paste0("^process_", spn, ".*\\.png$"),
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
