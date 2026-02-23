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


spn = 'SPN03'
purity=0.9
vcf_caller = "mutect2"
cna_caller = "ascat"
coverage = 100

spns = c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN06', 'SPN07')

coverage_list = c(50,100, 150)
purity_list = c(0.3, 0.6, 0.9)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04','SPN05', 'SPN06', 'SPN07')
spn_list = c('SPN02', 'SPN03')

tool = 'viber'

combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    # spn=spn_list,
                    cna_caller=cna_caller_list)

interpreted_driver=T

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

source(file.path(save_path, "utils_plots_final.R"))
source(file.path(save_path, "generate_table_main.R"))

for(spn in spn_list){
  
  for(i in 1:nrow(combs)){
    
    coverage = combs[i, "coverage"]
    purity = combs[i, "purity"]
    vcf_caller = combs[i, "vcf_caller"]
    cna_caller = combs[i, "cna_caller"]
    # spn = combs[i, "spn"]
    
    simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
    print(paste0(spn,'_', simulation_id))
    
    # Get interpreted table
    # final_table = readRDS(file.path(main_path, "validation_subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds"))) # process table in folder tables/
    
    final_table = tryCatch(
      readRDS(file.path(main_path, "validation_subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds"))),
      error = function(e) {
        message("Skipping simulation_id: ", simulation_id,
                " (", e$message, ")")
        return(NULL)
      }
    )
    
    if (is.null(final_table)) {
      next
    }
    nmi_complete = randnet::NMI(as.factor(final_table$cluster_id_tool),
                        as.factor(final_table$cluster_id_process))
    
    ari_complete = aricode::ARI(as.factor(final_table$cluster_id_tool), 
                as.factor(final_table$cluster_id_process))
  
    sample_names = sort(unique(final_table$sample_id))
    
    ### Plot tool with new labels ####
    color_palette_tool = c(
      "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00","#a65628",
      "#FFD700", "#000000", "#f781bf", # First 10 colors (Set1)  
      "#46f0f0", "#f032e6", "#bcf60c", "#fabed4", "#008080", "#e6beff",  
      "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",  
      "#000075", "#808080", "#d3a6f3", "#ff9cdd", "#73d7b0"  ) %>% setNames(str_sort(unique(final_table$cluster_id_tool), numeric=T))
    
    color_palette_tool["Subclonal"] = "#cccccc"
  
    # Scatterplot tool
    
    plot_tool = plot_scatter_tool(final_table, color_palette_tool, sample_names, type ='interpreted', vertical = T)
    # plot_tool
    
    plot_tool_raw = plot_scatter_tool(final_table, color_palette_tool, sample_names, type ='original', vertical = T)
    # plot_tool_raw
    
    plot_tool_interpreted_driver = plot_scatter_tool(final_table, color_palette_tool, sample_names, type ='interpreted_driver', vertical = T)
    plot_tool_interpreted_driver
    
    
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
  
    # a = final_table %>% mutate( cluster_id_process=  cluster_id_process_full)
    # color_palette_process = hue_pal()(length(unique(a$cluster_id_process))) %>% 
    #   setNames(str_sort(unique(a$cluster_id_process), numeric=T))
    # color_palette_process["Subclonal"] = "#cccccc"
    # scatter_process = plot_scatter_process(a, sample_names, color_palette_process, driver=T)
    # scatter_process
    
    scatter_process = plot_scatter_process(final_table, sample_names, color_palette_process, driver=T, vertical = T)
    # scatter_process
  
    nmi_interpreted = randnet::NMI(as.factor(final_table$cluster_id_tool_interpreted),
                                 as.factor(final_table$cluster_id_process))
    
    ari_interpreted = aricode::ARI(as.factor(final_table$cluster_id_tool_interpreted),
                                as.factor(final_table$cluster_id_process))
  
  if(interpreted_driver==T){
    width=40
    design="abcd"
    patch_t = patchwork::wrap_plots(
      plot_tool_raw +labs(title=paste0(tool, " blind clusters")),
      plot_tool +labs(title=paste0(tool, " interpreted clusters")),
      plot_tool_interpreted_driver+labs(title=paste0(tool, " interpreted driver clusters")),
      scatter_process + labs(title="Process clusters"),
      design=design, guides = 'collect') &
      patchwork::plot_annotation(tag_levels="a", 
                                 title = paste0(tool, "_", spn, "_", simulation_id),
                                 subtitle = paste0("NMI blind = ", nmi_complete, "\nNMI interpreted = ", nmi_interpreted)) &
      theme(plot.tag=element_text(size=12, face="bold"),
            plot.title = element_text(size=12, face="bold", hjust=0.5), legend.position = 'right')
    
    plot_name = 'multivariates_png/process'
    
  }else if(interpreted_driver==F){
      width=30
      design="abc"
      patch_t = patchwork::wrap_plots(
        plot_tool_raw +labs(title=paste0(tool, " blind clusters")),
        plot_tool +labs(title=paste0(tool, " interpreted clusters")),
        scatter_process + labs(title="Process clusters"),
        design=design, guides = 'collect') &
        patchwork::plot_annotation(tag_levels="a", 
                                   title = paste0(tool, "_", spn, "_", simulation_id),
                                   subtitle = paste0("NMI blind = ", nmi_complete, "\nNMI interpreted = ", nmi_interpreted)) &
        theme(plot.tag=element_text(size=12, face="bold"),
              plot.title = element_text(size=12, face="bold", hjust=0.5), legend.position = 'right')
      
      plot_name = 'multivariates_png/process'
  }
    
  
  if(spn=='SPN03'){
    height=50
  }else if(spn=='SPN01'| spn=='SPN05'){
    height=25
  }else if(spn=='SPN06' | spn=='SPN07'){
    # height=70
    height=70
  }else if(spn=='SPN02' | spn=='SPN04'){
    height=10
  }
  # print(plot_name)
  ggsave(get_plots_path(save_path, tool, spn, simulation_id, plot_name=plot_name), plot = patch_t,
         device="png", width=width, height=height, units="cm")
  
  # ggsave(get_plots_path_shared(main_path, tool, spn, simulation_id, plot_name=plot_name), plot = patch_t,
  #        device="png", width=width, height=height, units="cm")
  
   
}
  
  # Save here a unique pdf file with all the previous plots (i.e. at the end I want one pdf for each spn)
  png_dir <- dirname(
    get_plots_path(save_path, tool, spn, simulation_id, plot_name = plot_name)
  )
  
  pdf_dir <- file.path(save_path, "plots/multivariates/")
  # dir.create(pdf_dir, showWarnings = FALSE, recursive = TRUE)
  
  pdf_file = file.path(pdf_dir, paste0(tool, "_", spn, ".pdf"))
  
  # Open PDF device
  pdf(pdf_file, width = width, height = height)
  
  # Find only this SPN's images
  png_files <- list.files(
    png_dir,
    pattern = paste0("^process_", tool, "_", spn, ".*\\.png$"),
    full.names = TRUE
  )
  
  # Sort them
  png_files <- sort(png_files)
  
  # # Read and write to PDF
  # img_list <- image_read(png_files)
  # image_write(img_list, path = pdf_file, format = "pdf")
  
  for (f in png_files) {
    img <- png::readPNG(f)
    grid::grid.newpage()
    grid::grid.raster(img)
  }
  
  # Close PDF
  dev.off()
  
  message("Saved PDF: ", pdf_file)
}


