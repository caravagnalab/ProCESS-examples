library(tidyverse)
library(randnet)
library(scales)
library(mobster)
library(patchwork)
library(ggrepel)
library(RColorBrewer)

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

source(file.path(github_path, "validation/Subclonal_deconvolution/utils_plots_final.R"))

coverage_list = c(50,100,150)
purity_list = c(0.3,0.6,0.9)

vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04','SPN05', 'SPN06', 'SPN07')
# spn_list = c('SPN01')

# tool = 'mobster_univariate'
tool = 'mobster'
combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller=cna_caller_list)

i = 1

for(spn in spn_list){
  print(spn)
  for(i in 1:nrow(combs)){
    coverage = combs[i, "coverage"]
    purity = combs[i, "purity"]
    vcf_caller = combs[i, "vcf_caller"]
    cna_caller = combs[i, "cna_caller"]
    # spn = combs[i, "spn"]
    
    simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
    
    print(paste0(spn, "_", coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller))
    
    # Take mobster table
    tool = 'mobster_univariate'
    
    table = tryCatch(
      readRDS(file.path(main_path, "validation_subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds"))),
      error = function(e) {
        message("Skipping simulation_id: ", simulation_id,
                " (", e$message, ")")
        return(NULL)
      }
    )
    
    if (is.null(table)) {
      next
    }
    
    # table = readRDS(file.path(main_path, "validation_subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds")))
    # table = readRDS(file.path(save_path, "/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds")))
    
    table = table %>%
        group_by(cluster_id_process, sample_id) %>%
        mutate(ccf_process=mean(ccf_process, na.rm=TRUE)) %>%
        ungroup()

    table = table %>% group_by(cluster_id_process, sample_id) %>%
      mutate(is_clonal_process=replace(FALSE, ccf_process > 0.9, TRUE)) %>% ungroup() %>%
      mutate(cluster_id_process_full = cluster_id_process) %>%
      mutate(cluster_id_process=replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal'))
      
    sample_names = table$sample_id %>% unique()
    
    tool = 'mobster'
    main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
    path_m = file.path(main_path,spn, "tumourevo", simulation_id, "subclonal_deconvolution", tool, "SCOUT", spn)
    samples = list.dirs(path_m, recursive=F, full.names=F)  # get_sample_names(spn, base_path=main_path)
    
    print(samples)
    
    plots_mobster = lapply(samples, function(sample_name) {
      # sample_name = samples[[1]]
      print(file.path(path_m,
                      paste0(sample_name),
                      paste0("SCOUT_", spn, "_", sample_name, "_mobsterbest_fit.rds")))
      
      obj = readRDS(file.path(path_m,
                              paste0(sample_name),
                              paste0("SCOUT_", spn, "_", sample_name, "_mobster_best_fit.rds")))
      
      plot(obj)+ggtitle(sample_name)
    })
    
    set3 = brewer.pal(12, "Set3")
    extra_colors = c(
      "#8B5A2B",  # brown
      "#000000",  # black
      "#7F0000"   # dark red
    )
    set3_extended = c(set3, extra_colors)
    set3_extended[2] = 'goldenrod'
    
    color_palette = set3_extended[1:length(unique(table$cluster_id_process))] %>% 
      setNames(str_sort(unique(table$cluster_id_process), numeric=T))
    
    color_palette["Subclonal"] = "#cccccc"
    
    plots_process = lapply(samples, function(sample_name) {
      
      table_sample = table %>% filter(sample_id == sample_name)
      
      cluster_order = table_sample %>%
        dplyr::count(cluster_id_process) %>% # Count the number of data points per cluster
        arrange(desc(n)) %>% # Sort clusters by size (descending order)
        pull(cluster_id_process)
      
      table_sample = table_sample %>%
        mutate(cluster_id_process = factor(cluster_id_process, levels = cluster_order))
      
      # Here plot the marginal of process
      table_sample %>%
        filter(vaf_process>0) %>% 
        ggplot() +
        xlim(0,1)+
        geom_histogram(aes(x=vaf_process, fill=cluster_id_process), position="identity", alpha=1, bins=100) +
        ylab("") +
        ggtitle(sample_name)+
        scale_fill_manual(values=color_palette, name="Cluster") +
        theme(legend.position="bottom")#+
        # theme_minimal()
      
    })
    
    
    set1 = brewer.pal(8, "Set1")
    color_palette_interpreted = set1[1:length(unique(table$cluster_id_tool_interpreted))] %>% 
      setNames(str_sort(unique(table$cluster_id_tool_interpreted), numeric=T))
    
    color_palette_interpreted["Tail"] = "gainsboro"
    color_palette_interpreted["Other"] = "gray"
    color_palette_interpreted['C1'] = set1[1]
    color_palette_interpreted['C2'] = set1[2]
    color_palette_interpreted['C3'] = set1[3]
    
    plots_tool_interpreted = lapply(samples, function(sample_name) {
      
      table_sample = table %>% filter(sample_id == sample_name)
      
      cluster_order = table_sample %>%
        dplyr::count(cluster_id_tool_interpreted) %>% # Count the number of data points per cluster
        arrange(desc(n)) %>% # Sort clusters by size (descending order)
        pull(cluster_id_tool_interpreted)
      
      table_sample = table_sample %>%
        mutate(cluster_id_tool_interpreted = factor(cluster_id_tool_interpreted, levels = cluster_order))
      
      # Here plot the marginal of process
      table_sample %>%
        filter(vaf_tool>0) %>% 
        ggplot() +
        xlim(0,1)+
        geom_histogram(aes(x=vaf_tool, fill=cluster_id_tool_interpreted), position="identity", alpha=1, bins=100) +
        ylab("") +
        ggtitle(sample_name)+
        scale_fill_manual(values=color_palette_interpreted, name="Cluster") +
        scale_color_manual(values = color_palette_interpreted, guide = "none") +
        theme(legend.position="bottom")+
        
        ggrepel::geom_label_repel(
          data = subset(table_sample, is_driver_tool == TRUE),
          aes(x = vaf_tool, y = 0, label = driver_label_tool, color = cluster_id_tool_interpreted),
          size = 3,
          show.legend = FALSE,
          max.overlaps = Inf,
          nudge_y = 300,
          nudge_x = 0.5
        )
      
    })
    
    # color_palette_interpreted_driver = set1[1:length(unique(table$cluster_id_tool_interpreted_driver))+1] %>% 
    #   setNames(str_sort(unique(table$cluster_id_tool_interpreted_driver), numeric=T))
    
    color_palette_interpreted_driver = color_palette_interpreted
    color_palette_interpreted_driver['C1'] = set1[1]
    color_palette_interpreted_driver['C2'] = set1[2]
    color_palette_interpreted_driver['C3'] = set1[3]
    
    color_palette_interpreted_driver["Tail"] = "gainsboro"
    color_palette_interpreted_driver["Other"] = "gray"
    
    plots_tool_interpreted_driver = lapply(samples, function(sample_name) {
      
      table_sample = table %>% filter(sample_id == sample_name)
      
      cluster_order = table_sample %>%
        dplyr::count(cluster_id_tool_interpreted_driver) %>% # Count the number of data points per cluster
        arrange(desc(n)) %>% # Sort clusters by size (descending order)
        pull(cluster_id_tool_interpreted_driver)
      
      table_sample = table_sample %>%
        mutate(cluster_id_tool_interpreted_driver = factor(cluster_id_tool_interpreted_driver, levels = cluster_order))
      
      # Here plot the marginal of process
      table_sample %>%
        filter(vaf_tool>0) %>% 
        ggplot() +
        xlim(0,1)+
        geom_histogram(aes(x=vaf_tool, fill=cluster_id_tool_interpreted_driver), position="identity", alpha=1, bins=100) +
        ylab("") +
        ggtitle(sample_name)+
        scale_fill_manual(values=color_palette_interpreted_driver, name="Cluster") +
        scale_color_manual(values = color_palette_interpreted_driver, guide = "none") +
        theme(legend.position="bottom")+
        
        ggrepel::geom_label_repel(
          data = subset(table_sample, is_driver_process == TRUE),
          aes(x = vaf_tool, y = 0, label = driver_label_process, color = cluster_id_tool_interpreted),
          size = 3,
          show.legend = FALSE,
          max.overlaps = Inf,
          nudge_y = 300,
          nudge_x = 0.5
        )
      
    })
    
    ncol_max <- min(3, length(samples))
    
    patch_iter = (
      (wrap_plots(plots_mobster, ncol = ncol_max) +
          plot_annotation(title = "Mobster")) /
        (wrap_plots(plots_tool_interpreted, ncol = ncol_max) +
            plot_annotation(title = "Tool Interpreted")) /
        (wrap_plots(plots_tool_interpreted_driver, ncol = ncol_max) +
            plot_annotation(title = "Tool Interpreted Driver"))) +
      plot_annotation(
        title = paste0(
          spn, " | ",
          coverage, "x | ",
          purity, "p | ",
          vcf_caller, " | ",
          cna_caller
        )
      )
    
    png_file <- file.path(
      save_path,
      "plots/marginals_process_mobster/png/",
      paste0(spn, "_", simulation_id, ".png")
    )
    ggsave(
      filename = png_file,
      plot = patch_iter,
      # width = 4 * length(samples),
      width = 4*ncol_max,
      height = 8*ceiling(length(samples) / ncol_max),
      units = "in",
      dpi = 300
    )
  }


  png_files <- list.files(
    file.path(save_path,"plots/marginals_process_mobster/png/"),
    pattern = paste0(spn, ".*\\.png$"),
    full.names = TRUE
  )
  png_files <- sort(png_files)
  
  pdf_file = file.path(save_path, "plots/marginals_process_mobster/", paste0(tool, "_", spn, ".pdf"))
  
  
  pdf(pdf_file, width = 10, height = 12)
  
  for (f in png_files) {
    img <- png::readPNG(f)
    grid::grid.newpage()
    grid::grid.raster(img)
  }
  
  dev.off()
}


