library(tidyverse)
library(randnet)
library(scales)
library(mobster)
library(patchwork)

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")

source(file.path(save_path, "utils_plots_final.R"))

coverage_list = c(100)
purity_list = c(0.9)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN06', 'SPN07')
spn_list = c('SPN01')
# tool = 'mobster_univariate'
tool = 'mobster'
combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller=cna_caller_list,
                    spn=spn_list)


for(i in 1:nrow(combs)){
  coverage = combs[i, "coverage"]
  purity = combs[i, "purity"]
  vcf_caller = combs[i, "vcf_caller"]
  cna_caller = combs[i, "cna_caller"]
  spn = combs[i, "spn"]
  
  simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
  
  print(paste0(spn, "_", coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller))
  
  # tool = 'mobster_univariate'
  tool = 'mobster'
  table = readRDS(file.path(main_path, "subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds")))
  
  tool = 'mobster'
  table = table %>% 
      group_by(cluster_id_process, sample_id) %>%
      mutate(ccf_process=mean(ccf_process, na.rm=TRUE)) %>%
      ungroup()
    
  # SPN06, Clone 2
  table = table %>% group_by(cluster_id_process, sample_id) %>% 
    mutate(is_clonal_process=replace(FALSE, ccf_process > 0.9, TRUE)) %>% ungroup() %>% 
    mutate(cluster_id_process_full = cluster_id_process) %>% 
    mutate(cluster_id_process=replace(cluster_id_process_full, is_clonal_process==TRUE, 'Clonal'))
    
  sample_names = table$sample_id %>% unique()
  
  main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
  path_m = file.path(main_path,spn, "tumourevo", simulation_id, "subclonal_deconvolution", tool, "SCOUT", spn)
  samples = list.dirs(path_m, recursive=F, full.names=F)  # get_sample_names(spn, base_path=main_path)
  
  print(samples)
  
  plots_mobster = lapply(samples, function(sample_name) {
    # sample_name = samples[[1]]
    print(file.path(path_m,
                    paste0(sample_name),
                    paste0("SCOUT_", spn, "_", sample_name, "_mobsterh_st_best_fit.rds")))
    
    obj = readRDS(file.path(path_m,
                            paste0(sample_name),
                            paste0("SCOUT_", spn, "_", sample_name, "_mobsterh_st_best_fit.rds")))
    
    plot(obj)+ggtitle(sample_name)
  })
  
  color_palette = RColorBrewer::brewer.pal(n = max(3,length(unique(table$cluster_id_process))), name = "Dark2") %>% 
    setNames(str_sort(unique(table$cluster_id_process), numeric=T))
  color_palette['Subclonal'] = 'gainsboro'
  
  tool = 'viber'
  table = readRDS(file.path(save_path, "tables_interpreted_new_clusters", paste0(tool, "_", spn, "_", simulation_id, ".rds")))
  
  plots_process = lapply(samples, function(sample_name) {
    
    table_sample = table %>% filter(sample_id == sample_name)
    
    cluster_order = table_sample %>%
      dplyr::count(cluster_id_process) %>%       # Count the number of data points per cluster
      arrange(desc(n)) %>%           # Sort clusters by size (descending order)
      pull(cluster_id_process)
    
    table_sample = table_sample %>%
      mutate(cluster_id_process = factor(cluster_id_process, levels = cluster_order))
    
    # Here plot the marginal of process
    table_sample %>% 
      filter(vaf_process>0) %>% 
      ggplot() +
      geom_histogram(aes(x=vaf_process, fill=cluster_id_process), position="identity", alpha=1, bins=100) +
      ylab("") +
      ggtitle(sample_name)+
      scale_fill_manual(values=color_palette, name="Cluster") +
      theme(legend.position="bottom")#+
      # theme_minimal()
    
  })
  
  ncol_max <- min(3, length(samples))
  # patch_iter = wrap_plots(plots_mobster, ncol = length(samples)) /
  #   wrap_plots(plots_process, ncol = length(samples)) +
  patch_iter = wrap_plots(plots_mobster, ncol = ncol_max) /
    wrap_plots(plots_process, ncol = ncol_max) +
    plot_annotation(
      title = paste0(
        spn, " | ",
        coverage, "x | ",
        purity, "p | ",
        vcf_caller, " | ",
        cna_caller
      )
    )
  
  scatter_process = plot_scatter_process(table, sample_names, color_palette, driver=F)
  # scatter_process
  
  pairs <- combn(samples, 2, simplify = FALSE)
  ncol <- min(3, length(pairs))
  nrow <- ceiling(length(pairs) / ncol)
  
  ggsave(
    filename = file.path(
      save_path,
      "plots/marginals_process_mobster/scatter_multivariate/",
      paste0(spn, "_", simulation_id, ".png")
    ),
    plot = scatter_process,
    width = 3*ncol,
    height = 3*nrow,
    units = "in",
    dpi = 300
  )
  
  png_file <- file.path(
    save_path,
    "plots/marginals_process_mobster/",
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
  file.path(save_path,
            "plots/marginals_process_mobster/"),
  pattern = "\\.png$",
  full.names = TRUE
)

pdf_file <- file.path(save_path, "plots/marginals_process_mobster/", "all_simulations_univariate.pdf")

pdf(pdf_file, width = 10, height = 12)

for (f in png_files) {
  img <- png::readPNG(f)
  grid::grid.newpage()
  grid::grid.raster(img)
}

dev.off()

# Scatter multiavriate
png_files <- list.files(
  file.path(save_path,
            "plots/marginals_process_mobster/scatter_multivariate/"),
  pattern = "\\.png$",
  full.names = TRUE
)

pdf_file <- file.path(save_path, "plots/marginals_process_mobster/", "scatter_multivariate.pdf")

pdf(pdf_file, width = 10, height = 12)

for (f in png_files) {
  img <- png::readPNG(f)
  grid::grid.newpage()
  grid::grid.raster(img)
}

dev.off()



