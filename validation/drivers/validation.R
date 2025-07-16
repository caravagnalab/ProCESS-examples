.libPaths(new = '/u/area/vgazziero/R/rstudio_4_4')

# drivers validation
library(tidyverse)
library(ProCESS)

if (!require('devtools', quietly = T)) 
  install.packages('devtools', repos = 'https://cran.stat.unipd.it/')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = 'https://cran.stat.unipd.it/')
if (!require("ComplexHeatmap", quietly = TRUE)) 
  BiocManager::install("ComplexHeatmap")
if (!require('circlize', quietly = TRUE)) 
  install.packages('circlize', repos = 'https://cran.stat.unipd.it/')
if (!require("viridis", quietly = TRUE)) 
  install.packages("viridis", repos = 'https://cran.stat.unipd.it/')
if (!require('awtools', quietly = TRUE)) 
  devtools::install_github('awhstin/awtools')
if (!require(optparse, quietly = TRUE)) 
  install.packages('optparse', repos = 'https://cran.stat.unipd.it/')

library(ComplexHeatmap)
library(circlize)
library(viridis)
library(awtools)
library(cowplot)
library(grid)
library(optparse)

# command line arguments 

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN07'),
                    make_option(c("--purity"), type = "character", default = '0.9'),
                    make_option(c("--coverage"), type = "character", default = '50'), 
                    make_option(c("--caller"), type = "character", default = 'mutect2_ascat'), 
                    make_option(c("--process_example"), type = "character", default = "/orfeo/cephfs/scratch/cdslab/vgazziero/rRACES/ProCESS-examples"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# define the parameters
data_dir = '/orfeo/scratch/cdslab/shared/SCOUT'
spn_id = opt$spn_id
coverage = opt$coverage
purity = opt$purity
caller = opt$caller
workdir = opt$process_example
real_caller = str_split_1(caller, '_')[1]

# check what you have to do
res_path = paste0(data_dir,'/', spn_id, '/validation/drivers/', coverage, 'x_', purity, 'p/', real_caller)
dir.create(res_path, showWarnings = FALSE, recursive = TRUE)
files_done = list.files(res_path)
expected_files = c('drivers_validation_report.pdf', 'drivers_validation_difference.rds', 'real_drivers_vaf_comparison.rds')
to_do = setdiff(expected_files, files_done)

setwd(workdir)

# getters
source('getters/process_getters.R')
source('getters/sarek_getters.R')
source('validation/drivers/drivers_validation_functions.R')

# loading data if there is something to do!

if (length(to_do) > 0) {
  
  print('Loading data')
  
  samples = get_sample_names(spn = spn_id)
  
  # load data processed by tumourevo (and if they do not exist, kill the job)
  tumourevo_mutations = get_drivers_results(spn = spn_id, 
                                            samples = samples, 
                                            purity = purity, 
                                            coverage = coverage, 
                                            callers = caller, 
                                            cohort = 'SCOUT', 
                                            path = data_dir)
  
  # phylo_forest = get_phylo_forest(spn_id, base_path = data_dir)
  # phylo_forest = load_phylogenetic_forest(phylo_forest)
  
  # process_drivers = phylo_forest$get_driver_mutations()
  
  # load process data and select its drivers
  process_seq_res = get_process_seq(spn = spn_id, 
                                    purity = purity, 
                                    coverage = coverage,
                                    path = data_dir)
  
  
  # if the heatmap has not been done, create it
  if ('real_drivers_vaf_comparison.rds' %in% to_do) {
    # get the mutation id of the process drivers
    process_drivers_ids = get_process_drivers_ids(process_seq_res)
    
    # filter the sequencing results to keep only driver mutations
    process_drivers = get_process_drivers(process_seq_res)
    
    true_drivers_table = create_true_drivers_table(process_drivers,
                                                   tumourevo_mutations, 
                                                   process_drivers_ids)
    
    # plot the heatmap
    ht = plot_drivers_heatmap(true_drivers_table)
    saveRDS(ht, paste0(res_path, '/real_drivers_vaf_comparison.rds'))
  } else {
    print('Heatmap already saved')
  }
  
  # check how much is shared the driver labelling among the samples 
  
  # get tumourevo drivers ids
  if ('drivers_validation_difference.rds' %in% to_do){
    process_drivers_ids = get_process_drivers_ids(process_seq_res)
    tumourevo_drivers_ids = get_tumourevo_drivers_ids(tumourevo_mutations)
    all_drivers = c(process_drivers_ids, tumourevo_drivers_ids) %>% unique
    
    process_all_drs = get_all_drivers_process(process_seq_res, all_drivers)
    
    tumourevo_all_drs = get_all_drivers_tumourevo(tumourevo_mutations, all_drivers)
    
    all_drivers_table = merge_drivers(process_all_drs, tumourevo_all_drs)
    plt = plot_drivers(all_drivers_table, colors = colors)
    saveRDS(plt, paste0(res_path, '/drivers_validation_difference.rds'))
  } else {
    print('Plot already saved')
  }
  
  if ('drivers_validation_report.pdf' %in% to_do){
    
    ht = readRDS(paste0(res_path, '/real_drivers_vaf_comparison.rds'))
    plt = readRDS(paste0(res_path, '/drivers_validation_difference.rds'))
    # arrange everything in a single page 
    
    ht_grob <- grid.grabExpr(draw(ht))
    # gg_grob <- ggplotGrob(plt)
    
    
    pdf(paste0(res_path, '/drivers_validation_report.pdf'), width = 10, height = 20)
    
    print(patchwork::wrap_plots(list(ht_grob, plt), ncol = 1) + 
      patchwork::plot_annotation(title = paste0(
        "SPN: ", spn_id, '\n',
        "Purity: ", purity,'\n',
        "Coverage: ", coverage, 'x', 
        "Caller: ", real_caller
      )))
    dev.off()
  } else {
    print('Report already saved!')
  }
} else {
  
  print('Everything already done! Skipping validation for this SPN with this combination of caller, coverage and purity')
  
}
# counting how many false positive and false negative we have









# #see matches between the simulation and tumourevo
# all_tumourevo_drivers = lapply(names(drivers_annotated), function(x) {
#   drivers_annotated[[x]][[1]]$mutations %>% 
#     dplyr::mutate(mut_id = paste0(chr, ':', from, '_', ref, '>', alt)) %>% 
#     dplyr::filter(is_driver | mut_id %in% process_drivers_ids) %>% 
#     dplyr::mutate(sample = x) %>% 
#     dplyr::select(mut_id, driver_label, sample, is_driver) %>% 
#     dplyr::rename(is_driver_tumourevo = is_driver)
# })
# 
# names(all_tumourevo_drivers) = lapply(all_tumourevo_drivers, function(x) {
#   x$sample %>% unique
# }) %>% unlist
# 
# process_drivers = process_seq_res %>% 
#   dplyr::filter(classes == 'driver') %>% 
#   dplyr::mutate(chr = paste0('chr', chr)) %>% 
#   dplyr::mutate(mut_id = paste0(chr, ':', chr_pos, '_', ref, '>', alt))
# 
# 
# 
# 
# driver_comparison = lapply(names(drivers_tumourevo), function(x) {
#   
#   te = drivers_tumourevo[[x]]
#   te = te %>% 
#     dplyr::rename(tumourevo_NV = NV) %>% 
#     dplyr::rename(tumourevo_DP = DP) %>% 
#     dplyr::rename(tumourevo_VAF = VAF)
#   
#   pr = process_drivers %>% 
#     dplyr::select(chr, chr_pos, ref, alt, causes, starts_with(x)) %>% 
#     dplyr::mutate(chr = paste0('chr', chr))
#   
#   colnames(pr) = c('classes', 'chr', 'from', 'ref', 'alt', 'causes', 'process_NV', 'process_DP', 'process_VAF')
#   
#   full_join(te, pr)
#   
# })
# 
# names(driver_comparison) = names(drivers_tumourevo)

