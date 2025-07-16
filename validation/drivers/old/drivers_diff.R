rm(list=ls())
setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/rRACES/ProCESS-examples')
# drivers validation
library(tidyverse)
library(ProCESS)
library(ComplexHeatmap)
library(UpSetR)


# getters
source('getters/process_getters.R')
source('getters/sarek_getters.R')
source('validation/drivers/drivers_validation_functions.R')

# get samples and phylogenetic forest
samples = get_sample_names(spn = 'SPN03')
phylo_forest = get_phylo_forest('SPN03', base_path = '/orfeo/cephfs/scratch/cdslab/shared/SCOUT')
phylo_forest = load_phylogenetic_forest(phylo_forest)

# results from tumourevo --> annotated drivers
drivers_annotated = get_drivers_results(spn = 'SPN03', 
                                        samples = samples, 
                                        purity = '0.9', 
                                        coverage = '50x', 
                                        callers = 'mutect2_ascat', 
                                        cohort = 'SCOUT')

# laod results of process sequencing
process_seq_res = get_process_seq(spn = 'SPN03', coverage = '50x', purity = 0.9)

# filter the mutation tables to get only mutation ids that are labelled as driver in one of the two






  
dd = all_drivers %>% 
  dplyr::select(sample, driver_label, driver_class) %>% 
  tidyr::pivot_wider(names_from = driver_label, 
                     values_from = driver_class) %>% 
  tibble::column_to_rownames('sample')

dr_ht = Heatmap(as.matrix(dd), 
        col = colors, 
        rect_gp = gpar(col = "white", lwd = 3), 
        name = 'Drivers classes')
draw(dr_ht, heatmap_legend_side = "bottom")
  
all_driver_split = all_drivers %>% 
  dplyr::filter(driver_class != 'Not a driver') %>% 
  dplyr::select(mut_id,sample, is_driver_process, is_driver_tumourevo) %>% 
  tidyr::pivot_longer(names_to = 'origin', cols = c(is_driver_tumourevo, is_driver_process)) %>% 
  group_by(sample) %>% 
  group_split() 

names(all_driver_split) = lapply(all_driver_split, function(x) {
  x$sample %>% unique
}) %>% unlist

all_driver_split = lapply(all_driver_split, function(x) {
  x = x %>% 
    group_by(origin) %>% 
    group_split()
  names(x) = lapply(x, function(s) {
    s$origin %>% unique
  }) %>% unlist
  x = lapply(x, function(s) {
    s %>% 
      dplyr::filter(value == TRUE) %>% 
      dplyr::pull(mut_id)
  })
  return(x)
})
