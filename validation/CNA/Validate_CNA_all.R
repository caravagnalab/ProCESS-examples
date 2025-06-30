options(bitmapType='cairo')
library(dplyr)
library(ProCESS)
library(optparse)
library(tidyr)
library(ggplot2)
library(future.apply)
library(progressr)
source("../../getters/sarek_getters.R")
source("../../getters/process_getters.R")
source("utils.R")

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN06'),
                    make_option(c("--coverages"), type = "character", default = '50, 100'),
                    make_option(c("--purities"), default = '0.6, 0.9, 0.3')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
data_dir = '/orfeo/scratch/cdslab/shared/SCOUT/'
spn_id = opt$spn_id

cleaned <- gsub('^"|"$', '', opt$purities)
purity_list <- strsplit(cleaned, ",")[[1]]
PURITY <- trimws(purity_list)
cleaned <- gsub('^"|"$', '', opt$coverages)
coverage_list <- strsplit(cleaned, ",")[[1]]
COVERAGES <- trimws(coverage_list)
print(COVERAGES)
print(PURITY)

params_grid = expand.grid(COVERAGES, PURITY)
colnames(params_grid) = c("coverage", "purity")

input_dir <-  paste0(data_dir,spn_id,"/validation/cna/")
i = 1

df_metric = lapply(1:nrow(params_grid), function(i) {
  coverage = params_grid[i,]$coverage
  purity = params_grid[i,]$purity
  
  combination = paste0(coverage, "x_", purity, "p")
  
  results_folder_path = file.path(input_dir, spn_id, combination)
  samples_name = get_sample_names(spn_id)
  
  tmp_df = lapply(samples_name, FUN = function(sample){
    file_name = file.path(results_folder_path, sample, "metrics.rds")
    if (file.exists(file_name)){
      metrics = readRDS(file_name) %>% mutate(true_purity = as.numeric(true_purity),
                                              coverage = as.numeric(coverage))
    }
  }) %>% bind_rows()
}) %>% bind_rows()


plt <- df_metric %>% 
  ggplot() + 
  geom_point(aes(x = sample, y = as.numeric(true_purity) - as.numeric(purity), col = tool, size = fgs, alpha = 1-fga)) + # size = fgs
  geom_hline(aes(yintercept = 0)) +
  ylab('true_purity - inferred_purity') +
  theme_bw() +
  facet_grid(as.numeric(coverage) ~ as.numeric(true_purity)) +

df_metric %>% 
  ggplot() + 
  geom_point(aes(x = sample, y = as.numeric(true_ploidy) - as.numeric(ploidy), col = tool,  size = fgs, alpha = 1-fga)) + 
  geom_hline(aes(yintercept = 0)) +
  ylab('true_ploidy - inferred_ploidy') +
  theme_bw() +
  facet_grid(as.numeric(coverage)  ~ as.numeric(true_purity))  +
  
df_metric %>% 
  ggplot() + 
  geom_point(aes(x = sample, y = correctness_clonal, col = tool,  size = fgs, alpha = 1-fga)) + 
  ylab('% correctness') +
  theme_bw() +
  facet_grid(as.numeric(coverage)  ~ as.numeric(true_purity))  +
  
df_metric %>% 
  ggplot() + 
  geom_point(aes(x = sample, y = bp_distance, col = tool,  size = fgs, alpha = 1-fga)) + 
  ylab('breakpoint distance') +
  theme_bw() +
  facet_grid(as.numeric(coverage)  ~ as.numeric(true_purity)) + 
  plot_layout(nrow = 4) +
  plot_annotation(title = spn_id) + plot_layout(guides = 'collect')


plt



