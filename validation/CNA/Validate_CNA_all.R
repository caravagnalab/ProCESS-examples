options(bitmapType='cairo')
library(dplyr)
library(ProCESS)
library(optparse)
library(tidyr)
library(ggplot2)
library(patchwork)
source("../../getters/sarek_getters.R")
source("../../getters/process_getters.R")
source("utils.R")

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN01'),
                    make_option(c("--coverages"), type = "character", default = '50'),
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

df_bp = lapply(1:nrow(params_grid), function(i) {
  coverage = params_grid[i,]$coverage
  purity = params_grid[i,]$purity
  
  combination = paste0(coverage, "x_", purity, "p")
  
  results_folder_path = file.path(input_dir, spn_id, combination)
  samples_name = get_sample_names(spn_id)
  
  tmp_df = lapply(samples_name, FUN = function(sample){
    file_name = file.path(results_folder_path, sample, "metrics_bp.rds")
    if (file.exists(file_name)){
      metrics = readRDS(file_name) %>% mutate(true_purity = as.numeric(true_purity),
                                              coverage = as.numeric(coverage))
    }
  }) %>% bind_rows()
}) %>% bind_rows()

df_metric <- df_metric %>% mutate(fga = round(fga),
                                  fgs = round(fgs))

df_bp <- df_bp %>% mutate(fga = round(fga),
                                  fgs = round(fgs))

plt_data <- df_metric %>% 
  select(sample, fga,fgs) %>% 
  distinct() %>% 
  pivot_longer(cols = c(fga,fgs)) %>% 
  ggplot() + 
  geom_col(aes(x = sample, y = value,fill=name ),position = position_dodge()) +
  scale_fill_manual('', values = c('indianred', 'orange'))  +
  ylab('% of genome') +
  theme_minimal() 

plt_ploidy <- df_metric %>% 
  select(sample, true_ploidy) %>% 
  distinct() %>% 
  mutate(true_ploidy = round(true_ploidy,1)) %>% 
  ggplot(aes(x = sample, y = true_ploidy )) + 
  geom_col() +
  theme_minimal()

plt_bp <- df_bp %>% 
  filter(chr == 'genome') %>% 
  pivot_longer(cols = c(precision, recall, f1)) %>% 
  ggplot() +
  geom_boxplot(aes(x = as.factor(true_purity), y = value, col = tool)) +
  scale_color_manual(values = color_caller) +
  facet_grid(name~coverage)  +
  xlab('purity') + 
  theme_bw() 

plt <- plt_ploidy + 
  plt_data + 
  
  df_metric %>%
  mutate(delta_purity = as.numeric(true_purity) - as.numeric(purity)) %>%
  ggplot(aes(x = sample, y = delta_purity, fill = tool)) +
  geom_col(position = position_dodge(width = 0.35), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = color_caller) +
  ylab("true_purity - inferred_purity") +
  facet_grid(as.numeric(coverage) ~ as.numeric(true_purity)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 

  df_metric %>%
  mutate(delta_ploidy = as.numeric(true_ploidy) - as.numeric(ploidy)) %>%
  ggplot(aes(x = sample, y = delta_ploidy, fill = tool)) +
  geom_col(position = position_dodge(width = 0.35), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = color_caller) +
  ylab("true_ploidy - inferred_ploidy") +
  facet_grid(as.numeric(coverage) ~ as.numeric(true_purity)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  
df_metric %>% 
  ggplot() + 
  geom_col(aes(x = sample, y = correctness_clonal*100, fill = tool), position = position_dodge()) + 
  scale_fill_manual(values = color_caller) +
  ylab('% CN correctness') +
  theme_bw() +
  facet_grid(as.numeric(coverage)  ~ as.numeric(true_purity))  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  plt_bp + 
  
  plot_layout(design = 'AB\nCC\nDD\nEE\nFF') +
  plot_annotation(title = spn_id) + plot_layout(guides = 'collect')

ggsave(filename = file.path(data_dir, spn_id, 'validation/cna/report', paste0(spn_id, '.pdf')), plot = plt, dpi = 400, height = 12, width = 12, units = 'in')
plt
