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

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN04'),
                    make_option(c("--coverages"), type = "character", default = '50, 100'),
                    make_option(c("--purities"), default = '0.6, 0.9, 0.3')
)

colors = c('ascat' = 'coral2', 'sequenza' = 'darkslategray4', 'cnvkit' = 'maroon')

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

df_metric <- df_metric %>% mutate(fga = round(fga),
                                  fgs = round(fgs))

plt <- df_metric %>% 
  ggplot() + 
  geom_point(aes(x = sample, y = as.numeric(true_purity) - as.numeric(purity), col = tool)) +
  scale_color_manual(values = colors) +
  geom_hline(aes(yintercept = 0)) +
  ylab('true_purity - inferred_purity') +
  theme_bw() +
  facet_grid(as.numeric(coverage) ~ as.numeric(true_purity)) +

df_metric %>% 
  ggplot() + 
  geom_point(aes(x = sample, y = as.numeric(true_ploidy) - as.numeric(ploidy), col = tool)) + 
  scale_color_manual(values = colors) +
  geom_hline(aes(yintercept = 0)) +
  ylab('true_ploidy - inferred_ploidy') +
  theme_bw() +
  facet_grid(as.numeric(coverage)  ~ as.numeric(true_purity))  +
  
df_metric %>% 
  ggplot() + 
  geom_point(aes(x = sample, y = correctness_clonal, col = tool)) + 
  scale_color_manual(values = colors) +
  ylab('% correctness') +
  theme_bw() +
  facet_grid(as.numeric(coverage)  ~ as.numeric(true_purity))  +
  
df_metric %>% 
  ggplot() + 
  geom_point(aes(x = sample, y = bp_distance, col = tool)) + 
  scale_color_manual(values = colors) +
  ylab('breakpoint distance') +
  theme_bw() +
  facet_grid(as.numeric(coverage)  ~ as.numeric(true_purity)) + 
  plot_layout(nrow = 4) +
  plot_annotation(title = spn_id) + plot_layout(guides = 'collect')

ggsave(filename = file.path(data_dir, spn_id, 'validation/cna/report', paste0(spn_id, '.pdf')), plot = plt, dpi = 400, height = 12, width = 12, units = 'in')


p1 = df_metric %>%
  filter(tool != 'cnvkit') %>%
  group_by(sample) %>%
  dplyr::mutate(sample_fga = mean(fga, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::mutate(
    purity_err =  true_purity - purity,
    sample_fga_rank = rank(sample_fga)
  ) %>%
  ggplot(mapping = aes(x = true_purity, y = purity_err, col = sample_fga_rank, group = sample, shape = sample)) +
  geom_point(size = 3) +
  geom_line() +
  facet_grid(coverage~tool) +
  scale_fill_viridis_c(name = "FGA Rank") +
  scale_color_viridis_c(name = "FGA Rank") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "True Purity", y = "Delta Purity (True - Predicted)", fill = "FGA", title = 'FGA')

p2 = df_metric %>%
  filter(tool != 'cnvkit') %>%
  group_by(sample) %>%
  dplyr::mutate(sample_fgs = mean(fgs, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::mutate(
    purity_err =  true_purity - purity,
    sample_fgs_rank = rank(sample_fgs)
  ) %>%
  ggplot(mapping = aes(x = true_purity, y = purity_err, col = sample_fgs_rank, group = sample, shape = sample)) +
  geom_point(size = 3) +
  geom_line() +
  facet_grid(coverage~tool) +
  scale_fill_viridis_c(name = "FGS Rank") +
  scale_color_viridis_c(name = "FGS Rank") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "True Purity", y = "Delta Purity (True - Predicted)", fill = "FGS", title ='FGS')

p3 = df_metric %>%
  filter(tool != 'cnvkit') %>%
  group_by(sample) %>%
  dplyr::mutate(sample_fga = mean(fga, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::mutate(
    ploidy_err =  true_ploidy - ploidy,
    sample_fga_rank = rank(sample_fga)
  ) %>%
  ggplot(mapping = aes(x = true_purity, y = ploidy_err, col = sample_fga_rank, group = sample, shape = sample)) +
  geom_point(size = 3) +
  geom_line() +
  facet_grid(coverage~tool) +
  scale_fill_viridis_c(name = "FGA Rank") +
  scale_color_viridis_c(name = "FGA Rank") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "True Purity", y = "Delta Ploidy (True - Predicted)", fill = "FGA")

p4 = df_metric %>%
  filter(tool != 'cnvkit') %>%
  group_by(sample) %>%
  dplyr::mutate(sample_fgs = mean(fgs, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::mutate(
    ploidy_err =  true_ploidy - ploidy,
    sample_fgs_rank = rank(sample_fgs)
  ) %>%
  ggplot(mapping = aes(x = true_purity, y = ploidy_err, col = sample_fgs_rank, group = sample, shape = sample)) +
  geom_point(size = 3) +
  geom_line() +
  facet_grid(coverage~tool) +
  scale_fill_viridis_c(name = "FGS Rank") +
  scale_color_viridis_c(name = "FGS Rank") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "True Purity", y = "Delta Ploidy (True - Predicted)", fill = "FGS")

p5 = df_metric %>%
  group_by(sample) %>%
  dplyr::mutate(sample_fga = mean(fga, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::mutate(
    ploidy_err =  true_ploidy - ploidy,
    sample_fga_rank = rank(sample_fga)
  ) %>%
  ggplot(mapping = aes(x = true_purity, y = correctness_clonal, col = sample_fga_rank, group = sample, shape = sample)) +
  geom_point(size = 3) +
  geom_line() +
  facet_grid(coverage~tool) +
  scale_fill_viridis_c(name = "FGA Rank") +
  scale_color_viridis_c(name = "FGA Rank") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "True Purity", y = "% of correct CN", fill = "FGA")

p6 = df_metric %>%
  group_by(sample) %>%
  dplyr::mutate(sample_fgs = mean(fgs, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::mutate(
    ploidy_err =  true_ploidy - ploidy,
    sample_fgs_rank = rank(sample_fgs)
  ) %>%
  ggplot(mapping = aes(x = true_purity, y = correctness_clonal, col = sample_fgs_rank, group = sample, shape = sample)) +
  geom_point(size = 3) +
  geom_line() +
  facet_grid(coverage~tool) +
  scale_fill_viridis_c(name = "FGS Rank") +
  scale_color_viridis_c(name = "FGS Rank") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "True Purity", y = "% of correct CN", fill = "FGS")


pp <- p1 + p2 + p3 + p4 + p5 +theme(legend.position = 'none') + p6 +theme(legend.position = 'none') + plot_layout(ncol=2, guides = 'collect')
ggsave(filename = file.path(data_dir, spn_id, 'validation/cna/report', paste0(spn_id, '_v2.pdf')), plot = pp, dpi = 400, height = 12, width = 12, units = 'in')
