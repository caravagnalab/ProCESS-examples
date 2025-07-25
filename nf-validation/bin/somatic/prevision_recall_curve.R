
rm(list = ls())
options(bitmapType='cairo')
require(tidyverse)
library(optparse)

source("../../getters/process_getters.R")
source("utils/utils.R")
source("utils/plot_utils.R")

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN03'),
                    make_option(c("--purity"), type = "character", default = '0.6'),
                    make_option(c("--coverage"), type = "character", default = '100'))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
data_dir = '/orfeo/scratch/cdslab/shared/SCOUT/'
spn_id = opt$spn_id
coverage = opt$coverage
purity = opt$purity

data_dir = '/orfeo/scratch/cdslab/shared/SCOUT/'
spn_id = "SPN04"
coverage = "50"
purity = "0.6"

gender = get_process_gender(spn = spn_id)
if (gender=="XX"){
  chromosomes = c(paste0('chr',1:22), 'chrX')
} else {
  chromosomes = c(paste0('chr',1:22), 'chrX', 'chrY')
}

# INPUT PARAMATERS ####

callers = c("mutect2", "strelka", "freebayes")
min_vaf = .02
mut_types = c("INDEL", "SNV")
comb = list(PI = purity, COV = coverage)
samples = get_sample_names(spn = spn_id)

input_dir <-  paste0(data_dir,spn_id,"/validation/somatic/")
outdir <- paste0(data_dir,spn_id,"/validation/somatic/report")
dir.create(outdir, recursive = T, showWarnings = F)

# Preparing report
all_res = dplyr::tibble()
message("Parsing combination: purity=", purity, ", cov=", coverage)
for (caller in callers) {
  message("  Using caller: ", caller)
  for (sample_id in samples) {
    message("    Working with sample : ", sample_id)
    for (mut_type in mut_types) {
      message("      Considering only ", mut_type)
      
      # spn <- gsub(".*SCOUT/(SPN[0-9]+).*", "\\1", path_to_seq)
      # purity <- gsub(".*purity_([0-9.]+).*", "\\1", path_to_seq)
      # coverage <- gsub(".*coverage_([0-9]+).*", "\\1", path_to_seq)
      combination = paste0(coverage, "x_", purity, "p")
      process_folder_path <- file.path(input_dir, spn_id, combination, "process", sample_id, mut_type)
      caller_folder_path = file.path(input_dir, spn_id, combination, caller, sample_id, mut_type)
      
      # Get ground truth
      gt_res = lapply(chromosomes, function(chromosome){
        gt_path = file.path(process_folder_path, paste0(chromosome,".rds"))
        readRDS(gt_path)
      }) %>% do.call("bind_rows", .)
      
      # Get caller res
      caller_res = lapply(chromosomes, function(chromosome) {
        caller_path = file.path(caller_folder_path, paste0(chromosome,".rds"))
        readRDS(caller_path)  
      }) %>% do.call("bind_rows", .) %>% 
        dplyr::filter(!is.na(VAF))
      
      # Run the analysis with custom parameters
      results <- analyze_vaf_performance(seq_res_long = gt_res, caller_res = caller_res,
                                         only_pass = TRUE, 
                                         vaf_tolerance_pct = 5,  # Allow 10% relative difference
                                         min_vaf_threshold = 0.02) # Minimum VAF of 2%
      all_res = dplyr::bind_rows(
        all_res,
        cbind(results$performance_table, dplyr::tibble(caller=caller, sample_id=sample_id, mut_type=mut_type))
      )      
      
    }
  }
}

all_res = all_res %>% dplyr::mutate(caller = factor(caller, levels = c("freebayes", "strelka", "mutect2")))

all_res %>% 
  ggplot(mapping = aes(x = VAF_bin, y = sensitivity_pct, col = caller)) +
  geom_boxplot() +
  #geom_jitter() +
  labs(
    title = "Variant Caller Sensitivity by VAF Range",
    x = "VAF Range (Truth)",
    y = "Sensitivity (%)",
    col = "Variant caller"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("strelka" = "steelblue", "mutect2" = "coral", "freebayes" = "#8FBC8B")) +
  facet_wrap(~mut_type, nrow = 2) +
  scale_y_continuous(limits = c(0,100))

all_res %>% 
  ggplot(mapping = aes(x = VAF_bin, y = sensitivity_pct, col = caller)) +
  geom_boxplot() +
  #geom_jitter() +
  labs(
    title = "Variant Caller Sensitivity by VAF Range",
    x = "VAF Range (Truth)",
    y = "Sensitivity (%)",
    col = "Variant caller"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("strelka" = "steelblue", "mutect2" = "coral", "freebayes" = "#8FBC8B")) +
  facet_wrap(~mut_type, nrow = 2) +
  scale_y_continuous(limits = c(0,100))

all_res %>%
  mutate(VAF_bin_num = as.numeric(factor(VAF_bin))) %>%  # Convert VAF_bin to numeric for plotting
  ggplot(aes(x = VAF_bin_num, y = sensitivity_pct, col = caller)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "loess") +
  labs(
    title = "Variant Caller Sensitivity by VAF Range",
    x = "VAF Range (Truth)",
    y = "Sensitivity (%)",
    col = "Variant caller"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("strelka" = "steelblue", "mutect2" = "coral", "freebayes" = "#8FBC8B")) +
  facet_wrap(~mut_type, nrow = 2) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(
    breaks = unique(as.numeric(factor(all_res$VAF_bin))),
    labels = levels(factor(all_res$VAF_bin))
  )

all_res %>%
  mutate(VAF_bin_num = as.numeric(factor(VAF_bin))) %>%  # Convert VAF_bin to numeric for plotting
  dplyr::group_by(VAF_bin_num, caller, VAF_bin, mut_type) %>% 
  dplyr::summarise(y = mean(sensitivity_pct), s = sd(sensitivity_pct)) %>% 
  ggplot(aes(x = VAF_bin_num, y = y, ymin=y-s, ymax=y+s, col = caller)) +
  geom_pointrange() +
  geom_line() +
  labs(
    title = "Variant Caller Sensitivity by VAF Range",
    x = "VAF Range (Truth)",
    y = "Sensitivity (%)",
    col = "Variant caller"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("strelka" = "steelblue", "mutect2" = "coral", "freebayes" = "#8FBC8B")) +
  facet_wrap(~mut_type, nrow = 2) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(
    breaks = unique(as.numeric(factor(all_res$VAF_bin))),
    labels = levels(factor(all_res$VAF_bin))
  )


all_res %>% 
  ggplot(mapping = aes(x = VAF_bin, y = precision_pct, col = caller)) +
  geom_boxplot() +
  #geom_jitter() +
  labs(
    title = "Variant Caller Precision by VAF Range",
    x = "VAF Range (Truth)",
    y = "Precision (%)",
    col = "Variant caller"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("strelka" = "steelblue", "mutect2" = "coral", "freebayes" = "#8FBC8B")) +
  facet_wrap(~mut_type, nrow = 2) +
  scale_y_continuous(limits = c(0,100))

