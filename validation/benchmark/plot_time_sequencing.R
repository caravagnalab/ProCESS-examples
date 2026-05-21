library(dplyr)
library(hms)
library(ggplot2)
library(patchwork)
library(lubridate)
library(optparse)
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/benchmark/utils.R')
spns <- paste0('SPN0', 1:7)

all_time_sequencing <- tibble()
all_time_merging <- tibble()
all_time_samtools <- tibble()

sample_table <- tibble(sample = spns, N = c(3,2,4,2,3,5,5))
for (spn in spns){
  base <- paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/',spn,'/sequencing')
  if (dir.exists(base)){
    print(spn)
    tum <- parse_rds(base, type = 'tumour',merging = FALSE) %>%
      group_by(coverage, purity, tumour, type, chunk) %>%
      summarise(elapsed_time_mins = sum(elapsed_time_mins),
                cpu_time_secs = sum(cpu_time_secs),
                memory_used_MB = sum(memory_used_MB)) %>%
      mutate(SPN = spn,  N = sample_table %>% filter(sample == spn) %>% pull(N))
    tum_merging <- parse_rds(base, type = 'tumour',merging = TRUE) %>%
      mutate(SPN = spn,  N = sample_table %>% filter(sample == spn) %>% pull(N))

    nor <- parse_rds(base, 'normal', merging = FALSE)  %>%
      group_by(coverage, purity, tumour, type, chunk) %>%
      summarise(elapsed_time_mins = sum(elapsed_time_mins),
                cpu_time_secs = sum(cpu_time_secs),
                memory_used_MB = sum(memory_used_MB)) %>%
      mutate(SPN = spn,  N = sample_table %>% filter(sample == spn) %>% pull(N))
    nor_merging <- parse_rds(base, type = 'normal',merging = TRUE) %>%
      mutate(SPN = spn,  N = sample_table %>% filter(sample == spn) %>% pull(N))

    time_sequencing <- bind_rows(tum,nor)
    time_merging <- bind_rows(tum_merging,nor_merging)
    all_time_sequencing <- bind_rows(all_time_sequencing, time_sequencing)
    all_time_merging <- bind_rows(all_time_merging, time_merging)

    time_t <- parse_out(base, type = 'tumour') %>%
      mutate(SPN = spn, N = sample_table %>% filter(sample == spn) %>% pull(N))
    time_n <- parse_out(base, type = 'normal') %>%
      mutate(SPN = spn,  N = sample_table %>% filter(sample == spn) %>% pull(N))
    time <- bind_rows(time_t, time_n)
    all_time_samtools <- bind_rows(all_time_samtools, time)
  }
}

full_table <- bind_rows(all_time_merging %>% mutate(step = 'merge_rds') %>% filter(type == 'tumour') %>% select(cpu_time_secs, memory_used_MB, SPN,N, step),
                        all_time_sequencing %>% mutate(step = 'sequencing') %>% filter(type == 'tumour')  %>% select(cpu_time_secs, memory_used_MB, SPN,N, step),
                        all_time_samtools %>% filter(type == 'tumour') %>% mutate(cpu_time_secs = usr_time_sec,
                                                                                  memory_used_MB = memory_MB) %>% select(cpu_time_secs, memory_used_MB, SPN,N, step))
saveRDS(object = full_table, file = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_2/time_sequencing.rds')
# 
# memory <- full_table %>% 
#   ggplot(aes(x = SPN, y = memory_used_MB)) +
#   ylab('memory (MB)') + 
#   xlab('') +
#   geom_boxplot(aes(col = as.factor(N))) +
#   geom_jitter(aes(col = as.factor(N)), height = 0, width = 0.1, size = 1) +
#   scale_color_manual('N samples', values = c('#7CCAD5', '#A0A6BE', '#C481A7', '#454995')) + 
#   facet_grid(step~., scales = 'free_y') +
#   theme_custom() +
#   plot_annotation(caption = 'Memory for each lot (N = 40, coverage = 5X, purities = 0.3,0.6,0.9)') & theme(legend.position = 'bottom')
# 
# time <- full_table %>% 
#   ggplot(aes(x =SPN, y = hms::as_hms(cpu_time_secs))) +
#   ylab('time (H:M:S)') +
#   xlab('') + 
#   geom_boxplot(aes(col = as.factor(N))) +
#   geom_jitter(aes(col = as.factor(N)), height = 0, width = 0.1, size = 1) +
#   theme_custom() +
#   scale_color_manual('N samples', values = c('#7CCAD5', '#A0A6BE', '#C481A7', '#454995'))  +
#   facet_grid(step~., scales = 'free_y')  +
#   plot_layout(guides = 'collect') + 
#   plot_annotation(caption = 'Time for each lot (N = 40, coverage = 5X, purities = 0.3,0.6,0.9)') & theme(legend.position = 'bottom')
# 
# ggsave(filename = 'plot_time_sequencing.png', plot = time, width = 8, height = 8, units = 'in', dpi = 600)
# ggsave(filename = 'plot_mem_sequencing.png', plot = memory, width = 8, height = 8, units = 'in', dpi = 600)
# 
# 
# df_summary <- full_table %>% 
#   mutate(x =  hms::as_hms(cpu_time_secs)) %>% 
#   select(-memory_used_MB) %>% 
#   group_by(SPN, N, step) %>% 
#   summarise(
#     mean = mean(cpu_time_secs, na.rm = TRUE),
#     se   = sd(cpu_time_secs, na.rm = TRUE) / sqrt(n()),
#     .groups = "drop"
#   ) %>% 
#   mutate(lower = mean - se, upper = mean + se )
# 
# 
# df_summary <- df_summary %>% mutate(step = case_when(
#         step == 'fastq' ~ '4. Convert to FASTQ',
#         step == 'merge_rds' ~ '2. Merge RDS',
#         step == 'samtools_merge' ~ '2. Merge SAM',
#         step == 'samtools_split' ~ '3. Split SAM',
#         step == 'sequencing' ~ '1. Sequencing ProCESS'))
# 
# mean_plt <- ggplot(df_summary, aes(x = SPN, y =  hms::as_hms(mean), fill = as.factor(N))) +
#   geom_col(position = position_dodge(width = 0.9)) +
#   geom_errorbar(aes(ymin = hms::as_hms(lower), ymax = hms::as_hms(upper)),
#                 width = 0.2,
#                 position = position_dodge(width = 0.9)) +
#   scale_fill_manual('N samples', values = c('#7CCAD5', '#A0A6BE', '#C481A7', '#454995'))  +
#   theme_minimal(base_size = 13) +
#   facet_grid(step~., scales = 'free_y')  +
#   ylab('time (H:M:S)') +
#   xlab('SPN') +
#   plot_annotation(caption = 'Mean time over each lot (N = 40, coverage = 5X, purities = 0.3,0.6,0.9)') & theme(legend.position = 'bottom')
# ggsave(filename = 'plot_time_sequencing_v2.png', plot = mean_plt, width = 8, height = 10, units = 'in', dpi = 600)
# 
# 
