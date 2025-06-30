library(dplyr)
library(hms)
library(ggplot2)
library(patchwork)
library(lubridate)
library(optparse)
source('utils.R')
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
saveRDS(object = full_table, file = 'table_data.rds')
full_table <- readRDS('table_data.rds')

generation <- full_table %>% 
  filter(SPN != 'SPN07') %>% 
  ggplot(aes(x = SPN, y = memory_used_MB)) +
  ylab('memory (MB)') + 
  xlab('') +
  geom_boxplot(aes(col = as.factor(N))) +
  geom_jitter(aes(col = as.factor(N)), height = 0, width = 0.1, size = 1) +
  scale_color_manual('N samples', values = c('#7CCAD5', '#A0A6BE', '#C481A7', '#454995')) + 
  facet_grid(step~., scales = 'free_y') +
  theme_custom() + 
  full_table %>% 
  filter(SPN != 'SPN07') %>% 
  ggplot(aes(x =SPN, y = hms::as_hms(cpu_time_secs))) +
  ylab('time (H:M:S)') +
  xlab('') + 
  geom_boxplot(aes(col = as.factor(N))) +
  geom_jitter(aes(col = as.factor(N)), height = 0, width = 0.1, size = 1) +
  theme_custom()+ 
  scale_color_manual('N samples', values = c('#7CCAD5', '#A0A6BE', '#C481A7', '#454995'))  +
  facet_grid(step~., scales = 'free_y')  +
  plot_layout(guides = 'collect') + plot_annotation(caption = 'Time and memory for each lot (N = 40, coverage = 5X)')& theme(legend.position = 'bottom')

ggsave(filename = 'plot_benchmark.png', plot = generation, width = 12, height = 8, units = 'in', dpi = 600)


memory <- tibble(sample = spns,
                 N = c(3,2,4,2,3,5,5),
                 sample_forest = c('1.2', '0.8', '1.9', '0.7', NA, '2.4', '2.1'), #M
                 phylo_forest = c('2.1', '6', '1.8', '0.8', NA, '2.6', '2.8'), #G
                 sam = rep(NA, 7),
                 bam = c('0.9', '0.6', '1.3', '0.6', '1.3', '1.7', '1.7'), #T
                 fastq = c('0.7', '0.5', '0.9', '0.4', '0.7', '1.2', '1.2') #T
                 )

mem <- memory %>% 
  filter(sample != 'SPN05') %>% 
  ggplot() +
  geom_col(aes(x = sample, y = sample_forest, fill = as.factor(N))) +
  scale_fill_manual('N samples', values = c('#7CCAD5', '#A0A6BE', '#C481A7', '#454995'))  +
  ylab('MB') +
  ggtitle('Sample forest') +
  theme_minimal()  +
  
memory %>% 
  filter(sample != 'SPN05') %>% 
  ggplot() +
  geom_col(aes(x = sample, y = phylo_forest, fill = as.factor(N))) +
  scale_fill_manual('N samples', values = c('#7CCAD5', '#A0A6BE', '#C481A7', '#454995'))  +
  ylab('GB') +
  ggtitle('Phylogenetic forest') +
  theme_minimal()  +
  
memory %>% 
  ggplot() +
  geom_col(aes(x = sample, y = fastq, fill = as.factor(N))) +
  scale_fill_manual('N samples', values = c('#7CCAD5', '#A0A6BE', '#C481A7', '#454995'))  +
  ylab('TB') +
  ggtitle('fastq files') +
  theme_minimal()  +
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
ggsave(filename = 'plot_memory.png', plot = mem, width = 12, height = 4, units = 'in', dpi = 600)

samples_n <- memory %>% select(sample, N)
memory_sarek <- read.csv('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/sarek_memory', sep = '\t', header = F) %>% 
  tidyr::separate(V2, sep = '\\/', into = c('spn', 'sarek', 'comb')) %>% 
  tidyr::separate(comb, sep = '_', into = c('cov', 'purity')) %>% 
  mutate(cov = ifelse(cov=='normal', '30x', cov)) %>% 
  mutate(type = ifelse(cov=='30x', 'normal', 'tumor')) %>% 
  select(-sarek) %>% 
  dplyr::rename(GB = V1) %>% 
  mutate(GB = as.numeric(stringr::str_replace(GB, pattern = 'G', '')),
         cov = stringr::str_replace(cov, pattern = 'x', ''),
         purity = stringr::str_replace(purity, pattern = 'p', '')) %>% 
  filter(spn %in% spns) %>% 
  left_join(samples_n %>% dplyr::rename(spn = sample))


memory_sarek %>% 
  filter(type != 'normal') %>% 
  ggplot()+
  geom_col(aes(x = spn, y = GB, fill = as.factor(N))) +
  scale_fill_manual('N samples', values = c('#7CCAD5', '#A0A6BE', '#C481A7', '#454995'))  +
  ylab('GB') +
  facet_grid(cov~purity) +
  theme_minimal() 
  