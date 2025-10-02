library(dplyr)
library(hms)
library(ggplot2)
library(patchwork)
library(lubridate)
library(optparse)
source('utils.R')
spns <- paste0('SPN0', 1:7)
sample_table <- tibble(sample = spns, 
                       N = c(3,2,4,2,3,5,5), 
                       hypermutant = c(F,T,F,F,F,F,T), 
                       type = c('WGD', 'Hypermutant', NA, NA, NA, 'WGD', 'Hypermutant'))

memory_phylo <- read.csv('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/phylo_process', sep = '\t', header = F) %>% 
  tidyr::separate(V2, sep = '\\/', into = c('spn', 'name', 'sff')) %>% 
  select(-name, -sff) %>% 
  filter(spn %in% spns) %>% 
  dplyr::rename(MB = V1)  %>% 
  mutate(MB = as.numeric(MB)) %>% 
  left_join(sample_table %>% dplyr::rename(spn = sample))

memory_sample <- read.csv('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/sample_process', sep = '\t', header = F) %>% 
  tidyr::separate(V2, sep = '\\/', into = c('spn', 'name', 'sff')) %>% 
  select(-name, -sff) %>% 
  filter(spn %in% spns) %>% 
  dplyr::rename(KB = V1)  %>% 
  mutate(KB = as.numeric(KB)) %>% 
  left_join(sample_table %>% dplyr::rename(spn = sample))

file_memory_normal <- read.csv('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/fastq_normal', sep = '\t', header = F) %>% 
  tidyr::separate(V2, sep = '\\/', into = c('spn', 'seq', 'type', 'purity')) %>% 
  tidyr::separate(purity, sep = '_', into = c('name', 'purity')) %>% 
  select(-name, -seq) %>% 
  dplyr::rename(GB = V1)  %>% 
  mutate(GB = as.numeric(stringr::str_replace(GB, pattern = 'G', ''))) %>% 
  filter(spn %in% spns) %>% 
  mutate(N = 1)

file_memory_tumour <- read.csv('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/fastq_process', sep = '\t', header = F) %>% 
  tidyr::separate(V2, sep = '\\/', into = c('spn', 'seq', 'type', 'purity')) %>% 
  tidyr::separate(purity, sep = '_', into = c('name', 'purity')) %>% 
  select(-name, -seq) %>% 
  dplyr::rename(GB = V1) %>% 
  filter(spn %in% spns) %>% 
  mutate(GB = as.numeric(GB)) %>% 
  left_join(sample_table %>% dplyr::rename(spn = sample) %>% select(-type))

memory_process <- bind_rows(file_memory_tumour, file_memory_normal)

mem_fastq <- memory_process %>% 
  filter(spn != 'SPN05') %>% 
  mutate(N= ifelse(N==1, '1 - Normal', N)) %>% 
  ggplot() +
  geom_col(aes(x = spn, y = GB/1000, fill = as.factor(N))) +
  scale_fill_manual('N samples', values = c('2' = '#c9e4ca', '3' ='#87bba2', '4' ='#55828b', '5' ='#3b6064', '1 - Normal' = 'gray'))  +
  ylab('TB') + 
  xlab('SPN') +
  ggtitle('FASTQ - 200X x 3p') + 
  theme_minimal() 

mem_phylo <- memory_phylo %>% 
  filter(spn != 'SPN05') %>% 
  ggplot() +
  geom_col(aes(x = spn, y = MB/1000, fill = as.factor(type))) +
  #scale_fill_manual('Hypermutant', values = c('TRUE' = 'maroon', 'FALSE' = 'slategrey'))  +
  scale_fill_manual('', values = c('Hypermutant' = '#ce796b', 'NA' = 'slategrey', 'WGD' = '#577399')) + 
  ylab('GB') + 
  xlab('SPN') +
  ggtitle('Phylogenetic Forest') + 
  theme_minimal() 

mem_sample <- memory_sample %>% 
  filter(spn != 'SPN05') %>% 
  bind_rows(tibble(MB = 0, spn = 'SPN01', N=1)) %>% 
  mutate(N= ifelse(N==1, '1 - Normal', N)) %>% 
  ggplot() +
  geom_col(aes(x = spn, y = KB/1000, fill = as.factor(N))) +
  scale_fill_manual('N samples', values = c('2' = '#c9e4ca', '3' ='#87bba2', '4' ='#55828b', '5' ='#3b6064', '1 - Normal' = 'gray'))  +
  #scale_fill_manual('N samples', values = c('#c9e4ca', '#87bba2', '#55828b', '#3b6064'))  +
  ylab('MB') + 
  xlab('SPN') +
  ggtitle('Sample Forest') + 
  theme_minimal() 

mm <- mem_sample + mem_phylo + mem_fastq + plot_layout(guides = 'collect')
mm
ggsave(filename = 'plot_memory_ProCESS.png', plot = mm, width = 12, height = 3, units = 'in', dpi = 600)
ggsave(filename = 'plot_memory_ProCESS.pdf', plot = mm, width = 12, height = 3, units = 'in', dpi = 600)

memory = list('sample_forest' =  memory_sample %>% filter(spn != 'SPN05'), 'phylo_forest' = memory_phylo%>%filter(spn != 'SPN05'), 'fastq' = memory_process %>%filter(spn != 'SPN05'))
saveRDS(object = memory, file='memory_ProCESS.rds')
