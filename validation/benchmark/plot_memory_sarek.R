library(dplyr)
library(hms)
library(ggplot2)
library(patchwork)
library(lubridate)
library(optparse)
source('utils.R')
spns <- paste0('SPN0', 1:7)

sample_table <- tibble(sample = spns, N = c(3,2,4,2,3,5,5))

#du -hs SPN0*/sarek/*x_0.* > sarek_memory
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
  left_join(sample_table %>% dplyr::rename(spn = sample))


plt_sarek <- memory_sarek %>% 
  group_by(spn, cov, type, N) %>% 
  summarize(GB = sum(as.numeric(GB))) %>% 
  filter(type != 'normal') %>% 
  ggplot()+
  geom_col(aes(x = spn, y = GB/1000, fill = as.factor(N))) +
  scale_fill_manual('N samples', values = c('#7CCAD5', '#A0A6BE', '#C481A7', '#454995'))  +
  ylab('TB') +
  facet_grid(.~cov) +
  ggtitle('nf-core/sarek') +
  xlab('SPN') + 
  theme_minimal() +
  plot_annotation(caption = 'Data from 3 purities and normal')
ggsave(filename = 'plot_sarek_memory.png', plot = plt_sarek, width = 8, height = 4, units = 'in', dpi = 600)
ggsave(filename = 'plot_sarek_memory.pdf', plot = plt_sarek, width = 8, height = 4, units = 'in', dpi = 600)



plt_sarek_all <- memory_sarek %>% 
  group_by(spn, type, N) %>% 
  summarize(GB = sum(as.numeric(GB))) %>% 
  filter(type != 'normal') %>% 
  ggplot()+
  geom_col(aes(x = spn, y = GB/1000, fill = as.factor(N))) +
  scale_fill_manual('N samples', values = c('#7CCAD5', '#A0A6BE', '#C481A7', '#454995'))  +
  ylab('TB') +
  ggtitle('nf-core/sarek') +
  xlab('SPN') + 
  theme_minimal() +
  plot_annotation(caption = 'Data from 50x and 100x coverage, 3 purities and normal')
ggsave(filename = 'plot_sarek_memory_all.png', plot = plt_sarek_all, width = 5, height = 4, units = 'in', dpi = 600)
ggsave(filename = 'plot_sarek_memory_all.pdf', plot = plt_sarek_all, width = 5, height = 4, units = 'in', dpi = 600)

saveRDS(object = memory_sarek, file='memory_sarek.rds')
