library(dplyr)
library(hms)
library(ggplot2)
library(patchwork)
library(lubridate)
library(optparse)

spns <- paste0('SPN0', 1:7)
sample_table <- tibble(sample = spns, 
                       N = c(3,2,4,2,3,5,5), 
                       hypermutant = c(F,T,F,F,T,F,T), 
                       type = c('WGD', 'Hypermutant', NA, NA, 'Hypermutant', 'WGD', 'Hypermutant'))

base <- '../../SCOUT/time/'
table <- tibble()
for (spn in spns){
  if (file.exists(file.path(base, paste0(spn, '_tissue.rds'))) & file.exists(file.path(base, paste0(spn, '_mutations.rds')))){
    print(spn)
    tissue <- readRDS(file.path(base, paste0(spn, '_tissue.rds'))) 
    muts <- readRDS(file.path(base, paste0(spn, '_mutations.rds'))) 
    tmp <- tibble(sample = spn, 
                  time_sample_forest =  tissue$cpu_time_secs, 
                  mem_sample_forest = tissue$memory_used_MB,
                  time_phylo_forest =  muts$cpu_time_secs, 
                  mem_phylo_forest = muts$memory_used_MB)
    table <- bind_rows(table, tmp)
  }
}

table <- table %>% left_join(sample_table)

plt <- table %>% 
  ggplot() +
  geom_col(aes(x = sample, y = hms::as_hms(time_sample_forest), fill = as.factor(N))) +
  ylab('time (H:M:S)') +
  scale_fill_manual('N samples', values = c('#7CCAD5', '#A0A6BE', '#C481A7', '#454995'))   +
  theme_minimal() + 
  xlab('SPN') + 
  ggtitle('Simulate tissue') + 
table %>% 
  ggplot() +
  geom_col(aes(x = sample, y = hms::as_hms(time_phylo_forest), fill = as.factor(type))) +
  ylab('time (H:M:S)') +
  xlab('SPN') + 
  scale_fill_manual('', values = c('Hypermutant' = 'maroon', 'NA' = 'slategrey', 'WGD' = 'steelblue')) + 
  theme_minimal() + 
  ggtitle('Simulate Mutations')

plt

ggsave(filename = 'plot_time_ProCESS.png', plot = plt, width = 12, height = 4, units = 'in', dpi = 600)
ggsave(filename = 'plot_time_ProCESS.pdf', plot = plt, width = 12, height = 4, units = 'in', dpi = 600)
saveRDS(object = table, file='time_ProCESS.rds')
