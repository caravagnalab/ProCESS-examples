library(dplyr)
library(tidyverse)
library(GenomicRanges)
library(circlize)
library(ComplexHeatmap)
library(ProCESS)
library(ggrepel)
library(EnrichedHeatmap)

setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/cn_heatmap/')
source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/nf-validation/bin/cna/utils.R")
source('../../getters/process_getters.R')
source('utils.R')
source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/report/plotting/utils.R")

chr_level = paste0("chr", 1:22)
upper <- 0.89
lower <- 0.11

samples_levels <- lapply(c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN06', 'SPN07'), function(spn){
  samples <- get_sample_names(spn)
}) %>% unlist()

cna_data_all <- tibble()
for (spn in c('SPN01', 'SPN02', 'SPN04', 'SPN03', 'SPN06', 'SPN07')){ #'SPN02', 'SPN04', 'SPN03', 'SPN06', 'SPN07'
  samples <- get_sample_names(spn)
  cna_data <- lapply(samples, FUN = function(s){
    readRDS(get_process_cna(spn = spn, sample = s)) %>% mutate(sample = s)
  }) %>% bind_rows() 
  x <- cna_data %>% 
    dplyr::mutate(chr=paste0("chr",chr)) %>% 
    dplyr::rename("from"="begin") %>% 
    dplyr::rename("to"="end") %>% 
    dplyr::rename(Major=major) %>% 
    dplyr::rename(sample_id=sample) %>% 
    mutate(segment_id=paste(chr,from,to,sep=":")) %>% 
    mutate(CN_type = ifelse(ratio < upper & ratio > lower, 'sub-clonal', 'clonal')) %>% 
    filter(!(CN_type == 'sub-clonal' & ratio  > upper)) %>% 
    filter(!(CN_type == 'clonal' & ratio <= lower))  %>% 
    mutate(ratio = ifelse(ratio < upper & ratio > lower, ratio, 1)) %>% 
    mutate(ratio = round(ratio, digits=1)) %>% 
    filter(ratio != 0) 
  
  cna_data_all <- bind_rows(cna_data_all, x) 
}

out <- segment_fixed_windows_all(cna_data_all, window_size = 1e5) %>% 
  mutate(sample_id = factor(sample_id, levels = samples_levels))

out <-  out %>% 
  group_by(sample_id, segment_id) %>%
  filter(
    if (n_distinct(ratio) > 1) {
      # if there are different ratios, keep only Major:minor ≠ 1:1
      Major != 1 & minor != 1
    } else {
      # otherwise keep all rows
      TRUE
    }
  ) %>%
  ungroup() 


wide_df_kar_new <- out %>%
  filter(ratio == 1) %>% 
  select( -ratio, -segment_id) %>% 
  distinct() %>% 
  mutate(CN=paste(Major, minor, sep = ':' )) %>% 
  mutate(sample_id=as.character(sample_id)) %>% 
  select(chr, from, to, sample_id, CN) %>% # Select only relevant columns
  pivot_wider(names_from = sample_id, values_from = CN, values_fn = max) %>% 
  as.data.frame() %>% 
  replace(is.na(.), 'NA')%>% 
  filter(chr %in% paste0("chr", 1:22))

wide_df_kar_new <- wide_df_kar_new %>% 
  mutate(chr = factor(chr, levels = chr_level)) #%>% arrange(chr)


relative_to_absolute_coordinates_new = function(data)
{
  reference_genome = CNAqc:::get_reference('GRCh38')
  
  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr
  
  data %>%
    mutate(from = from + vfrom[chr],
           to = to + vfrom[chr])
}

# df <- out %>%
#   filter(ratio == 1) %>% 
#   select( -ratio, -segment_id) %>% 
#   distinct() %>% 
#   mutate(CN=paste(Major, minor, sep = ':' )) %>% 
#   mutate(sample_id=as.character(sample_id)) %>% 
#   select(chr, from, to, sample_id, CN) %>%
#   mutate(CN = ifelse(CN == '1:1', 0, 1)) %>% 
#   filter(CN != 0) %>% 
#   filter(!chr%in%c('chrX', 'chrY') )
# 
# df_abs <- relative_to_absolute_coordinates_new(df)
# 
# chrom_lengths <- df_abs %>%
#   mutate(chr = factor(chr, levels = chr_level)) %>% 
#   arrange(chr) %>% 
#   group_by(chr) %>%
#   summarise(chr_end = max(to), .groups = "drop",
#             chr_start = min(from)) %>%
#   mutate(
#     chr_mid   = (chr_start + chr_end)/2
#   ) %>% 
#   filter(!is.na(chr))
# 
# chrom_lengths <- chrom_lengths %>%
#   mutate(shade = as.factor(row_number() %% 2))
# 
# process_profile <- df_abs %>% 
#   ggplot() + 
#   geom_density(aes(x = from), bw = 1e7) + 
#   geom_vline(xintercept = chrom_lengths$chr_start[-1], color = "grey70", linetype = "dashed") +
#   scale_x_continuous(
#     breaks = chrom_lengths$chr_mid,
#     labels = chrom_lengths$chr
#   ) 

size = 1e5

df_ascat_list <- list()
df_battenberg_list <- list()
df_sequenza_list <- list()

df_ascat_all <- tibble()
df_battenberg_all <- tibble()
df_sequenza_all <- tibble()


for (c in c(100)){
  print(c)
  for (p in c(0.9)){
    print(p)
    
    i = 0
    for (spn in c('SPN01', 'SPN02', 'SPN04', 'SPN03', 'SPN06', 'SPN07')){ 
      samples = get_sample_names(spn)
      for (s in samples){
        print(s)
        if (file.exists(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION/',spn,'/cna/', c, 'x_',p,'/', s, '/data.rds'))){
          data = readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION/',spn,'/cna/', c, 'x_',p,'/', s, '/data.rds'))
          ascat = data$ascat_long %>% 
            select(-INFERRED_CN, INFERRED_ratio, -TRUE_CN, TRUE_ratio)  %>% 
            pivot_wider(names_from = Type, values_from = Value) %>% 
            select(chromosome, from, to, caller, is_match, type, 
                   TRUE_ratio, TRUE_Major1, TRUE_minor1, 
                   INFERRED_Major1, INFERRED_minor1,
                   TRUE_Major2, TRUE_minor2, 
                   INFERRED_Major2, INFERRED_minor2, INFERRED_ratio) %>% 
            mutate(sample_id = s)
          ascat_rel = relative_to_absolute_coordinates(ascat %>% dplyr::rename(chr = chromosome, start = from, end = to))
          ascat_new = segment_fixed_windows(ascat_rel %>% dplyr::rename(chromosome = chr, from = start, to = end), window_size = size)
          ascat_new_cn = ascat_new %>% mutate(is_match_new = ifelse(TRUE_Major1 ==1 & TRUE_minor1 ==1, '1:1', is_match)) %>% mutate(purity = p, coverage = c)
          match = ascat_new_cn$is_match_new
          df_ascat_all = bind_rows(df_ascat_all, ascat_new_cn)
          print('ASCAT done')
          
          ascat_chr = ascat %>% group_by(chromosome) %>% summarise(start = min(from), end = max(to))
          
          battenberg = data$Batt %>% 
            select(-INFERRED_CN, INFERRED_ratio, -TRUE_CN, TRUE_ratio)  %>% 
            pivot_wider(names_from = Type, values_from = Value) %>% 
            select(chromosome, from, to, caller, is_match, type, 
                   TRUE_ratio, TRUE_Major1, TRUE_minor1, 
                   INFERRED_Major1, INFERRED_minor1,
                   TRUE_Major2, TRUE_minor2, 
                   INFERRED_Major2, INFERRED_minor2, INFERRED_ratio) %>% 
            mutate(sample_id = s)
          battenberg_rel = relative_to_absolute_coordinates(battenberg %>% dplyr::rename(chr = chromosome, start = from, end = to))
          battenberg_new = segment_fixed_windows(battenberg_rel %>% dplyr::rename(chromosome = chr, from = start, to = end),window_size = size)
          battenberg_new_cn = battenberg_new %>% mutate(is_match_new = ifelse(TRUE_Major1 ==1 & TRUE_minor1 ==1, '1:1', is_match))
          battenberg_new_cn = battenberg_new_cn %>% right_join(ascat_new %>% select(chr, from, to)) %>% mutate(purity = p, coverage = c)
          match_battenberg = battenberg_new_cn$is_match_new
          df_battenberg_all = bind_rows(df_battenberg_all, battenberg_new_cn)
          print('Battenberg done')
          
          sequenza = data$sequenza_long %>%
            select(-INFERRED_CN, INFERRED_ratio, -TRUE_CN, TRUE_ratio)  %>%
            pivot_wider(names_from = Type, values_from = Value) %>%
            select(chromosome, from, to, caller, is_match, type,
                   TRUE_ratio, TRUE_Major1, TRUE_minor1,
                   INFERRED_Major1, INFERRED_minor1,
                   TRUE_Major2, TRUE_minor2,
                   INFERRED_Major2, INFERRED_minor2, INFERRED_ratio) %>%
            mutate(sample_id = s)
          
          
          sequenza_fixed <- sequenza %>%
            # Add chromosome bounds
            left_join(ascat_chr, by = "chromosome") %>%
            
            # Clamp to chromosome limits
            mutate(
              from = pmax(from, start),
              to   = pmin(to, end)
            ) %>%
            
            # Enforce valid interval ordering BEFORE grouping
            mutate(
              from_tmp = pmin(from, to),
              to_tmp   = pmax(from, to)
            ) %>%
            
            # Apply ordering per chromosome
            group_by(chromosome) %>%
            arrange(from_tmp, .by_group = TRUE) %>%
            
            # Force chromosome extremes
            mutate(
              from_tmp = if_else(row_number() == 1, start, from_tmp),
              to_tmp   = if_else(row_number() == n(), end, to_tmp)
            ) %>%
            
            ungroup() %>%
            
            # Final safety pass (guards against extreme-edge inversions)
            mutate(
              from = pmin(from_tmp, to_tmp),
              to   = pmax(from_tmp, to_tmp)
            ) %>%
            
            # Cleanup
            select(-start, -end, -from_tmp, -to_tmp)
          
          
          
          if (sum(!paste0('chr', 1:22) %in% unique(sequenza_fixed$chromosome)) != 0){
            value_columns <- c(
              "TRUE_ratio",
              "TRUE_Major1",
              "TRUE_minor1",
              "INFERRED_Major1"
            )
            
            sequenza_fixed <- add_missing_chromosomes(
              sequenza = sequenza_fixed,
              chr_reference = ascat_chr,
              value_cols = value_columns
            ) %>%
              mutate(sample_id = s)
          }
          
          sequenza_rel = relative_to_absolute_coordinates(sequenza_fixed %>% dplyr::rename(chr = chromosome, start = from, end = to))
          sequenza_new = segment_fixed_windows(sequenza_rel %>% dplyr::rename(chromosome = chr, from = start, to = end), window_size = size)
          sequenza_new_cn = sequenza_new %>% mutate(is_match_new = ifelse(TRUE_Major1 ==1 & TRUE_minor1 ==1, '1:1', is_match)) %>% mutate(purity = p, coverage = c)
          match_sequenza = sequenza_new_cn$is_match_new
          df_sequenza_all = bind_rows(df_sequenza_all, sequenza_new_cn)
          print('Sequenza done')
          
          if (i == 0){
            df_ascat = tibble(chr = ascat_new_cn$chr, from = ascat_new_cn$from, to = ascat_new_cn$to)
            df_battenberg = tibble(chr = battenberg_new_cn$chr, from = battenberg_new_cn$from, to = battenberg_new_cn$to)
            df_sequenza = tibble(chr = sequenza_new_cn$chr, from = sequenza_new_cn$from, to = sequenza_new_cn$to)
          }
          
          df_ascat[[`s`]] = match
          df_battenberg[[`s`]] = match_battenberg
          df_sequenza[[`s`]] = match_sequenza
          i = i + 1
        }
      }
    }
    tmp = paste(c, p, sep = ':')
    df_sequenza_list[[tmp]] = df_sequenza
    df_ascat_list[[tmp]] = df_ascat
    df_battenberg_list[[tmp]] = df_battenberg
  }
}

count_CN <- wide_df_kar_new %>%
  mutate(n_CN_event = rowSums(
    across(starts_with("SPN"), ~ {
      tmp <- as.character(.x)
      tmp != '1:1' & tmp != 'NA'
    }),
    na.rm = TRUE
  )) %>%
  ungroup() %>%
  select(chr, from, to, n_CN_event) #%>%
  #mutate(pos = 1:n())

i=0
for (tmp in names(df_ascat_list)){
  
  tmp_df <- df_ascat_list[[tmp]] %>% 
    rowwise() %>%
    mutate(n_match = sum(c_across(-c(chr, from, to)) == "match")) %>%
    ungroup() #%>% 
    #mutate(rel_value = n_match/count_CN$n_CN_event)
  
  if (i == 0){
    ascat_match = tibble(chr = tmp_df$chr, from = tmp_df$from, to = tmp_df$to)
  }
  ascat_match[[tmp]] <- tmp_df$n_match
  i = i+1
  
}

i=0
for (tmp in names(df_sequenza_list)){
  
  tmp_df <- df_sequenza_list[[tmp]] %>% 
    rowwise() %>%
    mutate(n_match = sum(c_across(-c(chr, from, to)) == "match")) %>%
    ungroup() #%>% 
    #mutate(rel_value = n_match/count_CN$n_CN_event)
  
  if (i == 0){
    sequenza_match = tibble(chr = tmp_df$chr, from = tmp_df$from, to = tmp_df$to)
  }
  sequenza_match[[tmp]] <- tmp_df$n_match
  i = i+1
}

i=0
for (tmp in names(df_battenberg_list)){
  
  tmp_df <- df_battenberg_list[[tmp]] %>% 
    rowwise() %>%
    mutate(n_match = sum(c_across(-c(chr, from, to)) == "match")) %>%
    ungroup()
  
  if (i == 0){
    battenberg_match = tibble(chr = tmp_df$chr, from = tmp_df$from, to = tmp_df$to)
  }
  battenberg_match[[tmp]] <- tmp_df$n_match
  i = i+1
  
}

battenberg_match <- battenberg_match %>% 
  rowwise() %>%
  mutate(mean_match = mean(c_across(-c(chr, from, to)))) 
  
sequenza_match <- sequenza_match %>% 
  rowwise() %>%
  mutate(mean_match = mean(c_across(-c(chr, from, to)))) 

ascat_match <- ascat_match %>% 
  rowwise() %>%
  mutate(mean_match = mean(c_across(-c(chr, from, to)))) 

#saveRDS(object = battenberg_match, 'battenberg_match.rds')
#saveRDS(object = sequenza_match, 'sequenza_match.rds')
#saveRDS(object = ascat_match, 'ascat_match.rds')

count_CN <- wide_df_kar_new %>%
  mutate(n_CN_event = rowSums(
    across(starts_with("SPN"), ~ {
      tmp <- as.character(.x)
      tmp != '1:1' & tmp != 'NA'
    }),
    na.rm = TRUE
  )) %>%
  ungroup() %>%
  select(chr, from, to, n_CN_event) #%>%
#mutate(pos = 1:n())


#saveRDS(object = wide_df_kar_new, 'wide_df_kar_new.rds')

join1 <- full_join(count_CN, ascat_match %>% select(chr, from, to, mean_match)) %>% 
  full_join(sequenza_match %>% select(chr, from, to, mean_match), by = join_by(chr, from, to), suffix = c('_ascat', '_sequenza')) %>% 
  full_join(battenberg_match %>% select(chr, from, to, mean_match), by = join_by(chr, from, to))  %>% 
  mutate(chr = factor(chr, levels = chr_level)) %>% 
  arrange(chr) %>% 
  mutate(pos = 1:n()) 


chrom_lengths <- join1 %>%
  mutate(chr = factor(chr, levels = chr_level)) %>% 
  arrange(chr) %>% 
  mutate(pos = 1:n()) %>% 
  group_by(chr) %>%
  summarise(chr_end = max(pos), .groups = "drop",
            chr_start = min(pos)) %>%
  mutate(
    chr_mid   = (chr_start + chr_end)/2
  ) 


chrom_lengths <- chrom_lengths %>%
  mutate(shade = as.factor(row_number() %% 2))


join = join1 %>% 
  dplyr::rename(ProCESS = n_CN_event, ASCAT = mean_match_ascat, Sequenza = mean_match_sequenza, Battenberg = mean_match) %>% 
  pivot_longer(cols = c(ASCAT, Sequenza,Battenberg)) %>% 
  mutate(rel = (value/ProCESS)*100) %>% 
  filter(rel <= 100)


         
rect_count <- count_CN %>%
  arrange(chr) %>% 
  mutate(pos = 1:n()) %>% 
  ggplot() +
  geom_rect(data = chrom_lengths,
            aes(xmin = chr_start, xmax = chr_end , ymin = 0, ymax = 13, fill = shade), alpha = 0.1, inherit.aes = FALSE) +
  scale_fill_manual(values = c("grey50", "white"), guide = "none") +
  geom_col(aes(x = pos, y = n_CN_event), color = "grey50") +
  geom_vline(xintercept = chrom_lengths$chr_start[-1], color = "grey70", linetype = "dashed") +
  scale_x_continuous(
    breaks = chrom_lengths$chr_mid,
    labels = chrom_lengths$chr
  ) +
  #ylim(0, 13) +
  labs(x = "", y = "# of CN events") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x  = element_blank(),
    plot.margin = margin(t = 2, r = 5, b = 2, l = 5)
  )  

df <- out %>%
  filter(ratio == 1) %>%
  select( -ratio, -segment_id) %>%
  distinct() %>%
  mutate(CN=paste(Major, minor, sep = ':' )) %>%
  mutate(sample_id=as.character(sample_id)) %>%
  select(chr, from, to, sample_id, CN) %>%
  mutate(CN = ifelse(CN == '1:1', 0, 1)) %>%
  filter(CN != 0) %>%
  filter(!chr%in%c('chrX', 'chrY') )

df_abs <- relative_to_absolute_coordinates_new(df)

chrom_lengths_abs <- df_abs %>%
  mutate(chr = factor(chr, levels = chr_level)) %>%
  arrange(chr) %>%
  group_by(chr) %>%
  summarise(chr_end = max(to), .groups = "drop",
            chr_start = min(from)) %>%
  mutate(
    chr_mid   = (chr_start + chr_end)/2
  ) %>%
  filter(!is.na(chr))

chrom_lengths_abs <- chrom_lengths_abs %>%
  mutate(shade = as.factor(row_number() %% 2))

density_count <- df_abs %>%
  ggplot() +
  geom_rect(data = chrom_lengths_abs,
            aes(xmin = chr_start, xmax = chr_end , ymin = 0, ymax = 6e-10, fill = shade), alpha = 0.1, inherit.aes = FALSE) +
  scale_fill_manual(values = c("grey50", "white"), guide = "none") +
  geom_density(aes(x = from), bw = 1e7) +
  geom_vline(xintercept = chrom_lengths_abs$chr_start[-1], color = "grey70", linetype = "dashed") +
  scale_x_continuous(
    breaks = chrom_lengths_abs$chr_mid,
    labels = chrom_lengths_abs$chr
  ) + 
  labs(x = "", y = "CN events") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x  = element_blank(),
    plot.margin = margin(t = 2, r = 5, b = 2, l = 5)
  )  

shade_profile <- join %>% 
  ggplot() + 
  geom_tile(
    aes(x = pos, y = name, fill = rel),
    height = 0.9
  ) +
  geom_vline(xintercept = chrom_lengths$chr_start[-1], color = "grey70", linetype = "dashed") +
  scale_x_continuous(
    breaks = chrom_lengths$chr_mid,
    labels = chrom_lengths$chr
  ) +
  scale_fill_gradientn(
    colors = c("peachpuff", "coral", "palegreen4"),
    #c("#f7fbff", "#6baed6", "#08306b")
    name = "Mean % Correct"
  ) + 
  labs(x = "Chromosome", 
       y = "Tool") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x  = element_blank(),
    plot.margin = margin(t = 2, r = 5, b = 2, l = 5)
  )  

v1 <- rect_count + shade_profile + plot_layout(nrow = 2, design = 'A\nB\nB')
v2 <- density_count + shade_profile + plot_layout(nrow = 2, design = 'A\nB\nB')


join2 <- join1 %>% 
  dplyr::rename(ProCESS = n_CN_event, ASCAT = mean_match_ascat, Sequenza = mean_match_sequenza, Battenberg = mean_match) %>% 
  mutate(ASCAT = ASCAT - 0.1) %>% 
  mutate(Sequenza = Sequenza - 0.2) %>% 
  mutate(Battenberg = Battenberg - 0.3) %>% 
  pivot_longer(cols = c(ProCESS, ASCAT, Sequenza,Battenberg))

line_profile <- join2 %>% 
  ggplot() + 
  geom_rect(data = chrom_lengths,
            aes(xmin = chr_start, xmax = chr_end , ymin = 0, ymax = 6e-10, fill = shade), alpha = 0.1, inherit.aes = FALSE) +
  scale_fill_manual(values = c("grey50", "white"), guide = "none") +
  geom_vline(xintercept = chrom_lengths$chr_start[-1], color = "grey70", linetype = "dashed") +
  geom_line(
    aes(x = pos, y = value, col = name),
    height = 0.9
  ) +
  scale_color_manual('Tool', values = c('ASCAT' = 'mediumpurple3', 'Sequenza' = 'darkorange2', 'Battenberg' = 'steelblue3', 'ProCESS' = 'gray'))+ 
  scale_x_continuous(
    breaks = chrom_lengths$chr_mid,
    labels = chrom_lengths$chr
  ) +
  labs(x = "Chromosome", 
       y = "# of CN events") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x  = element_blank(),
    plot.margin = margin(t = 2, r = 5, b = 2, l = 5)
  )  


line_profile2 <- join2 %>% 
  ggplot() + 
  # geom_rect(data = chrom_lengths,
  #           aes(xmin = chr_start, xmax = chr_end , ymin = 0, ymax = 6e-10, fill = shade), alpha = 0.1, inherit.aes = FALSE) +
  # scale_fill_manual(values = c("grey50", "white"), guide = "none") +
  geom_vline(xintercept = chrom_lengths$chr_start[-1], color = "grey70", linetype = "dashed") +
  geom_rect(
    aes(xmin = pos, xmax = pos, ymin = 0, ymax = value, color = name),
    fill = 'white'
  ) +
  scale_color_manual('Tool', values = c('ASCAT' = 'mediumpurple3', 'Sequenza' = 'darkorange2', 'Battenberg' = 'steelblue3', 'ProCESS' = 'gray'))+ 
  scale_x_continuous(
    breaks = chrom_lengths$chr_mid,
    labels = chrom_lengths$chr
  ) +
  labs(x = "Chromosome", 
       y = "# of CN events") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x  = element_blank(),
    plot.margin = margin(t = 2, r = 5, b = 2, l = 5)
  )  

v3 <- line_profile
v4 <- all_density
v5 <- upset
v6 <- line_profile2

ggsave(plot = v1,
       filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/segm_cna/plot_cn/v1.pdf',
       width = 11, height = 4)

ggsave(plot = v1,
       filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/segm_cna/plot_cn/v1.png',
       width = 11, height = 4, dpi = 400)

ggsave(plot = v2,
       filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/segm_cna/plot_cn/v2.pdf',
       width = 11, height = 4)

ggsave(plot = v3,
       filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/segm_cna/plot_cn/v3.pdf',
       width = 11, height = 3)

ggsave(plot = v4,
       filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/segm_cna/plot_cn/v4.pdf',
       width = 11, height = 3)

ggsave(plot = v5,
       filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/segm_cna/plot_cn/v5.pdf',
       width = 4, height = 3)

ggsave(plot = v6,
       filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/segm_cna/plot_cn/v6.pdf',
       width = 11, height = 3)
