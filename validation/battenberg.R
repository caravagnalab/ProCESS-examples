library(dplyr)
library(ProCESS)
library(optparse)
library(tidyr)
library(ggplot2)
library(patchwork)

setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/nf-validation/')

source("bin/getters/process_getters.R")
source("bin/getters/sarek_getters.R")
source("bin/cna/utils.R")

ascat_df_rel = readRDS('bin/data/ascat.rds') %>% dplyr::rename(chr=V1)
ascat_df_abs = absolute_to_relative_coordinates(ascat_df_rel)

data_dir = "/orfeo/scratch/cdslab/shared/SCOUT/"
spn_id = "SPN03"
coverage = "50"
purity = "0.9"
sample_id = "SPN03_2.1"

chromosomes = c(paste0('chr',1:22))
ProCESS_output = read_ProCESS(spn_id,sample_id,coverage,purity)
CNA_ProCESS = ProCESS_output[["CNA"]] 

CNA_ProCESS = CNA_ProCESS %>% 
  filter(ratio > 0.09) %>% 
  group_by(seg_id) %>% 
  mutate(ratio = ratio/sum(ratio)) %>% 
  #mutate(ratio = round(ratio,1)) %>% 
  filter(ratio != 0)  %>% 
  ungroup()

final_CNA_ProCESS = lapply(chromosomes, FUN = function(c){
  df = CNA_ProCESS %>% 
    filter(chr == c) 
  
  segments = unique(df$seg_id)
  new_df = df %>% filter(seg_id == segments[1]) %>% arrange(ratio)
  
  if (length(segments)>=2){
    for (s in 2:length(segments)){
      curr_seg = segments[s]
      curr_df = df %>% filter(seg_id == curr_seg)  %>% arrange(ratio)
      
      prev_seg = segments[s-1]
      prev_df = new_df %>% filter(seg_id == prev_seg) %>% arrange(ratio)
      
      if (identical(prev_df %>% select(major, minor, ratio), curr_df %>% select(major, minor, ratio))){
        tmp_df = curr_df %>% mutate(from = prev_df$from)
        new_df = new_df  %>% filter(seg_id != prev_seg)
        new_df = new_df %>% bind_rows(tmp_df)
      } else {
        new_df = new_df %>% bind_rows(curr_df)
      }
    }
  }
  return(new_df)
}) %>% bind_rows() %>% select(-seg_id)

# Compute real ploidy
ploidy = compute_true_ploidy(CNA_ProCESS)

# Compute FGA 
tot_genome = final_CNA_ProCESS  %>% 
  filter(!(chr %in% c('chrX', 'chrY')))  %>% 
  select(chr, from, to) %>% 
  distinct() %>%  
  mutate(len=to-from) %>% 
  pull(len) %>% 
  unique() %>% 
  sum()

altered = final_CNA_ProCESS %>% 
  filter(!(chr %in% c('chrX', 'chrY'))) %>% 
  mutate(len = to-from, CN = paste(major, minor, sep=':')) %>% 
  filter(ratio < 1 | CN !='1:1') %>% 
  select(-ratio, -CN, -major, -minor) %>% 
  distinct() %>% 
  pull(len) %>%
  unique() %>% 
  sum()
fga = (altered/tot_genome)*100

subclonal = final_CNA_ProCESS %>% 
  filter(!(chr %in% c('chrX', 'chrY'))) %>% 
  mutate(len = to-from, CN = paste(major, minor, sep=':')) %>% 
  filter(ratio < 1) %>% 
  select(-ratio, -CN, -major, -minor) %>% 
  distinct() %>% 
  pull(len) %>% 
  unique() %>% 
  sum()

fgs = (subclonal/tot_genome)*100
CNA_ProCESS  = final_CNA_ProCESS %>% mutate(type = ifelse(ratio != 1, 'subclonal', 'clonal'))

#### battenberg data
message("Reading Battenberg data")

file_cn <- file.path(data_dir, spn_id, 'sarek', paste0(coverage, 'x_', purity, 'p/variant_calling/battenberg/'), paste0(sample_id, '_vs_normal_sample/', paste0(sample_id, '_vs_normal_sample_subclones.txt')))
file_info <- file.path(data_dir, spn_id, 'sarek', paste0(coverage, 'x_', purity, 'p/variant_calling/battenberg/'), paste0(sample_id, '_vs_normal_sample/', paste0(sample_id, '_vs_normal_sample_rho_and_psi.txt')))
purity_ploidy_battenberg <- read.table(file_info, header = T) %>% filter(is.best == TRUE) %>% dplyr::rename(purity = rho) %>% select(purity, ploidy)
CNA_battenberg <- read.table(file_cn, header = T) %>% 
  select(chr, startpos, endpos, nMaj1_A, nMin1_A, frac1_A, nMaj2_A, nMin2_A, frac2_A) %>% 
  dplyr::rename(from = startpos, to = endpos) %>% 
  dplyr::rename(Major1 = nMaj1_A, minor1 = nMin1_A, Major2 = nMaj2_A, minor2 = nMin2_A)


message("Create joint table ProCESS and Battenberg calls") 
joint_segmentation_batt = create_joint_segmentation(CNA_ProCESS = final_CNA_ProCESS, CNA_target=CNA_battenberg, caller='battenberg', chromosomes = chromosomes)
joint_segmentation_batt_long = joint_segmentation_batt[['joint_segmentation_long']] %>% 
  filter(to - from > 1) %>% 
  left_join(ascat_df_abs %>% dplyr::rename(chromosome = chr)) %>% 
  filter(from >= start) %>%  
  filter(to <= end) %>% 
  select(-start, -end)

joint_segmentation_batt[['joint_segmentation']] = joint_segmentation_batt[['joint_segmentation']] %>% 
  filter(to - from > 1) %>% 
  left_join(ascat_df_abs %>% dplyr::rename(chromosome = chr)) %>% 
  filter(from >= start) %>%  
  filter(to <= end) %>% 
  select(-start, -end)

joint_segmentation_batt$joint_segmentation %>% filter(!is.na(TRUE_Major2))
joint_segmentation_batt$joint_segmentation_long %>% filter(type == 'subclonal')




CNA_ProCESS_rel = absolute_to_relative_coordinates(final_CNA_ProCESS %>% dplyr::rename(start = from, end = to))%>% 
  mutate(type = ifelse(ratio != 1, 'sub-clonal', 'clonal')) %>% 
  mutate(c = ifelse(ratio < 1 & ratio > 0.5, '1', '2')) %>% 
  mutate(c = ifelse(ratio == 1, 1, c)) %>% 
  pivot_longer(cols = c(major, minor)) %>% 
  mutate(name = ifelse(name == 'major', 'Major', 'minor')) %>% 
  mutate(name = paste0(name, c))  %>% 
  filter(chr %in% chromosomes)

process_annotations <- CNA_ProCESS_rel %>%
  filter(type == "sub-clonal") %>% 
  mutate(
    ratio_label = paste0(round(ratio, 1)),
    xpos = (start + end) / 2,
    ypos = .5  # Adjust as needed based on your y-axis scale
  ) %>% 
  group_by(chr, start, end, xpos, ypos) %>% 
  summarise(ratio_label = max(ratio_label))

ProCESS_plt =  CNAqc:::blank_genome(chromosomes = chromosomes) +
  geom_rect(data = CNA_ProCESS_rel, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill=type), alpha = .1) +
  geom_segment(data = CNA_ProCESS_rel %>% 
                 mutate(value = case_when(
                   name == 'Major1' ~ value + .04,
                   name == 'minor1' ~ value - .04,
                   name == 'Major2' ~ value + .10,
                   name == 'minor2' ~ value - .10,
                 )), aes(x=start, xend=end, y=value, col = name), size=1) +
  geom_text(data = process_annotations, aes(x = xpos, y = ypos, label = ratio_label), size = 3) +
  scale_color_manual('', values = color_process) + 
  scale_fill_manual('', values = color_type) + 
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.text=element_text(size=10)) + ggtitle('ProCESS')

file_cn <- file.path(data_dir, spn_id, 'sarek', paste0(coverage, 'x_', purity, 'p/variant_calling/battenberg/'), paste0(sample_id, '_vs_normal_sample/', paste0(sample_id, '_vs_normal_sample_subclones.txt')))
file_info <- file.path(data_dir, spn_id, 'sarek', paste0(coverage, 'x_', purity, 'p/variant_calling/battenberg/'), paste0(sample_id, '_vs_normal_sample/', paste0(sample_id, '_vs_normal_sample_rho_and_psi.txt')))
pur_ploidy <- read.table(file_info, header = T) %>% filter(is.best == TRUE) %>% dplyr::rename(purity = rho) %>% select(purity, ploidy)
batt_subclone <- read.table(file_cn, header = T) %>% select(chr, startpos, endpos, nMaj1_A, nMin1_A, frac1_A, nMaj2_A, nMin2_A, frac2_A) %>% 
  dplyr::rename(start = startpos, end = endpos)

batt_subclone_rel <- absolute_to_relative_coordinates(batt_subclone) %>%  
  dplyr::rename(Major1 = nMaj1_A, minor1 = nMin1_A, Major2 = nMaj2_A, minor2 = nMin2_A) %>% 
  pivot_longer(cols = c(Major1, minor1, Major2, minor2))  %>% 
  mutate(type = ifelse(frac1_A < 1, 'sub-clonal', 'clonal'))

batt_annotations <- batt_subclone_rel %>%
  filter(type == "sub-clonal") %>%
  filter(!is.na(frac1_A), !is.na(frac2_A), frac1_A != frac2_A) %>%
  mutate(
    ratio_label = paste0(round(frac1_A, 1), ":", round(frac2_A, 1)),
    xpos = (start + end) / 2,
    ypos = .5
  )

plt_battenberg <- CNAqc:::blank_genome(chromosomes = chromosomes) +
  geom_rect(data = batt_subclone_rel, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill=type), alpha = .1) +
  geom_segment(data = batt_subclone_rel %>% 
                 mutate(value = case_when(
                   name == 'Major1' ~ value + .04,
                   name == 'minor1' ~ value - .04,
                   name == 'Major2' ~ value + .10,
                   name == 'minor2' ~ value - .10,
                 )), aes(x=start, xend=end, y=value, col = name), size=1) +
    geom_text(data = batt_annotations, aes(x = xpos, y = ypos, label = ratio_label), size = 3) +
  scale_color_manual('', values = color_process) + 
  scale_fill_manual('', values = color_type) + 
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.text=element_text(size=10)) + ggtitle('Battenberg')

ProCESS_plt + plt_battenberg + plot_layout(nrow = 2)
