library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
library(rstatix)

source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation/SCOUT/colors.R')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils.R")

indir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/assing_signature/"

plot_exposure <- function(df, sig_type, tool, col, type = '', signature = ''){
  plot <- df %>% 
    ggplot(aes(fill=Signature, y=Exposure, x=as.factor(Samples))) + 
    geom_bar(position="fill", stat="identity")+
    coord_flip()+
    scale_fill_manual(values=col)+
    theme_bw()+
    theme(legend.position = "bottom")+
    xlab("Cluster")+
    ylab("Exposures")+
    ggtitle(label = paste(sig_type, type, signature))
  return(plot)
}

spn='SPN07'
cov = 50
pur = 0.9
tool = 'viber'
sign_type = 'SBS'
sig_tool = 'BASCULE'
signature_tool = 'BASCULE'
s_type = 'SBS96'
colors = sbs_colors
type = 'int'
mut_caller='mutect2'
cna_caller = 'ascat'

process = readRDS(paste0(indir, spn, '/process/',cov, 'x_', pur, 'p/exposure_', sign_type, '.rds'))
process_df = process %>% 
  select(cluster_id_process, causes, nmuts_cause, tot_nmuts, exposure) %>% 
  dplyr::rename(Signature=causes,
                Exposure=exposure,
                Samples=cluster_id_process)

name_file = paste0('bascule_fit_',type,'.rds')
data = readRDS(paste0(indir, spn, '/', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '/', tool, '_', signature_tool, '_', type, '/', s_type, '/',name_file))

df = data[['nmf']][[sign_type]][['exposure']] %>% 
  dplyr::rename(Signature=sigs,
                Exposure=value,
                Samples=samples)

counts = data[['input']][[sign_type]][['counts']] %>% 
  group_by(samples) %>% 
  summarise(n = sum(value)) %>% 
  dplyr::rename(Samples = samples)

sp = spn
t = tool
result_int <- readRDS('/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/metrics_tables/metrics_drivers_clonal_vs_subclonal.rds') %>% 
  filter(purity == pur, 
         coverage == cov, 
         spn == sp, 
         tool == t)
tp_cluster = result_int$TP_c_list %>% unlist()
df_tp = df %>% filter(Samples %in% tp_cluster)

mutations <-  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal/tables_interpreted/', tool, '_', spn, '_', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '.rds'))
annot_table <- mutations %>% 
  # filter(cluster_id_tool %in% tp_cluster,
  #        is_driver_tool == T | is_driver_process == T) %>%
  select(cluster_id_tool, driver_label_tool, driver_label_process, cluster_id_process,is_driver_tool, is_driver_process) %>% 
  unique() %>% 
  tidyr::separate(driver_label_tool, into = c('driver_tool', 'tmp_tool'), sep = '_') %>% 
  tidyr::separate(driver_label_process, into = c('driver_process', 'tmp_process'), sep = ' ')


tool_cluster = annot_table %>% filter(driver_tool %in% annot_table$driver_process) %>% select(cluster_id_tool, driver_tool, driver_process)#%>% filter(is_driver_tool == T) %>% select(cluster_id_tool, driver_tool, driver_process)
process_cluster = annot_table %>% filter(driver_process %in% tool_cluster$driver_process) %>% filter(is_driver_process == T) %>% select(cluster_id_process, driver_process)

exp_process = process_df %>% filter(Samples %in% process_cluster$cluster_id_process) %>% left_join(process_cluster %>% dplyr::rename(Samples = cluster_id_process))
exp_tool = df %>% filter(Samples %in%  tool_cluster$cluster_id_tool) %>% left_join(tool_cluster %>% mutate(cluster_id_tool = as.character(cluster_id_tool)) %>%  dplyr::rename(Samples = cluster_id_tool))

color_palette_driver =   RColorBrewer::brewer.pal(n = max(3,length(unique(exp_process$driver_process))), name = "Dark2") %>%
  setNames(str_sort(unique(exp_process$driver_process), numeric=T))
color_palette_driver['Subclonal'] = 'gainsboro'
process_plot <- plot_exposure(exp_process %>% select(-Samples) %>% dplyr::rename(Samples = driver_process), sig_type = paste0(sign_type, ' ProCESS'), tool = 'ProCESS', col = colors)

color_palette_driver =   RColorBrewer::brewer.pal(n = max(3,length(unique(exp_process$driver_process))), name = "Dark2") %>%
  setNames(str_sort(unique(exp_process$driver_process), numeric=T))
color_palette_driver['Subclonal'] = 'gainsboro'
tool_plot <- plot_exposure(exp_tool %>% select(-Samples) %>% dplyr::rename(Samples = driver_tool), sig_type = paste0(sign_type, ' VIBER'), tool = 'Tool', col = colors)

plt_exp = process_plot + tool_plot + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

sp = spn
t = tool
result_int <- readRDS('/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/metrics_tables/metrics_drivers_clonal_vs_subclonal.rds') %>% 
  filter(purity == pur, 
         coverage == cov, 
         spn == sp, 
         tool == t)
tp_cluster = result_int$TP_c_list %>% unlist()

mutations <-  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal/tables_interpreted/', tool, '_', spn, '_', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '.rds'))
samples <- unique(mutations$sample_id)
pairs <- combn(samples, 2, simplify = FALSE)
annot_table <- mutations %>% 
  filter(cluster_id_tool %in% tp_cluster, 
         is_driver_tool == T | is_driver_process == T) %>% 
  select(cluster_id_tool, driver_label_tool, driver_label_process, cluster_id_process,is_driver_tool, is_driver_process) %>% 
  unique() %>% 
  tidyr::separate(driver_label_tool, into = c('driver_tool', 'tmp_tool'), sep = '_') %>% 
  tidyr::separate(driver_label_process, into = c('driver_process', 'tmp_process'), sep = ' ')


tool_cluster = annot_table %>% filter(driver_tool %in% annot_table$driver_process) %>% filter(is_driver_tool == T) %>% select(cluster_id_tool, driver_tool, driver_process)
process_cluster = annot_table %>% filter(driver_process %in% tool_cluster$driver_process) %>% filter(is_driver_process == T) %>% select(cluster_id_process, driver_process)

dups <- mutations %>% 
  mutate(cluster_id_process = ifelse(cluster_id_process %in% process_cluster$cluster_id_process, cluster_id_process, 'NA')) %>% 
  select(patient_id, mutation_id, sample_id, cluster_id_process, vaf_process, is_driver_process) %>% 
  left_join(process_cluster) %>% 
  replace(is.na(.), 'Other') %>% 
  dplyr::summarise(n = dplyr::n(), .by = c(patient_id, mutation_id, cluster_id_process, is_driver_process, driver_process,
                                           sample_id)) %>% 
  dplyr::filter(n > 1L)


mutations_process <- mutations %>% 
  mutate(cluster_id_process = ifelse(cluster_id_process %in% process_cluster$cluster_id_process, cluster_id_process, 'NA')) %>% 
  select(patient_id, mutation_id, sample_id, cluster_id_process, vaf_process, is_driver_process) %>% 
  left_join(process_cluster) %>% 
  replace(is.na(.), 'Other')  %>% 
  filter(!mutation_id %in% dups$mutation_id) %>% 
  pivot_wider(values_from = vaf_process, names_from = sample_id) %>% 
  replace(is.na(.), 0) 




mutations_tool <- mutations %>% 
  mutate(cluster_id_tool_interpreted = ifelse(cluster_id_tool_interpreted %in% tool_cluster$cluster_id_tool, cluster_id_tool_interpreted, 'NA')) %>% 
  select(patient_id, mutation_id, sample_id, cluster_id_tool_interpreted, vaf_tool, is_driver_tool, driver_label_tool) %>% 
  dplyr::rename(cluster_id_tool = cluster_id_tool_interpreted) %>% 
  left_join(tool_cluster %>% mutate(cluster_id_tool = as.character(cluster_id_tool))) %>% 
  replace(is.na(.), 'Other')  %>% 
  filter(!mutation_id %in% dups$mutation_id) %>% 
  pivot_wider(values_from = vaf_tool, names_from = sample_id) %>% 
  replace(is.na(.), 0) 


join = tool_cluster %>% left_join(process_cluster) %>% filter(!is.na(cluster_id_process))

tp_tool_muts <- mutations_tool %>% filter(cluster_id_tool %in% unique(tool_cluster$cluster_id_tool))
tp_process_muts <- mutations_process  %>% filter(cluster_id_process %in% unique(process_cluster$cluster_id_process))

muts_tool = tp_tool_muts %>% select(mutation_id, driver_tool)
muts_process = tp_process_muts %>% select(mutation_id, driver_process)



library(dplyr)
# Remove duplicated mutation_id mappings
mut_links <- mutations %>%
  distinct(mutation_id, cluster_id_tool, cluster_id_process) %>% 
  filter(cluster_id_tool %in% tp_cluster) %>% 
  filter(cluster_id_tool != 'C10') %>% 
  mutate(cluster_id_tool = case_when(
    cluster_id_tool == 'C1' ~ 'TP53',
    cluster_id_tool == 'C4' ~ 'ATRX',
    cluster_id_tool == 'C5' ~ 'NF1'
  )) %>% 
  mutate(cluster_id_process = case_when(
    cluster_id_process == 'Clone 3' ~ 'NF1',
    cluster_id_process == 'Clone 4' ~ 'ATRX',
    cluster_id_process == 'Clone 6' ~ 'TP53',
    cluster_id_process == 'Clonal' ~ 'Clonal',
    cluster_id_process == 'Clone 5' ~ 'Clone 5',
    cluster_id_process == 'Subclonal' ~ 'Subclonal'
  ))

sankey_edges <- mut_links %>%
  count(cluster_id_tool, cluster_id_process, name = "value")

ggplot(sankey_edges,
       aes(axis2 = cluster_id_tool,
           axis1 = cluster_id_process,
           y = value)) +
  geom_alluvium(aes(fill = cluster_id_tool), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("ProCESS cluster", "VIBER cluster"),
                   expand = c(.05, .05)) +
  labs(y = "Number of mutations") + 
  my_ggplot_theme() +
  scale_fill_manual('Cluster', values = c('dodgerblue4', 'salmon3', 'palegreen4')) +
  scale_color_manual(values = c('dodgerblue4', 'salmon3', 'palegreen4')) + 
  theme(legend.position = 'bottom') 

