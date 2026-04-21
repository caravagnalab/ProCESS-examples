setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')
source("../figures/figure3/utils_plot.R")
source("../figures/figure3/utils.R")
source("../validation/SCOUT/colors.R")


indir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/assign_signature/"
indir_process= "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/assign_signature/"

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN02'),
                    make_option(c("--purity"), type = "double", default = 0.3),
                    make_option(c("--coverage"), type = "integer", default = 150),
                    make_option(c("--cna_caller"), type = "character", default = 'sequenza'),
                    make_option(c("--vcf_caller"), type = "character", default = 'mutect2'),
                    make_option(c("--signature"), type = "character", default = 'SigProfiler'),
                    make_option(c("--tool"), type = "character", default = 'pyclonevi')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

spn = opt$spn_id
cov = opt$coverage
pur = opt$purity
signature_tool =  opt$signature
cna_caller = opt$cna_caller
mut_caller = opt$vcf_caller


colors_cluster_mobster = c('hotpink4',
                           'indianred',
                           'indianred2', 
                           'aquamarine3',
                           'forestgreen',
                           'coral',
                           'hotpink4',
                           'lightgoldenrod',
                           'deepskyblue4',
                           'burlywood',
                           'cornsilk4')

names(colors_cluster_mobster) = paste0('C',0:10)
tail = c('Tail' = 'gray60')
colors_cluster_mobster = c(colors_cluster_mobster, tail)


colors_cluster = c('indianred', 
                   'steelblue4',  
                   'forestgreen', 
                   'goldenrod', 
                   'darkorange2', 
                   'palevioletred', 
                   'mediumpurple3', 
                   'palevioletred', 
                   'olivedrab3', 
                   'forestgreen',  
                   'cornflowerblue',
                   'aquamarine3',
                   'saddlebrown',
                   'deeppink2',
                   'cornflowerblue',
                   'black')
names(colors_cluster) = paste0('C',0:15)

tool = 'viber'
spn = 'SPN02'
cov = 150
pur = 0.9
mut_caller = 'mutect2'
cna_caller = 'ascat'

mutations <-  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal_new//tables_interpreted/', tool, '_', spn, '_', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '.rds')) %>% 
  separate(cluster_id_process, sep = '_', into = c('cluster_id_process', 'tmp')) %>% 
  select(-tmp)
samples <- paste(spn, get_sample_names(spn), sep= '_')
mutations <- mutations %>% mutate(driver_label = sub("_.*", "", driver_label_tool))

driver_tool <- mutations %>% select(patient_id, mutation_id, is_driver_tool, driver_label, is_driver_process)
muts_tool_raw <- mutations %>% 
  select(patient_id, mutation_id, cluster_id_tool, sample_id, vaf_tool) %>% 
  pivot_wider(values_from = vaf_tool, names_from = sample_id) %>% 
  replace(is.na(.), 0) %>% 
  mutate(cluster_id_tool = as.character(cluster_id_tool))
muts_tool_raw <- muts_tool_raw %>% left_join(driver_tool) %>% distinct()

muts_process <- mutations %>% 
  select(patient_id, mutation_id, sample_id, vaf_process, cluster_id_process) %>%
  pivot_wider(values_from = vaf_process, names_from = sample_id) %>% 
  replace(is.na(.), 0) 

driver_process <- mutations %>% select(patient_id, mutation_id, is_driver_process, driver_label)
muts_process <- muts_process %>% left_join(driver_process) %>% distinct()

if (tool == 'pyclonevi'){
  muts_tool_raw$cluster_id_tool = paste0('C', muts_tool_raw$cluster_id_tool)
}

s1 <- samples[[1]]
s2 <- samples[[2]]

color_palette_process = RColorBrewer::brewer.pal(n = max(3,length(unique(muts_process$cluster_id_process))), name = "Dark2") %>%
  setNames(str_sort(unique(muts_process$cluster_id_process), numeric=T))
color_palette_process['Subclonal'] = 'gainsboro'


p_tool = muts_tool_raw %>%
  ggplot(aes(x =.data[[s1]], y = .data[[s2]], color=cluster_id_tool))+
  geom_point( alpha=0.2, size = .5)+
  xlim(0,1)+
  ylim(0,1)+
  xlab(s1) +
  ylab(s2)+
  scale_color_manual('Tool clusters', values = colors_cluster)+
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1) ) )

p_tool = p_tool + ggrepel::geom_label_repel(
  data = muts_tool_raw %>% filter(is_driver_tool == TRUE),
  aes(
    x = .data[[s1]],
    y = .data[[s2]],
    label = driver_label,
    colour = cluster_id_tool, 
  ),
  show.legend = F,
  inherit.aes = FALSE,
  size = 2,
  min.segment.length = 0,
  box.padding = 1,
  max.overlaps = 50)

base = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/interpretation_', mut_caller, "_", cna_caller, '/')
s_tool = 'SigProfiler'
path = paste0(base, '/res_int/', cov, 'x_', pur, 'p_', s_tool,'_',tool)
sp =spn
score = readRDS(paste0(path, '_df.rds')) %>% filter(spn == sp)

sz = 4
plt <- score %>%
  mutate(cluster = factor(cluster, levels = paste0('C', 1:15))) %>% 
  pivot_longer(cols = c(score_driver, score_all, score_tail, score_no_driver, score_sign)) %>%
  mutate(name = factor(name, levels = c('score_driver', 'score_tail', 'score_sign','score_no_driver','score_all'))) %>%
  ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.9, ymax = 1, 
           fill = "palegreen4", alpha = 0.2) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.55, ymax = .9, 
           fill = "goldenrod", alpha = 0.2) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.2, ymax = 0.55, 
           fill = "salmon1", alpha = 0.2) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.2, ymax = 0, 
           fill = "gainsboro", alpha = 0.2) +
  geom_point(aes(x = name, y = value, col = cluster, shape = contains_driver_process), size = 3) +
  geom_line(data = ~ filter(.x, !is.na(value)), aes(x = name, y = value, col = cluster, group = cluster), linewidth = .7)  +
  #geom_text(data = ~ filter(.x, is_clonal_tool & name == 'score_all'), aes(x = name, y = value+0.07, col = cluster, group = cluster, label = 'Clonal'), size = sz) + 
  #geom_text(data = ~ filter(.x, contains_driver_tool & name == 'score_all'), aes(x = name, y = value+0.03, col = cluster, group = cluster, label = driver_label_process), size = sz) + 
  theme_minimal() +
  scale_shape_manual('Contains True Driver', values = c(4, 20)) +
  ylab('Score') +
  xlab('')+
  scale_color_manual('Cluster',
                     values = colors_cluster) +
  scale_x_discrete(labels = c('score_driver' = 'Driver',
                              'score_tail'   = 'Tail',
                              'score_sign'   = 'Signature',
                              'score_no_driver' = 'Signature\nTail',
                              #'score_no_tail' = 'Driver\nSignature',
                              #'score_no_sign' = 'Driver\nTail',
                              'score_all'    = 'All')) +
  my_ggplot_theme() 



plot_exposure <- function(df, sig_type, tool, col, type = '', signature = ''){
  plot <- df %>% 
    ggplot(aes(fill=Signature, y=Exposure, x=Samples)) + 
    geom_bar(position="fill", stat="identity")+
    coord_flip()+
    scale_fill_manual(values=col)+
    theme_bw()+
    theme(legend.position = "bottom")+
    xlab("Cluster")+
    ylab("Exposures")
  return(plot)
}

s_type = 'SBS96'
sign_type = 'SBS'
type = 'raw'
if (signature_tool == 'BASCULE'){
  name_file = paste0('bascule_fit_',type,'.rds')
} else{
  name_file = paste0('Assignment_Solution/Activities/Assignment_Solution_Activities.txt')
}
data = read.table(file = paste0(indir, spn, '/', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '/', tool, '_', signature_tool, '_', type, '/', s_type, '/',name_file),header = T) 
df = data %>% pivot_longer(
  cols = starts_with(sign_type),
  names_to = "Signature",
  values_to = "Exposure")
df <- df %>% 
  mutate(Samples = factor(Samples, levels = paste0('C', 1:12))) %>% 
  filter(Samples %in% score$cluster)
plot_sbs <- plot_exposure(df = df, sig_type = sign_type, tool = tool, type = tool, signature = signature_tool, col = c(sbs_colors, 'SBS58' ='darkseagreen'))

s_type = 'ID83'
sign_type = 'ID'
data = read.table(file = paste0(indir, spn, '/', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '/', tool, '_', signature_tool, '_', type, '/', s_type, '/',name_file),header = T) 
df = data %>% pivot_longer(
  cols = starts_with(sign_type),
  names_to = "Signature",
  values_to = "Exposure")
df <- df %>% 
  mutate(Samples = factor(Samples, levels = paste0('C', 1:12))) %>% 
  filter(Samples %in% score$cluster)
plot_id <- plot_exposure(df, sig_type = sign_type, tool = tool, type = tool, signature = signature_tool, col = id_colors)

# sz = 7
# all <- free(p_tool + theme(legend.position = 'none')) +  
#   plot_sbs + theme(legend.position = 'none') + 
#   plot_id  + theme(legend.position = 'none') + 
#   free(plt) + 
#   plot_layout(design = 'AABDDD\nAACDDD')
# ggsave(plot = all, width = 10, height = 2.8, units = 'in',
#        filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/spn02/spn02.png')
# ggsave(plot = all, width = 10, height = 2.8, units = 'in',
#        filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/spn02/spn02.pdf')


# sz = 7
# all <- free(p_tool + theme(legend.position = 'none')) +  
#   # plot_sbs + theme(legend.position = 'none') + 
#   # plot_id  + theme(legend.position = 'none') + 
#   free(plt) + 
#   plot_layout(design = 'AABBB')

# ggsave(plot = all, width = 7.5, height = 2.8, units = 'in',
#        filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/spn02/spn02_v2.png')
# ggsave(plot = all, width = 7.5, height = 2.8, units = 'in',
#        filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/spn02/spn02_v2.pdf')



# 
# tool = 'mobster_univariate'
# table_univariate = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal_new/tables_interpreted/', tool, '_', spn, '_', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '.rds')
# mutations <-  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal_new//tables_interpreted/', tool, '_', spn, '_', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '.rds')) %>% 
#   separate(cluster_id_process, sep = '_', into = c('cluster_id_process', 'tmp')) %>% 
#   select(-tmp)
# samples <- paste(spn, get_sample_names(spn), sep= '_')
# mutations <- mutations %>% mutate(driver_label = sub("_.*", "", driver_label_tool))
# 
# driver_tool <- mutations %>% select(patient_id, mutation_id, is_driver_tool, driver_label, is_driver_process)
# 
# id_muts <- mutations %>% 
#   dplyr::summarise(n = dplyr::n(), .by = c(patient_id, mutation_id, cluster_id_tool,
#                                            sample_id)) |>
#   dplyr::filter(n > 1L) 
# 
# muts_tool_raw <- mutations %>% 
#   filter(!mutation_id %in% id_muts$mutation_id) %>% 
#   select(patient_id, mutation_id, cluster_id_tool, sample_id, vaf_tool) %>% 
#   pivot_wider(values_from = vaf_tool, names_from = sample_id) %>% 
#   replace(is.na(.), 0) %>% 
#   mutate(cluster_id_tool = as.character(cluster_id_tool))
# muts_tool_raw <- muts_tool_raw %>% left_join(driver_tool) %>% distinct()
# 
# 
# cluster_plots_samples <- list()
# for (i in 1:length(samples)){
#   s1 <- samples[[i]]
#   if (s1 %in% colnames(muts_tool_raw)){
#     p_tool_raw = muts_tool_raw %>%
#       filter(.data[[s1]] != 0) %>% 
#       ggplot(aes(x =.data[[s1]], fill=cluster_id_tool))+
#       geom_histogram( alpha=0.6, binwidth = 0.005, position = "identity")+
#       xlim(-0.01,1.01)+
#       xlab(paste0('VAF ', s1)) +
#       ylab('') + 
#       scale_fill_manual('MOBSTER clusters', values = colors_cluster_mobster)+
#       #theme_bw() +
#       theme_void()
#     guides(fill = guide_legend(override.aes = list(size = 3, alpha = 1)))
#     cluster_plots_samples[[s1]] <- p_tool_raw
#   }
# }
# 
# uni <- wrap_plots(cluster_plots_samples$SPN02_SPN02_1.1 + theme(legend.position = 'none'),
#                   cluster_plots_samples$SPN02_SPN02_1.2 + theme(legend.position = 'none'), 
#                   guides = 'collect')
# ggsave(plot = uni, width = 3, height = .8, units = 'in',
#        filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/spn02/univariate_spn02.pdf')
# 
