library(tidyverse)
library(patchwork)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")
source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/interpretation/mseq/utils.R")

base = "/orfeo/cephfs/scratch/area/lvaleriani/te_scripts/results_mseq_last/"
patient = 'Set06'
tool = 'viber'

# viber
subclonal_viber = readRDS(paste0(base, '/subclonal_deconvolution/viber/MSeq/', 
                         patient, '/MSeq_', patient, '_viber_best_fit.rds'))

ctree_viber = readRDS(paste0(base, '/subclonal_deconvolution/ctree/MSeq/', 
                                 patient, '/MSeq_', patient, '_ctree_VIBER.rds'))

data_ctree <- ctree_viber[[1]]$CCF %>% 
  select(cluster, is.clonal, is.driver, nMuts)

data_viber = subclonal_viber$data %>% 
  bind_cols(subclonal_viber$labels) %>% 
  mutate(mutation_id = paste(chr, from, alt, sep = ':')) %>% 
  mutate(cluster = cluster.Binomial) %>% 
  select(-chr, -from, -ref, -alt ) 

samples = list.dirs(paste0(base, '/subclonal_deconvolution/mobster/MSeq/', 
                           patient), full.names = F)

data_mobster = lapply(samples, FUN = function(s){
  if (s != ''){
    data = readRDS(paste0(base, '/subclonal_deconvolution/mobster/MSeq/',
                          patient, '/', s, '/MSeq_', patient, '_', s, '_mobster_best_fit.rds'))
    fit = data$data %>% select(chr, from, to, ref, alt, NV, DP, VAF, cluster) %>%
      mutate(mutation_id = paste(chr, from, alt, sep = ':'))  %>%
      mutate(sample = s)
  }
}) %>% bind_rows()

final_table_subclonal <- data_viber %>%
  left_join(data_mobster, by = join_by(mutation_id), suffix = c('_viber', '_mobster')) %>%
  group_by(mutation_id) %>%
  mutate(never_tail = all(cluster_mobster != "Tail")) %>%
  group_by(cluster_viber) %>%
  summarise(n_never_tail = sum(never_tail, na.rm = T)/n())  %>%
  full_join(data_ctree %>% dplyr::rename(cluster_viber = cluster))

samples <- subclonal_viber$y %>% colnames()
samples <- samples[2:ncol(subclonal_viber$y)]
pairs <- combn(samples[samples != ''], 2, simplify = FALSE)
cluster_plots_pair <- list()

for (i in 1:length(pairs)){
  s1 <- paste0('VAF.',pairs[[i]][1])
  s2 <- paste0('VAF.',pairs[[i]][2])
  
  plot = data_viber %>%
    ggplot(aes(x =.data[[s1]], y = .data[[s2]], color=cluster))+
    geom_point( alpha=0.2, size = .5)+
    xlim(0,1)+
    ylim(0,1)+
    xlab(s1) +
    ylab(s2)+
    scale_color_manual('Tool clusters', values = colors_cluster)+
    theme_bw() +
    ggtitle(tool) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1) ) )
  
  plot = plot + ggrepel::geom_label_repel(
    data = data_viber %>% filter(driver == TRUE),
    aes(
      x = .data[[s1]],
      y = .data[[s2]],
      label = gene,
      colour = cluster,
    ),
    show.legend = F,
    inherit.aes = FALSE,
    size = 3,
    min.segment.length = 0,
    box.padding = 1,
    max.overlaps = 50)
  
  cluster_plots_pair[[i]] <- plot
}

plt_multivariate <- wrap_plots(cluster_plots_pair, ncol = 3, guides = 'collect')


signature <- read.table(file = paste0(base, 'assign_signature/assign_cluster/MSeq/', patient, '/MSeq_', patient, '_viber/SBS/Assignment_Solution/Activities/Assignment_Solution_Activities.txt'),header = T) 
df = signature %>%
  pivot_longer(
    cols = starts_with('SBS'),
    names_to = "Signature",
    values_to = "Nmuts") %>%
  group_by(Samples) %>%
  mutate(Ntot = sum(Nmuts),
         Exposure = Nmuts/Ntot)

names(signature_colors) = unique(df$Signature)
plt_signature <- plot_exposure(df = df, sig_type = 'SBS', colors = signature_colors)


clonal <- final_table_subclonal %>% filter(is.clonal == T) %>% pull(cluster_viber)
df <- df %>%
  mutate(is_clonal = ifelse(Samples==clonal, T, F))

clonal_sig <- df %>% filter(is_clonal == T) %>% filter(Nmuts > 0) %>% dplyr::rename(n_sig = Nmuts)
other_sig <- df %>% filter(is_clonal == F) %>% filter(Nmuts > 0) %>% dplyr::rename(n_sig = Nmuts)

background = other_sig %>%
  group_by(Signature) %>%
  summarise(n_sig = sum(n_sig)) %>%
  mutate(n = sum(n_sig)) %>%
  mutate(Exposure = n_sig/n)

tmp_final_table <- lapply(unique(other_sig$Samples), FUN = function(c){
  tmp <- other_sig %>% filter(Samples == c)
  result_table <- compare_signatures(df1 = tmp, df2 = clonal_sig)
  result <- result_table$df %>% mutate(match = result_table$match,
                                       cs_exp = result_table$cs_exp,
                                       cs_nsig = result_table$cs_n,
                                       n_diff = result_table$n_diff,
                                       n_tot = result_table$n_tot) %>%
    mutate(cluster = as.character(c))
  return(result)
}) %>% bind_rows()


tmp_final_table_background <- lapply(unique(other_sig$Samples), FUN = function(c){
  tmp <- other_sig %>% filter(Samples == c)
  result_table <- compare_signatures(df1 = tmp, df2 = background)
  result <- result_table$df %>% mutate(match = result_table$match,
                                       cs_exp = result_table$cs_exp,
                                       cs_nsig = result_table$cs_n,
                                       n_diff = result_table$n_diff,
                                       n_tot = result_table$n_tot) %>%
    mutate(cluster = as.character(c))
  return(result)
}) %>% bind_rows()

if (nrow(tmp_final_table) > 0){

  f_background <- tmp_final_table_background %>%
    mutate(n_rel = n_diff/n_tot) %>%
    select(cluster, match, cs_exp, n_rel) %>%
    distinct() %>%
    dplyr::rename(bg_match = match, bg_cs_exp =cs_exp, bg_n_rel = n_rel) %>%
    filter(!is.na(Samples))

  final_table_signature <- tmp_final_table %>%
    mutate(n_rel = n_diff/n_tot) %>%
    select(cluster, match, cs_exp, n_rel) %>%
    distinct() %>%
    left_join(f_background)

}


final_table <- final_table_subclonal %>%
  left_join(final_table_signature %>% dplyr::rename(cluster_viber= cluster) %>% select(-Samples)) %>%
  filter(nMuts > 100) %>%
  mutate(driver = ifelse(is.driver == F, 0, 1),
         bg_cs = 1 - bg_cs_exp,
         cs_sign = 1 - cs_exp,
         cs_sign = ifelse(is.na(cs_exp) & is.clonal == T, 1, cs_sign),
         bg_sign = ifelse(is.na(bg_cs) & is.clonal == T, 1, bg_cs),
         n_rel = ifelse(is.na(cs_exp) & is.clonal == T, 1, n_rel)) %>%
  rowwise() %>%
  mutate(
         score_all = (driver + n_never_tail + ((cs_sign+n_rel+bg_sign)/3))/3,
         score_driver = driver,
         score_sign = (cs_sign+n_rel+bg_sign)/3,
         score_tail = n_never_tail,
         score_no_tail = (driver + ((cs_sign+n_rel+bg_sign)/3))/2,
         score_no_driver = (((cs_sign+n_rel+bg_sign)/3) + n_never_tail)/2,
         score_no_sign = (driver + n_never_tail)/2)

plt <- final_table %>%
  pivot_longer(cols = c(score_driver, score_all, score_tail, score_no_driver, score_no_tail, score_no_sign, score_sign)) %>%
  mutate(name = factor(name, levels = c('score_driver', 'score_tail', 'score_sign','score_no_driver', 'score_no_tail', 'score_no_sign', 'score_all'))) %>%
  ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.9, ymax = 1,
           fill = "palegreen4", alpha = 0.2) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.55, ymax = .9,
           fill = "goldenrod", alpha = 0.2) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.2, ymax = 0.55,
           fill = "salmon1", alpha = 0.2) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.2, ymax = 0,
           fill = "gainsboro", alpha = 0.2) +
  geom_point(aes(x = name, y = value, col = cluster_viber, shape = is.driver), size = 4) +
  geom_line(data = ~ filter(.x, !is.na(value)), aes(x = name, y = value, col = cluster_viber, group = cluster_viber), linewidth = .6)  +
  geom_text(data = ~ filter(.x, is.clonal & name == 'score_all'), aes(x = name, y = value+0.07, col = cluster_viber, group = cluster_viber, label = 'Clonal')) +
  theme_minimal() +
  scale_shape_manual('Contains Driver', values = c(4, 20)) +
  xlab('')+
  scale_color_manual('Cluster',
                     values = colors_cluster) +
  scale_x_discrete(labels = c('score_driver' = 'Driver',
                              'score_tail'   = 'Tail',
                              'score_sign'   = 'Signature',
                              'score_no_driver' = 'Signature\nTail',
                              'score_no_tail' = 'Driver\nSignature',
                              'score_no_sign' = 'Driver\nTail',
                              'score_all'    = 'All'))

# plt <- ggplot()
# plt_signature <- ggplot()
if (patient == 'Set06'){
  final_plt <- wrap_plots(plt_multivariate, plt, plt_signature) + plot_layout(design = 'AA\nAA\nAA\nAA\nAA\nAA\nAA\nAA\nBC')
  h = 18
} else {
  final_plt <- wrap_plots(plt_multivariate, plt, plt_signature) + plot_layout(design = 'AA\nAA\nAA\nBC')
  h = 10
}
ggsave(filename = paste0('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/interpretation/mseq/MSeq_',patient,'_viber.png'),
       plot = final_plt,
       width = 10,
       height = h)

