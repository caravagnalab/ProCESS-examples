library(VIBER)
library(patchwork)
library(tidyverse)

spn = 'SPN03'
source('/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/getters/process_getters.R')
source('/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/getters/tumourevo_getters.R')
samples <- get_sample_names(spn)
combs_df <- as.data.frame(t(combn(samples, 2)))

viber <- readRDS("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/subclonal/metrics_tables/viber_SPN03_100x_0.9p_mutect2_ascat.rds")

table_wide <- viber %>% 
  select(patient_id, sample_id, mutation_id, cluster_id_tool_interpreted, vaf_tool, is_driver_tool, driver_label_tool) %>%
  mutate(gene=driver_label_tool) %>% 
  pivot_wider(values_from="vaf_tool", names_from="sample_id")

table_wide <- table_wide %>%
  filter(!is.na(cluster_id_tool_interpreted))

table_wide <- table_wide %>%
  mutate(across(starts_with("Spn"), ~replace_na(., 0.0))) %>% 
  dplyr::rename(
    SPN03_1.1 = SPN03_SPN03_1.1,
    SPN03_2.1 = SPN03_SPN03_2.1,
    SPN03_3.1 = SPN03_SPN03_3.1,
    SPN03_4.1 = SPN03_SPN03_4.1
  )  %>% 
  separate(col = driver_label_tool, into = c('gene', 'other'), '_')

table_wide$cluster_id_tool_interpreted <- factor(
  table_wide$cluster_id_tool_interpreted,
  levels = sort(unique(table_wide$cluster_id_tool_interpreted))
)


plt <- c()
for (n in 1:nrow(combs_df)){
  s1 <- combs_df[n, ]$V1
  s2<- combs_df[n, ]$V2
  
  plt[[paste0('s_',n)]] <- ggplot() + 
    geom_point(data = table_wide, 
               aes(x = .data[[s1]], y =.data[[s2]], col = .data$cluster_id_tool_interpreted), 
               size = .5, alpha = .2) +
    scale_color_manual('Cluster', values = c(colors, 'Subclonal' = "steelblue4")) +
    theme_minimal() + 
    ylim(0,1) +
    xlim(0,1) +
    ggtitle('') +
    xlab(s1) +
    ylab(s2) + 
    geom_point(
      data = subset(table_wide, is_driver_tool == TRUE),
      aes(x = .data[[s1]], y = .data[[s2]]),
      color = "black", size = 1, shape = 15
    ) +
    ggrepel::geom_label_repel(
      data = subset(table_wide, is_driver_tool == TRUE),
      aes(x = .data[[s1]], y = .data[[s2]], label =.data$gene, color = .data$cluster_id_tool_interpreted),
      #color = 'black',
      size = 3,
      nudge_y = 0,
      nudge_x = 0,
      show.legend = FALSE
    )  +
    guides(color = guide_legend(override.aes = list(size = 3, alpha=1)))
}
p <- wrap_plots(plt, guides = 'collect') 
#p <- wrap_plots(p, tree, design = 'AAAAB')& theme(legend.position = 'bottom')
ggsave(filename = '/orfeo/cephfs/scratch/area/lvaleriani/tesi/tumourevo/subclonal_interpreted.png', 
       plot = p, width = 8, height = 5.5, units = 'in')



# signature
viber_int <- readRDS("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/subclonal/metrics_tables/viber_SPN03_100x_0.9p_mutect2_ascat.rds")
viber_true <- readRDS("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN03/tumourevo/100x_0.9p_mutect2_ascat/subclonal_deconvolution/viber/SCOUT/SPN03/SCOUT_SPN03_viber_best_st_fit.rds") 

data_table <- bind_cols(cluster = viber_true$x$cluster.Binomial, viber_true$data) %>% 
  select(cluster, chr, from, ref, alt) %>% 
  distinct() %>% 
  separate(col = chr, sep = 'chr', into = c('tmp', 'chr')) %>% 
  mutate(mutation_id = paste(spn, chr, from, alt, sep = ':')) %>% 
  select(-cluster, -tmp)

table <- viber_int %>% 
  left_join(data_table) %>% 
  select(chr,  from, ref, alt, cluster_id_tool_interpreted) %>% 
  rename(cluster = cluster_id_tool_interpreted) %>% 
  dplyr::mutate(Project = spn, Genome = 'GRCh38', mut_type = 'SNP', Type = 'SOMATIC', ID = cluster, Sample = cluster) %>%
  dplyr::rename(chrom = chr, pos_start = from) %>%
  rowwise() %>%
  dplyr::mutate(pos_end = pos_start + abs(str_count(ref) - str_count(alt))) %>%
  dplyr::select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref,alt, Type) %>%
  filter(ref != alt) %>% 
  distinct()


reticulate::use_python("/orfeo/scratch/area/lvaleriani/myconda/bin")
reticulate::py_config()
library(SigProfilerMatrixGeneratorR)
library(SigProfilerAssignmentR)

out_data ="/orfeo/cephfs/scratch/area/lvaleriani/tesi/tumourevo/sigprof/"
write.table(table, file = paste0(out_data, spn,'.txt'), quote = F, sep = '\t', row.names = F)
SigProfilerMatrixGeneratorR(project = spn, genome = "GRCh38", matrix_path = out_data, plot=T)

data_signature <- get_tumourevo_signatures(spn = spn, 
                                           coverage = 100, 
                                           purity = 0.9, 
                                           tool = 'SigProfiler',
                                           vcf_caller = 'mutect2',
                                           cna_caller = 'ascat',
                                           context = 'SBS96')
fit <- read.table(data_signature$COSMIC_signatures, header = T)
signature <- colnames(fit)[2:ncol(fit)]  

cosmic <- read.table('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/GRCh38/SBS_signatures.txt', header = T)
subset_cosmic <- cosmic[c('Type', signature)]
write.table(subset_cosmic, file = paste0(out_data, '/subset_cosmic.txt'), quote = F, sep = '\t', row.names = F)


data_signature <- get_tumourevo_signatures(spn = spn, 
                                           coverage = 100, 
                                           purity = 0.9, 
                                           tool = 'SigProfiler',
                                           vcf_caller = 'mutect2',
                                           cna_caller = 'ascat',
                                           context = 'ID83')
fit <- read.table(data_signature$COSMIC_signatures, header = T)
signature <- colnames(fit)[2:ncol(fit)]  

cosmic <- read.table('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/GRCh38/indel_signatures.txt', header = T)
subset_cosmic <- cosmic[c('Type', signature)]
write.table(subset_cosmic, file = paste0(out_data, '/subset_cosmic_ID.txt'), quote = F, sep = '\t', row.names = F)



data_id = paste0(out_data, '/output/SBS/', spn, '.SBS96.all')
cosmic_fit(samples = data_id, 
           output = out_data, 
           input_type='matrix', 
           context_type='SBS',
           collapse_to_SBS96=F, 
           cosmic_version=3.3, 
           genome_build="GRCh38", 
           signature_database=paste0(out_data, '/subset_cosmic.txt'),
           export_probabilities=TRUE,
           make_plots=TRUE)

activity <- read.table("/orfeo/cephfs/scratch/area/lvaleriani/tesi/tumourevo/sigprof/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
                       header = T) %>% 
  pivot_longer(
    cols = starts_with("SBS"),
    names_to = "Signature",
    values_to = "Exposure")

source('/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/validation/SCOUT/colors.R')
sbs <- activity %>% 
  mutate(Exposure = Exposure/sum(Exposure)) %>% 
  ggplot() +
  geom_col(aes(y = Samples, x = Exposure, fill = Signature)) + 
  scale_fill_manual('SBS', values = sbs_colors) + 
  ylab('Cluster') + 
  theme_minimal()+
  theme(text = element_text(size = 12))+ 
  xlim(0,1)




data_id = paste0(out_data, '/output/ID/', spn, '.ID83.all')
cosmic_fit(samples = data_id, 
           output = out_data, 
           input_type='matrix', 
           context_type='ID',
           collapse_to_SBS96=F, 
           cosmic_version=3.3, 
           genome_build="GRCh38", 
           signature_database=paste0(out_data, '/subset_cosmic_ID.txt'),
           export_probabilities=TRUE,
           make_plots=TRUE)

activity_id <- read.table("/orfeo/cephfs/scratch/area/lvaleriani/tesi/tumourevo/sigprof/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
                          header = T) %>% 
  pivot_longer(
    cols = starts_with("ID"),
    names_to = "Signature",
    values_to = "Exposure") %>% 
  mutate(Exposure = Exposure/sum(Exposure))

source('/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/validation/SCOUT/colors.R')
id <- activity_id %>% 
  ggplot() +
  geom_col(aes(y = Samples, x = Exposure, fill = Signature)) + 
  scale_fill_manual('ID', values = id_colors) + 
  ylab('Cluster') + 
  theme_minimal() +
  theme(text = element_text(size = 12)) + 
  xlim(0,1)

pp <- sbs + id + plot_layout(guides = 'collect')
ggsave(filename = '/orfeo/cephfs/scratch/area/lvaleriani/tesi/tumourevo/sign_interpreted.png', plot = pp,
       width = 8, height = 4, dpi = 400)

ggsave(filename = '/orfeo/cephfs/scratch/area/lvaleriani/tesi/tumourevo/sign_interpreted.pdf', plot = pp,
       width = 8, height = 4)