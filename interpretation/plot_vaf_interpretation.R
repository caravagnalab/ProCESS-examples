setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
library(patchwork)
source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

indir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/assing_signature/"

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN02'),
                    make_option(c("--purity"), type = "double", default = 0.9),
                    make_option(c("--coverage"), type = "integer", default = 100),
                    make_option(c("--cna_caller"), type = "character", default = 'ascat'),
                    make_option(c("--vcf_caller"), type = "character", default = 'mutect2'),
                    make_option(c("--signature"), type = "character", default = 'BASCULE'),
                    make_option(c("--tool"), type = "character", default = 'viber')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

spn = opt$spn_id
cov = opt$coverage
pur = opt$purity
signature_tool =  opt$signature
cna_caller = opt$cna_caller
mut_caller = opt$vcf_caller


colors_cluster = c('indianred', 
                   'steelblue', 
                   'forestgreen', 
                   'goldenrod', 
                   'darkorange3', 
                   'palevioletred', 
                   'mediumpurple', 
                   'cornsilk4', 
                   'olivedrab3', 
                   'steelblue4', 
                   'indianred4',
                   'aquamarine3',
                   'saddlebrown',
                   'deeppink2',
                   'cornflowerblue',
                   'black')
names(colors_cluster) = paste0('C',0:15)
tail = c('Tail' = 'gray60')
colors_cluster = c(colors_cluster, tail)


# univariate
tool = 'mobster_univariate'
mutations <-  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal/tables_interpreted/', tool, '_', spn, '_', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '.rds')) %>% 
  separate(cluster_id_process, sep = '_', into = c('cluster_id_process', 'tmp')) %>% 
  select(-tmp)
samples <- paste(spn, get_sample_names(spn), sep= '_')
mutations <- mutations %>% mutate(driver_label = sub("_.*", "", driver_label_tool))

driver_tool <- mutations %>% select(patient_id, mutation_id, is_driver_tool, driver_label, is_driver_process)

id_muts <- mutations %>% 
  dplyr::summarise(n = dplyr::n(), .by = c(patient_id, mutation_id, cluster_id_tool,
                                           sample_id)) |>
  dplyr::filter(n > 1L) 

muts_tool_raw <- mutations %>% 
  filter(!mutation_id %in% id_muts$mutation_id) %>% 
  select(patient_id, mutation_id, cluster_id_tool, sample_id, vaf_tool) %>% 
  pivot_wider(values_from = vaf_tool, names_from = sample_id) %>% 
  replace(is.na(.), 0) %>% 
  mutate(cluster_id_tool = as.character(cluster_id_tool))
muts_tool_raw <- muts_tool_raw %>% left_join(driver_tool) %>% distinct()

muts_process <- mutations %>% 
  filter(!mutation_id %in% id_muts$mutation_id) %>% 
  select(patient_id, mutation_id, sample_id, vaf_process, cluster_id_process) %>%
  pivot_wider(values_from = vaf_process, names_from = sample_id) %>% 
  replace(is.na(.), 0) 

driver_process <- mutations %>%   filter(!mutation_id %in% id_muts$mutation_id) %>% select(patient_id, mutation_id, is_driver_process, driver_label)
muts_process <- muts_process %>% left_join(driver_process) %>% distinct()

cluster_plots_samples <- list()
for (i in 1:length(samples)){
  s1 <- samples[[i]]
  if (s1 %in% colnames(muts_tool_raw)){
    
    color_palette_process = RColorBrewer::brewer.pal(n = max(3,length(unique(muts_process$cluster_id_process))), name = "Dark2") %>%
      setNames(str_sort(unique(muts_process$cluster_id_process), numeric=T))
    color_palette_process['Subclonal'] = 'gray70'
    
    p_process = muts_process %>%
      filter(.data[[s1]] != 0) %>% 
      ggplot(aes(x =.data[[s1]], fill = cluster_id_process))+
      geom_histogram( alpha=0.6, binwidth = 0.01, position = "identity")+
      xlim(-0.01,1.01)+
      xlab(paste0('VAF ', s1)) +
      theme_bw() +
      scale_fill_manual('ProCESS clusters', values = color_palette_process)+
      ggtitle('ProCESS') +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))  +
      ylab('')
    
    p_process = p_process + ggrepel::geom_label_repel(
      data = muts_process %>% filter(is_driver_process == TRUE) %>% filter(.data[[s1]] != 0),
      aes(
        x = .data[[s1]],
        y = 10,
        label = driver_label,
        colour = cluster_id_process, 
      ),
      show.legend = F,
      inherit.aes = FALSE,
      size = 3,
      min.segment.length = 0,
      box.padding = 1) +
      scale_color_manual('ProCESS clusters', values = color_palette_process)
    
  
    p_tool_raw = muts_tool_raw %>%
      filter(.data[[s1]] != 0) %>% 
      ggplot(aes(x =.data[[s1]], fill=cluster_id_tool))+
      geom_histogram( alpha=0.6, binwidth = 0.01, position = "identity")+
      xlim(-0.01,1.01)+
      xlab(paste0('VAF ', s1)) +
      theme_bw() +
      ylab('') + 
      scale_fill_manual('MOBSTER clusters', values = colors_cluster)+
      theme_bw() +
      ggtitle('MOBSTER') +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
    
    
    p_tool_raw = p_tool_raw + ggrepel::geom_label_repel(
      data = muts_tool_raw %>% filter(is_driver_tool == TRUE) %>% filter(.data[[s1]] != 0),
      aes(
        x = .data[[s1]],
        y = 10,
        label = driver_label,
        colour = cluster_id_tool, 
      ),
      show.legend = F,
      inherit.aes = FALSE,
      size = 3,
      min.segment.length = 0,
      box.padding = 1,
      max.overlaps = 50) +
      scale_color_manual('MOBSTER clusters', values = colors_cluster)
      
    cluster_plots_samples[[s1]] <- p_process + p_tool_raw + plot_layout(guides = 'collect')
  }
}

heigh <- rep(1, length(samples))
heigh <- c(heigh, 2)

# multivariate
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


tool = opt$tool
plot_sign_all <- list()
plot_sign <- list()

sign_type='SBS'
for (sign_type in c('SBS', 'ID')){
  if (sign_type == 'SBS'){
    s_type = 'SBS96'
    colors = sbs_colors
  } else {
    s_type = 'ID83'
    colors = id_colors
  }
  
  
  process = readRDS(paste0(indir, spn, '/process/',cov, 'x_', pur, 'p/exposure_', sign_type, '.rds'))
  process_df = process %>% 
    select(cluster_id_process, causes, nmuts_cause, tot_nmuts, exposure) %>% 
    dplyr::rename(Signature=causes,
                  Exposure=exposure,
                  Samples=cluster_id_process)
  
  process_plot <- plot_exposure(process_df, sig_type = sign_type, tool = 'ProCESS', col = colors, type = 'ProCESS')
  
  type = 'raw'
  if (signature_tool == 'BASCULE'){
    name_file = paste0('bascule_fit_',type,'.rds')
  } else{
    name_file = paste0('Assignment_Solution/Activities/Assignment_Solution_Activities.txt')
  }
  
  if (signature_tool == 'BASCULE'){
    
    
    data = readRDS(paste0(indir, spn, '/', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '/', tool, '_', signature_tool, '_', type, '/', s_type, '/',name_file))
    
    df = data[['nmf']][[sign_type]][['exposure']] %>% 
      dplyr::rename(Signature=sigs,
                    Exposure=value,
                    Samples=samples)
    if (tool == 'pyclonevi'){
      df = df %>% mutate(Samples = str_replace(Samples, "X", "C"))
    }
  } else {
    data = read.table(file = paste0(indir, spn, '/', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '/', tool, '_', signature_tool, '_', type, '/', s_type, '/',name_file),header = T) 
    df = data %>% pivot_longer(
      cols = starts_with(sign_type),
      names_to = "Signature",
      values_to = "Exposure")
    
  }

  tool_plot <- plot_exposure(df, sig_type = sign_type, tool = tool, type = tool, signature = signature_tool, col = colors)
  plot_sign_all[[sign_type]]  <- wrap_plots(process_plot,tool_plot) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
}

mutations <-  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal/tables_interpreted/', tool, '_', spn, '_', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '.rds'))
samples <- unique(mutations$sample_id)
pairs <- combn(samples, 2, simplify = FALSE)
mutations <- mutations %>% mutate(driver_label = sub("_.*", "", driver_label_tool))


if (spn == 'SPN07'){
  driver_tool <- mutations %>% select(patient_id, mutation_id, is_driver_tool, driver_label)
  dups_muts <- mutations %>% dplyr::summarise(n = dplyr::n(), .by = c(patient_id, mutation_id, cluster_id_tool, sample_id)) %>% dplyr::filter(n > 1L) 
  muts_tool_raw <- mutations %>% 
    select(patient_id, mutation_id, cluster_id_tool, sample_id, vaf_tool) %>% 
    filter(!(mutation_id %in% dups_muts$mutation_id)) %>% 
    pivot_wider(values_from = vaf_tool, names_from = sample_id) %>% 
    replace(is.na(.), 0) %>% 
    mutate(cluster_id_tool = as.character(cluster_id_tool))
  muts_tool_raw <- muts_tool_raw %>% left_join(driver_tool) %>% distinct()
  
  muts_process <- mutations %>% 
    select(patient_id, mutation_id, sample_id, vaf_process, cluster_id_process) %>%
    filter(!(mutation_id %in% dups_muts$mutation_id)) %>% 
    pivot_wider(values_from = vaf_process, names_from = sample_id) %>% 
    replace(is.na(.), 0) 
  
  driver_process <- mutations %>% select(patient_id, mutation_id, is_driver_process, driver_label)
  muts_process <- muts_process %>% left_join(driver_process) %>% distinct()
  
} else {
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
}

if (tool == 'pyclonevi'){
  muts_tool_raw$cluster_id_tool = paste0('C', muts_tool_raw$cluster_id_tool)
}

cluster_plots_pair <- list()

for (i in 1:length(pairs)){
  s1 <- pairs[[i]][1]
  s2 <- pairs[[i]][2]
  
  color_palette_process = RColorBrewer::brewer.pal(n = max(3,length(unique(muts_process$cluster_id_process))), name = "Dark2") %>%
    setNames(str_sort(unique(muts_process$cluster_id_process), numeric=T))
  color_palette_process['Subclonal'] = 'gainsboro'
  
  p_process = muts_process %>%
    ggplot(aes(x =.data[[s1]], y = .data[[s2]], color=cluster_id_process))+
    geom_point( alpha=0.2, size = .5)+
    xlim(0,1)+
    ylim(0,1)+
    xlab(pairs[[i]][1]) +
    ylab(pairs[[i]][2])+
    theme_bw() +
    scale_color_manual('ProCESS clusters', values = color_palette_process)+
    ggtitle('ProCESS') +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1) ) ) 
  
  p_process = p_process + ggrepel::geom_label_repel(
    data = muts_process %>% filter(is_driver_process == TRUE),
    aes(
      x = .data[[s1]],
      y = .data[[s2]],
      label = driver_label,
      colour = cluster_id_process, 
    ),
    show.legend = F,
    inherit.aes = FALSE,
    size = 3,
    min.segment.length = 0,
    box.padding = 1)
  
  
  p_tool = muts_tool_raw %>%
      ggplot(aes(x =.data[[s1]], y = .data[[s2]], color=cluster_id_tool))+
      geom_point( alpha=0.2, size = .5)+
      xlim(0,1)+
      ylim(0,1)+
      xlab(pairs[[i]][1]) +
      ylab(pairs[[i]][2])+
      scale_color_manual('Tool clusters', values = colors_cluster)+
      theme_bw() +
      ggtitle(tool) +
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
        size = 3,
        min.segment.length = 0,
        box.padding = 1,
        max.overlaps = 50)
      
  cluster_plots_pair[[i]] <- p_process + p_tool
}

samples <- names(cluster_plots_samples)
heigh_mobster <- rep(.5, length(samples))
heigh <- c(heigh_mobster,length(pairs),1.5)
p_final <- wrap_plots(cluster_plots_samples, 
                       ncol = 2, 
                       guides = 'collect') +  
  wrap_plots(cluster_plots_pair, 
             nrow = length(pairs), 
             guides = 'collect')   + 
  wrap_plots(plot_sign_all$SBS,
             plot_sign_all$ID, 
             nrow = 2,  
             guides = 'collect') + 
  plot_layout(ncol = 1, 
              heights = heigh)

name_file = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/plot/', spn, '/',spn,'_',cov, 'x_', pur, 'p_', tool,'_',signature_tool,'.png')
ggsave(plot = p_final, filename = name_file, dpi = 400, width = 7, height = 3*(length(pairs)+2)+(length(samples)+2))
