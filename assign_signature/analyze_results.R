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

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN07'),
                    make_option(c("--purity"), type = "double", default = 0.9),
                    make_option(c("--coverage"), type = "integer", default = 50),
                    make_option(c("--cna_caller"), type = "character", default = 'ascat'),
                    make_option(c("--vcf_caller"), type = "character", default = 'mutect2'),
                    make_option(c("--signature"), type = "character", default = 'BASCULE'),
                    make_option(c("--tool"), type = "character", default = 'viber')
)

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

plot_stat <- function(df, sig_type, tool, col, type = '', signature = ''){
  if (tool == 'ProCESS'){
    plot <- df %>% 
      select(-Exposure, -nmuts_cause, -Signature) %>% 
      distinct() %>% 
      ggplot(aes(y=Samples, x=tot_nmuts, fill = Samples)) + 
      geom_col(show.legend = F)+
      scale_fill_manual(values=col)+
      theme_bw()+
      theme(legend.position = "bottom")+
      scale_x_continuous(labels = scales::label_scientific()) +
      ylab("")+
      xlab("Nmuts") +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  } else{
    plot <- df %>% 
      group_by(Samples) %>% 
      summarise(tot_nmuts = sum(Exposure)) %>% 
      ggplot(aes(y=Samples, x=tot_nmuts, fill = Samples)) + 
      geom_col(show.legend = F)+
      scale_fill_manual(values=col)+
      scale_x_continuous(labels = scales::label_scientific()) +
      theme_bw()+
      theme(legend.position = "bottom")+
      ylab("")+
      xlab("N muts") +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  return(plot)
}

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

spn = opt$spn_id
cov = opt$coverage
pur = opt$purity
signature_tool =  opt$signature
cna_caller = opt$cna_caller
mut_caller = opt$vcf_caller

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
  
  color_palette_process = RColorBrewer::brewer.pal(n = max(3,length(unique(process_df$Samples))), name = "Dark2") %>%
    setNames(str_sort(unique(process_df$Samples), numeric=T))
  color_palette_process['Subclonal'] = 'gainsboro'
  
  process_plot <- plot_exposure(process_df, sig_type = sign_type, tool = 'ProCESS', col = colors)
  process_stat <- plot_stat(process_df, sig_type = sign_type, tool = 'ProCESS', col = color_palette_process)
  
  plot_process <- process_plot + process_stat + plot_layout(design = 'aaab')
  
  for (type in c('raw', 'int')){

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
    } else {
      data = read.table(file = paste0(indir, spn, '/', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '/', tool, '_', signature_tool, '_', type, '/', s_type, '/',name_file),header = T) 
      df = data %>% pivot_longer(
        cols = starts_with(sign_type),
        names_to = "Signature",
        values_to = "Exposure")
    }
    
    if (length(unique(df$Samples)) <= 12){
      color_palette_tool = RColorBrewer::brewer.pal(n = max(3,length(unique(df$Samples))), name = "Paired") %>%
        setNames(str_sort(unique(df$Samples), numeric=T))
      color_palette_tool['Subclonal'] = 'gainsboro'
    }else{
      miss = length(unique(df$Samples)) - 12
      color_palette_tool = c(RColorBrewer::brewer.pal(n = 12, name = "Paired"), RColorBrewer::brewer.pal(n = miss, name = "Dark2"))%>%
        setNames(str_sort(unique(df$Samples), numeric=T))
      color_palette_tool['Subclonal'] = 'gainsboro'
    }

    
    tool_plot <- plot_exposure(df, sig_type = sign_type, tool = tool, type = type, signature = signature_tool, col = colors)
    tool_stat <- plot_stat(df %>% mutate(Samples = as.character(Samples)), sig_type = sign_type, tool = tool, type = type, signature = signature_tool, col = color_palette_tool)
  
    plot_sign[[type]] <- tool_plot + tool_stat + plot_layout(design = 'aaab')
  }
  
  plot_sign_all[[sign_type]]  <- wrap_plots(plot_process,plot_sign$raw,plot_sign$int) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
}

plot_sign_all$SBS
# sp = spn
# t = tool
# metrics <- readRDS('/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/metrics_tables/table_clusters_metrics.rds') %>% 
#   filter(spn == sp & coverage == cov & purity == pur & vcf_caller == mut_caller & cna_caller == cna_caller & tool == t) %>% 
#   distinct()

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
  
  muts_tool_int <- mutations %>% 
    select(patient_id, mutation_id, cluster_id_tool_interpreted, sample_id, vaf_tool) %>% 
    filter(!(mutation_id %in% dups_muts$mutation_id)) %>% 
    pivot_wider(values_from = vaf_tool, names_from = sample_id) %>% 
    replace(is.na(.), 0) %>% 
    mutate(cluster_id_tool_interpreted = as.character(cluster_id_tool_interpreted))
  muts_tool_int <- muts_tool_int %>% left_join(driver_tool) %>% distinct()
  
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
  
  muts_tool_int <- mutations %>% 
    select(patient_id, mutation_id, cluster_id_tool_interpreted, sample_id, vaf_tool) %>% 
    pivot_wider(values_from = vaf_tool, names_from = sample_id) %>% 
    replace(is.na(.), 0) %>% 
    mutate(cluster_id_tool_interpreted = as.character(cluster_id_tool_interpreted))
  muts_tool_int <- muts_tool_int %>% left_join(driver_tool) %>% distinct()
  
  muts_tool_int_driver <- mutations %>% 
    select(patient_id, mutation_id, cluster_id_tool_interpreted_driver, sample_id, vaf_tool) %>% 
    pivot_wider(values_from = vaf_tool, names_from = sample_id) %>% 
    replace(is.na(.), 0) %>% 
    mutate(cluster_id_tool = as.character(cluster_id_tool_interpreted_driver))
  muts_tool_int_driver <- muts_tool_int_driver %>% left_join(driver_tool) %>% distinct()
  
  muts_process <- mutations %>% 
    select(patient_id, mutation_id, sample_id, vaf_process, cluster_id_process) %>%
    pivot_wider(values_from = vaf_process, names_from = sample_id) %>% 
    replace(is.na(.), 0) 
  
  driver_process <- mutations %>% select(patient_id, mutation_id, is_driver_process, driver_label)
  muts_process <- muts_process %>% left_join(driver_process) %>% distinct()
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
  
  for (type in c('raw', 'int', 'driver')){
    if (type == 'raw'){
      if (length(unique(muts_tool_raw$cluster_id_tool)) <= 12){
        color_palette_tool = RColorBrewer::brewer.pal(n = max(3,length(unique(muts_tool_raw$cluster_id_tool))), name = "Paired") %>%
          setNames(str_sort(unique(muts_tool_raw$cluster_id_tool), numeric=T))
        color_palette_tool['Subclonal'] = 'gainsboro'
      }else{
        miss = length(unique(muts_tool_raw$cluster_id_tool)) - 12
        color_palette_tool = c(RColorBrewer::brewer.pal(n = 12, name = "Paired"), RColorBrewer::brewer.pal(n = miss, name = "Dark2"))%>%
          setNames(str_sort(unique(muts_tool_raw$cluster_id_tool), numeric=T))
        color_palette_tool['Subclonal'] = 'gainsboro'
      }
      
      p_tool_raw = muts_tool_raw %>%
        ggplot(aes(x =.data[[s1]], y = .data[[s2]], color=cluster_id_tool))+
        geom_point( alpha=0.2, size = .5)+
        xlim(0,1)+
        ylim(0,1)+
        xlab(pairs[[i]][1]) +
        ylab(pairs[[i]][2])+
        scale_color_manual('Tool clusters Raw', values = color_palette_tool)+
        theme_bw() +
        ggtitle('Raw') +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1) ) )
      
    } else{
      if (length(unique(muts_tool_int$cluster_id_tool_interpreted)) <= 12){
        color_palette_tool = RColorBrewer::brewer.pal(n = max(3,length(unique(muts_tool_int$cluster_id_tool_interpreted))), name = "Paired") %>%
          setNames(str_sort(unique(muts_tool_int$cluster_id_tool_interpreted), numeric=T))
        color_palette_tool['Subclonal'] = 'gainsboro'
      }else{
        miss = length(unique(muts_tool_int$cluster_id_tool_interpreted)) - 12
        color_palette_tool = c(RColorBrewer::brewer.pal(n = 12, name = "Paired"), RColorBrewer::brewer.pal(n = miss, name = "Dark2"))%>%
          setNames(str_sort(unique(muts_tool_int$cluster_id_tool_interpreted), numeric=T))
        color_palette_tool['Subclonal'] = 'gainsboro'
      }

      p_tool_int = muts_tool_int %>%
        ggplot(aes(x =.data[[s1]], y = .data[[s2]], color=cluster_id_tool_interpreted), alpha=0.2, size = .1)+
        geom_point(alpha=0.2, size = .5)+
        xlim(0,1)+
        ylim(0,1)+
        xlab(pairs[[i]][1]) +
        ylab(pairs[[i]][2])+
        scale_color_manual('Tool clusters Int', values = color_palette_tool)+
        theme_bw() +
        ggtitle('Interpreted') +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1) ) ) 
      
      p_tool_int = p_tool_int + ggrepel::geom_label_repel(
        data = muts_tool_int %>% filter(is_driver_tool == TRUE),
        aes(
          x = .data[[s1]],
          y = .data[[s2]],
          label = driver_label,
          colour = cluster_id_tool_interpreted, 
        ),
        show.legend = F,
        inherit.aes = FALSE,
        size = 3,
        min.segment.length = 0,
        box.padding = 1,
        max.overlaps = 50)
      
    }
  }
  
  cluster_plots_pair[[i]] <- p_process + p_tool_raw + p_tool_int
}

heigh <- rep(1, length(pairs))
heigh <- c(heigh,2)
p_final <- wrap_plots(cluster_plots_pair, nrow = length(pairs)+1, guides = 'collect')   + wrap_plots(plot_sign_all$SBS,plot_sign_all$ID, nrow = 2,  guides = 'collect') + plot_layout(heights = heigh)
ggsave(plot = p_final, filename = paste0('assign_plot/', spn, '/',spn,'_',cov, 'x_', pur, 'p_', tool,'_',signature_tool,'.png'), dpi = 400, width = 9.5, height = 2.5*(length(pairs)+2))
#ggsave(plot = p_final, filename = paste0(spn,'_',cov, 'x_', pur, 'p_', tool,'_',signature_tool,'.pdf'), width = 9.5, height = 2.5*(length(pairs)+2))

pp <- wrap_plots(cluster_plots_pair, nrow = length(pairs)+1, guides = 'collect')
ggsave(plot = pp, filename = paste0('assign_plot/',spn, '/', spn,'_',cov, 'x_', pur, 'p_', tool,'_multivariate.png'), dpi = 400, width = 10, height = 2.5*(length(pairs)+2))

