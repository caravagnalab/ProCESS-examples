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
                    make_option(c("--tool"), type = "character", default = 'pyclonevi')
)


compare_exposure <- function(process_df, tool_df,
                             sig_col = "Signature",
                             exp_col = "Exposure") {
  
  sig_col <- rlang::sym(sig_col)
  exp_col <- rlang::sym(exp_col)
  
  # Aggregate in case multiple rows per signature exist
  proc <- process_df %>%
    group_by(!!sig_col) %>%
    summarise(exp_proc = sum(!!exp_col), .groups = "drop")
  
  tool <- tool_df %>%
    group_by(!!sig_col) %>%
    summarise(exp_tool = sum(!!exp_col), .groups = "drop")
  
  # Align signatures
  merged <- full_join(proc, tool, by = rlang::as_string(sig_col)) %>%
    mutate(
      exp_proc = replace_na(exp_proc, 0),
      exp_tool = replace_na(exp_tool, 0)
    ) %>%
    arrange(!!sig_col)
  
  p <- merged$exp_proc
  t <- merged$exp_tool
  
  # ---- Metrics ----
  rmse <- sqrt(mean((p - t)^2))
  
  cosine_sim <- if (sum(p^2) == 0 || sum(t^2) == 0) {
    NA_real_
  } else {
    sum(p * t) / (sqrt(sum(p^2)) * sqrt(sum(t^2)))
  }
  
  tibble(
    RMSE = rmse,
    CosineSimilarity = cosine_sim
  )
}

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


plot_sign_all <- list()
plot_sign <- list()
df_comparison <- list()

tool=opt$tool
sign_type = 'SBS'

for (sign_type in c('SBS', 'ID')){
  print(sign_type)
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
  
  type = 'int'
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
    df = data %>%
      rowwise() %>%
      mutate(
        total = sum(c_across(-Samples)),
        across(-c(Samples, total), ~ ifelse(total > 0, .x / total, 0))
      ) %>%
      ungroup() %>%
      select(-total) %>%
      pivot_longer(
        cols = -Samples,
        names_to = "Signature",
        values_to = "Exposure"
      )
  }
  
  
  if (type == 'int'){
    sp = spn
    t = tool
    result_int <- readRDS('/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/metrics_tables/metrics_drivers_clonal_vs_subclonal.rds') %>% 
      filter(purity == pur, 
             coverage == cov, 
             spn == sp, 
             tool == t)
    tp_cluster = result_int$TP_c_list %>% unlist()
    if (signature_tool == 'BASCULE' & tool == 'pyclonevi'){
      df_tp = df %>% filter(Samples %in% paste0('X',tp_cluster))
    } else{
      df_tp = df %>% filter(Samples %in% tp_cluster)
    }
    
    mutations <-  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal/tables_interpreted/', tool, '_', spn, '_', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '.rds'))
    annot_table <- mutations %>% 
      filter(cluster_id_tool %in% tp_cluster, 
             is_driver_tool == T | is_driver_process == T) %>% 
      select(cluster_id_tool, driver_label_tool, driver_label_process, cluster_id_process,is_driver_tool, is_driver_process) %>% 
      unique() %>% 
      tidyr::separate(driver_label_tool, into = c('driver_tool', 'tmp_tool'), sep = '_') %>% 
      tidyr::separate(driver_label_process, into = c('driver_process', 'tmp_process'), sep = ' ')
    
    
    tool_cluster = annot_table %>% filter(driver_tool %in% annot_table$driver_process) %>% filter(is_driver_tool == T) %>% select(cluster_id_tool, driver_tool, driver_process)
    process_cluster = annot_table %>% filter(driver_process %in% tool_cluster$driver_process) %>% filter(is_driver_process == T) %>% select(cluster_id_process, driver_process)
    
    exp_process = process_df %>% filter(Samples %in% process_cluster$cluster_id_process) %>% left_join(process_cluster %>% dplyr::rename(Samples = cluster_id_process))
    if (signature_tool == 'BASCULE' & tool == 'pyclonevi'){
      exp_tool = df %>% filter(Samples %in% paste0('X', tool_cluster$cluster_id_tool)) %>% left_join(tool_cluster %>% mutate(cluster_id_tool = paste0('X', as.character(cluster_id_tool))) %>%  dplyr::rename(Samples = cluster_id_tool))
    } else{
      exp_tool = df %>% filter(Samples %in%  tool_cluster$cluster_id_tool) %>% left_join(tool_cluster %>% mutate(cluster_id_tool = as.character(cluster_id_tool)) %>%  dplyr::rename(Samples = cluster_id_tool))
    }
    
    df_comparison[[sign_type]]  <- lapply(unique(process_cluster$driver_process), FUN = function(c){
      tmp_process = exp_process %>% filter(driver_process == c)
      tmp_tool = exp_tool %>% filter(driver_tool == c)
      compare_exposure(process_df = tmp_process, tool_df = tmp_tool) %>% mutate(cluster = c,
                                                                                cluster_process = unique(tmp_process$Samples))
    }) %>% bind_rows() %>% mutate(type = sign_type)
    
    
    color_palette_driver =   RColorBrewer::brewer.pal(n = max(3,length(unique(exp_process$driver_process))), name = "Dark2") %>%
      setNames(str_sort(unique(exp_process$driver_process), numeric=T))
    color_palette_driver['Subclonal'] = 'gainsboro'
    process_plot <- plot_exposure(exp_process %>% select(-Samples) %>% dplyr::rename(Samples = driver_process), sig_type = sign_type, tool = 'ProCESS', col = colors)
    process_stat <- plot_stat(exp_process %>% select(-Samples) %>% dplyr::rename(Samples = driver_process), sig_type = sign_type, tool = 'ProCESS', col = color_palette_driver)
    
    plot_process <- process_plot + process_stat + plot_layout(design = 'aaab')
    
    
    color_palette_driver =   RColorBrewer::brewer.pal(n = max(3,length(unique(exp_process$driver_process))), name = "Dark2") %>%
      setNames(str_sort(unique(exp_process$driver_process), numeric=T))
    color_palette_driver['Subclonal'] = 'gainsboro'
    tool_plot <- plot_exposure(exp_tool %>% select(-Samples) %>% dplyr::rename(Samples = driver_tool), sig_type = sign_type, tool = 'Tool', col = colors)
    tool_stat <- plot_stat(exp_tool %>% select(-Samples) %>% dplyr::rename(Samples = driver_tool), sig_type = sign_type, tool = 'Tool', col = color_palette_driver)
    
    plot_process <- process_plot + process_stat + plot_layout(design = 'aaab')
    plot_tool <- tool_plot + tool_stat + plot_layout(design = 'aaab')
    
    
    plot_sign_all[[sign_type]]  <- wrap_plots(plot_process,plot_tool) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
  }
}


mutations <-  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal/tables_interpreted/', tool, '_', spn, '_', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '.rds'))
samples <- unique(mutations$sample_id)
pairs <- combn(samples, 2, simplify = FALSE)

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

cluster_plots_pair <- list()

for (i in 1:length(pairs)){
  s1 <- pairs[[i]][1]
  s2 <- pairs[[i]][2]
  
  p_process = mutations_process %>%
    ggplot(aes(x =.data[[s1]], y = .data[[s2]], color=driver_process))+
    geom_point( alpha=0.2, size = .5)+
    xlim(0,1)+
    ylim(0,1)+
    xlab(pairs[[i]][1]) +
    ylab(pairs[[i]][2])+
    theme_bw() +
    scale_color_manual('ProCESS clusters', values = color_palette_driver)+
    ggtitle('ProCESS') +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1) ) ) 
  
  p_process = p_process + ggrepel::geom_label_repel(
    data = mutations_process %>% filter(is_driver_process == TRUE),
    aes(
      x = .data[[s1]],
      y = .data[[s2]],
      label = driver_process,
      colour = driver_process, 
    ),
    show.legend = F,
    inherit.aes = FALSE,
    size = 3,
    min.segment.length = 0,
    box.padding = 1)
  
  
  p_tool = mutations_tool %>%
    ggplot(aes(x =.data[[s1]], y = .data[[s2]], color=driver_tool))+
    geom_point( alpha=0.2, size = .5)+
    xlim(0,1)+
    ylim(0,1)+
    xlab(pairs[[i]][1]) +
    ylab(pairs[[i]][2])+
    theme_bw() +
    scale_color_manual(paste0(tool, ' clusters'), values = color_palette_driver)+
    ggtitle(tool) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1) ) ) 
  
  p_tool = p_tool + ggrepel::geom_label_repel(
    data = mutations_tool %>% filter(is_driver_tool == TRUE, cluster_id_tool != 'NA') %>% 
      filter(str_detect(driver_label_tool, driver_tool)),
    aes(
      x = .data[[s1]],
      y = .data[[s2]],
      label = driver_tool,
      colour = driver_tool, 
    ),
    show.legend = F,
    inherit.aes = FALSE,
    size = 3,
    min.segment.length = 0,
    box.padding = 1)
  
  cluster_plots_pair[[i]] <- p_process + p_tool
}

heigh <- rep(1, length(pairs))
heigh <- c(heigh,2)

p_final <- wrap_plots(cluster_plots_pair, nrow = length(pairs)+1, guides = 'collect') + 
  wrap_plots(plot_sign_all$SBS,plot_sign_all$ID, nrow = 2,  guides = 'collect') + 
  plot_layout(heights = heigh)
ggsave(plot = p_final, filename = paste0('assign_res/', spn, '/',spn,'_',cov, 'x_', pur, 'p_',tool,'_',signature_tool,'.png'), dpi = 400, width = 7.5, height = 2.5*(length(pairs)+2))
#ggsave(plot = p_final, filename = paste0(spn,'_',tool,'_',signature_tool,'_res.pdf'), width = 9.5, height = 2.5*(length(pairs)+2))


results <- df_comparison %>% 
  bind_rows()  %>% 
  mutate(cov = cov, pur = pur, tool = tool, sig_tool = signature_tool, cna_caller = cna_caller, mut_caller = mut_caller, spn = spn)
saveRDS(object = results,file = paste0('assign_res/',spn, '/',spn,'_',cov, 'x_', pur, 'p_',tool,'_',signature_tool,'.rds'))


