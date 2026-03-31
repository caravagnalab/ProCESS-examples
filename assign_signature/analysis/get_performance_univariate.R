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
                    make_option(c("--signature"), type = "character", default = 'BASCULE')
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


tool = 'mobster_univariate'
plot_sign_all <- list()
samples <- get_sample_names(spn)
df_comparison <- list()
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
    
    for (s in samples){
    
      process = readRDS(paste0(indir, spn, '/process_univariate/', cov, 'x_', pur, 'p/', s, '_exposure_', sign_type, '.rds'))
      process_df = process %>% 
        select(sample_id, cluster_id_process, causes, nmuts_cause, tot_nmuts, exposure) %>% 
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
        
        
        data = readRDS(paste0(indir, spn, '/', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '/', s, '/', tool, '_', signature_tool, '_', type, '/', s_type, '/',name_file))
        
        df = data[['nmf']][[sign_type]][['exposure']] %>% 
          dplyr::rename(Signature=sigs,
                        Exposure=value,
                        Samples=samples)
        
        counts = data[['input']][[sign_type]][['counts']] %>% 
          group_by(samples) %>% 
          summarise(n = sum(value)) %>% 
          dplyr::rename(Samples = samples)
      } else {
        data = read.table(file = paste0(indir, spn, '/', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '/', s, '/', tool, '_', signature_tool, '_', type, '/', s_type, '/',name_file),header = T) 
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
        
        counts = data %>%  mutate(n = rowSums(across(-Samples))) %>% select(Samples, n)
      }
      
      
      if (type == 'int'){
        sp = spn
        t = tool
        result_int <- readRDS('/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/metrics_tables/metrics_drivers_clonal_vs_subclonal.rds') %>% 
          filter(purity == pur, 
                 coverage == cov, 
                 spn == sp, 
                 tool == 'mobster') %>% 
          filter(sample == paste0(spn, '_', s))
        tp_cluster = result_int$TP_c_list %>% unlist() %>% unique()
        df_tp = df %>% filter(Samples %in% tp_cluster)
        
        mutations <-  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal/tables_interpreted/', tool, '_', spn, '_', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '.rds'))
        annot_table <- mutations %>% 
          filter(sample_id == paste0(spn, '_', s)) %>% 
          filter(cluster_id_tool %in% tp_cluster, 
                 is_driver_tool == T | is_driver_process == T) %>% 
          select(cluster_id_tool, driver_label_tool, driver_label_process, cluster_id_process,is_driver_tool, is_driver_process) %>% 
          unique() %>% 
          tidyr::separate(driver_label_tool, into = c('driver_tool', 'tmp_tool'), sep = '_') %>% 
          tidyr::separate(driver_label_process, into = c('driver_process', 'tmp_process'), sep = ' ')
        
        
        tool_cluster = annot_table %>% filter(driver_tool %in% annot_table$driver_process) %>% filter(is_driver_tool == T) %>% select(cluster_id_tool, driver_tool, driver_process)
        process_cluster = annot_table %>% filter(driver_process %in% tool_cluster$driver_process) %>% filter(is_driver_process == T) %>% select(cluster_id_process, driver_process)
        
        exp_process = process_df %>% filter(Samples %in% process_cluster$cluster_id_process) %>% left_join(process_cluster %>% dplyr::rename(Samples = cluster_id_process))
        exp_tool = df %>% filter(Samples %in%  tool_cluster$cluster_id_tool) %>% left_join(tool_cluster %>% mutate(cluster_id_tool = as.character(cluster_id_tool)) %>%  dplyr::rename(Samples = cluster_id_tool))

        
        df_comparison[[s]][[sign_type]]  <- lapply(unique(process_cluster$driver_process), FUN = function(c){
          tmp_process = exp_process %>% filter(driver_process == c)
          tmp_tool = exp_tool %>% filter(driver_tool == c) %>% left_join(counts) %>% filter(!is.na(driver_process))
          compare_exposure(process_df = tmp_process, tool_df = tmp_tool) %>% mutate(cluster = c,
                                                                                    cluster_process = unique(tmp_process$Samples),
                                                                                    nmuts = unique(tmp_tool$n))
        }) %>% bind_rows() %>% mutate(type = sign_type, sample = s)
        
        # analyze subclonal vs tail 
        tmp_process = process_df %>% filter(Samples == 'Subclonal')
        tmp_tool = df %>% filter(Samples == 'Tail')
        compare = compare_exposure(process_df = tmp_process, tool_df = tmp_tool) %>% 
          mutate(cluster = 'Tail-Subclonal',
                 cluster_process = 'Subclonal', 
                 type = sign_type, 
                 sample = s)
        
        df_comparison[[s]][[sign_type]] <- df_comparison[[s]][[sign_type]] %>% bind_rows(compare)
        
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
        
        plot_sign_all[[s]][[sign_type]]  <- wrap_plots(plot_process,plot_tool) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
    }
  }
}



sp = spn
t = tool
result_int <- readRDS('/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/metrics_tables/metrics_drivers_clonal_vs_subclonal.rds') %>% 
  filter(purity == pur, 
         coverage == cov, 
         spn == sp, 
         tool == 'mobster')
tp_cluster = result_int$TP_c_list %>% unlist() %>% unique()

mutations <-  readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT//validation_subclonal/tables_interpreted/', tool, '_', spn, '_', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '.rds'))

jaccard_df <- tibble()
for (s in samples){
  annot_table <- mutations %>% 
    filter(sample_id == paste0(spn, '_', s)) %>% 
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
  
  
  sample = paste(spn, s, sep = '_')
  mutations_process <- mutations %>% 
    mutate(cluster_id_process = ifelse(cluster_id_process %in% process_cluster$cluster_id_process, cluster_id_process, 'NA')) %>% 
    select(patient_id, mutation_id, sample_id, cluster_id_process, vaf_process, is_driver_process) %>% 
    left_join(process_cluster) %>% 
    replace(is.na(.), 'Other')  %>% 
    filter(!mutation_id %in% dups$mutation_id) %>% 
    pivot_wider(values_from = vaf_process, names_from = sample_id) %>% 
    replace(is.na(.), 0) %>% 
    filter(.data[[sample]] != 0)
  
  mutations_tool <- mutations %>% 
    mutate(cluster_id_tool_interpreted = ifelse(cluster_id_tool_interpreted %in% tool_cluster$cluster_id_tool, cluster_id_tool_interpreted, 'NA')) %>% 
    select(patient_id, mutation_id, sample_id, cluster_id_tool_interpreted, vaf_tool, is_driver_tool, driver_label_tool) %>% 
    dplyr::rename(cluster_id_tool = cluster_id_tool_interpreted) %>% 
    left_join(tool_cluster %>% mutate(cluster_id_tool = as.character(cluster_id_tool))) %>% 
    replace(is.na(.), 'Other')  %>% 
    filter(!mutation_id %in% dups$mutation_id) %>% 
    pivot_wider(values_from = vaf_tool, names_from = sample_id) %>% 
    replace(is.na(.), 0)  %>% 
    filter(.data[[sample]] != 0)
  
  join = tool_cluster %>% left_join(process_cluster) %>% filter(!is.na(cluster_id_process))
  
  tp_tool_muts <- mutations_tool %>% filter(cluster_id_tool %in% unique(tool_cluster$cluster_id_tool))
  tp_process_muts <- mutations_process  %>% filter(cluster_id_process %in% unique(process_cluster$cluster_id_process))
  jaccard <- function(a, b) {
    length(intersect(a, b)) / length(union(a, b))
  }
  
  for (d in unique(join$driver_tool)){
    p = join %>% filter(driver_tool == d) %>% pull(cluster_id_process)
    t = join %>% filter(driver_tool == d)%>% pull(cluster_id_tool)
    
    id_p = tp_process_muts %>% filter(cluster_id_process == p) %>% pull(mutation_id)
    id_t = tp_tool_muts %>% filter(cluster_id_tool == t) %>% pull(mutation_id)
    jaccard_df <- bind_rows(jaccard_df, tibble(cluster = d, jaccard = jaccard(id_p, id_t), sample = s))
  }
}

results <- lapply(names(df_comparison), FUN = function(s){
  df_comparison[[s]] %>% bind_rows()
}) %>% bind_rows()

results <- lapply(names(df_comparison), FUN = function(s){
  df_comparison[[s]] %>% bind_rows()
  }) %>% bind_rows() %>% 
  mutate(cov = cov, pur = pur, tool = tool, sig_tool = signature_tool, cna_caller = cna_caller, mut_caller = mut_caller, spn = spn) %>% 
  left_join(jaccard_df)

saveRDS(object = results,file = paste0('assign_res_univariate/',spn, '/',spn,'_',cov, 'x_', pur, 'p_',tool,'_',signature_tool,'.rds'))

