pkgs <- c("ProCESS", "tidyverse", "ggplot2", 
          "caret", "ggtext", "reshape2", 
          "lsa", "Metrics", "ggalluvial", 
          "patchwork")
sapply(pkgs, require, character.only = TRUE)

source("../../getters/process_getters.R")
source("../../getters/tumourevo_getters.R")
source("utils/utils_getters.R")
source("utils/utils_validation.R")
source("utils/utils.R")
source("utils/utils_plots.R")


##### Use getters to get a pairs of RACES and tumourevo data #####
### Get ProCESS exposure data ###
base_path <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
coverage_list <- c(50)
purity_list <- c(0.9)
dataset <- "SCOUT"
spn <- c("SPN01", "SPN02", "SPN03")
context <- 'ID'
if (context == 'ID'){
  context_sig = 'ID83'
}
vcf_caller <- "mutect2"
cna_caller <- "ascat"

# gt
ground_truth_nested <- list()
for (spn_id in spn) {
  print(spn_id)
  
  if (is.null(ground_truth_nested[[spn_id]])) {
    ground_truth_nested[[spn_id]] <- list()
  }
  
  for (cov in coverage_list) {
    print(cov)
    coverage_key <- paste0("coverage_", cov)
    
    if (is.null(ground_truth_nested[[spn_id]][[coverage_key]])) {
      ground_truth_nested[[spn_id]][[coverage_key]] <- list()
    }
    
    for (pur in purity_list) {
      print(pur)
      
      purity_key <- paste0("purity_", pur)
      gt_exposure <-  get_process_exposures(spn = spn_id, 
                                            coverage = cov,
                                            purity = pur)
      ground_truth_nested[[spn_id]][[coverage_key]][[purity_key]] <- gt_exposure[[context]]  %>%
        tibble::column_to_rownames("Sample_ID") %>% 
        as.matrix()
    }
  }
}


### Load tumourevo signature data ###
tumourevo_signature_res <- list()
for (spn_id in spn) {
  tumourevo_signature_res[[spn_id]] <- list()
  
  for (cov in coverage_list) {
    tumourevo_signature_res[[spn_id]][[paste0("coverage_", cov)]] <- list()
    
    for (pur in purity_list) {
      message("Processing ", spn_id, " | Coverage: ", cov, " | Purity: ", pur)
      
      result <- tryCatch({
        # Get SigProfiler paths
        sigprofiler <- get_tumourevo_signatures(
          spn = spn_id,
          coverage = cov,
          purity = pur,
          vcf_caller = vcf_caller,
          cna_caller = cna_caller,
          tool = "SigProfiler",
          context = context_sig,  # ID83!
          base_path = base_path
        )
        
        # Combine and load paths
        paths <- c(
          sigprofiler$COSMIC_exposure,
          sigprofiler$COSMIC_signatures,
          sigprofiler$denovo_exposure,
          sigprofiler$denovo_signatures
        )
        
        data <- load_signature_data(paths)
        
        names(data) <- c(
          "SigProfiler_COSMIC_exposure",
          "SigProfiler_COSMIC_signatures",
          "SigProfiler_denovo_exposure",
          "SigProfiler_denovo_signatures"
        )
        
        data
      }, error = function(e) {
        message("Failed for ", spn, " | cov: ", cov, " | pur: ", pur)
        message("Reason: ", e$message)
        NULL
      })
      
      # Store result if successful
      if (!is.null(result)) {
        tumourevo_signature_res[[spn_id]][[paste0("coverage_", cov)]][[paste0("purity_", pur)]] <- result
      }
    }
  }
}


sigprof_aligned <- align_sigprofiler_res(tumourevo_signature_res)

### Compare exposures of estimated and true signatures  ###
sankey_df <- prepare_sankey_data_id(ground_truth_nested, sigprof_aligned)

plots <- list()
for (spn_id in spn){
  for (cov in coverage_list){
    for (pur in purity_list){
      sankey_plots <- generate_sankey_id(sankey_df = sankey_df, spn_id = spn_id, cov = cov, pur = pur)
      plots[[paste0('coverage_',cov, "_purity_", pur)]][[spn_id]] <- sankey_plots
    }
  }
}

plots_all <- list()
for (cov in coverage_list){
  for (pur in purity_list){
    plots_all[[paste0('coverage_',cov, "_purity_", pur)]] <- wrap_plots(plots[[paste0('coverage_',cov, "_purity_", pur)]], nrow = length(spn)) 
  }
}
plots_all$coverage_50_purity_0.9

## Exposure validation ##
metrics_sample <- tibble()
metrics_spn <- tibble()
cosine_mse <- tibble()
for (spn_id in spn) {
  
  for (cov in coverage_list) {
    c <- paste0("coverage_", cov)
    
    for (pur in purity_list) {
      p <- paste0("purity_", pur)
      inf_sig_sigprofiler <- sigprof_aligned[[spn_id]][[c]][[p]] 
      sim_sig <- ground_truth_nested[[spn_id]][[c]][[p]]  %>% as.data.frame()
      

      sample_metrics_sigprofiler <- per_sample_metrics(inferred_df = inf_sig_sigprofiler, simulated_df = sim_sig) %>% 
        mutate(caller = 'SigProfiler', spn = spn_id, coverage = cov, purity = pur)
      summary_sigprofiler <- summary_stats(sample_metrics_sigprofiler) %>% mutate(caller = 'SigProfiler', spn = spn_id, coverage = cov, purity = pur)
      
      cosine_mse_sigprofiler <- compute_cosine_mse(inferred = inf_sig_sigprofiler, simulated = sim_sig) %>% 
        mutate(caller = 'SigProfiler', spn = spn_id, coverage = cov, purity = pur)

      metrics_sample <- bind_rows(metrics_sample, sample_metrics_sigprofiler)
      metrics_spn <- bind_rows(metrics_spn, summary_sigprofiler)
      
      cosine_mse <- bind_rows(cosine_mse, cosine_mse_sigprofiler)
    }
  }
}

col_spn <- c('SPN01' = 'steelblue', 'SPN02' ='seagreen', 'SPN03' ='goldenrod', 
             'SPN04' ='coral', 'SPN06' ='palevioletred', 'SPN07' ='indianred3')

### Generate final plots ###
cosine_mse_plot <- cosine_mse %>% 
  #pivot_longer(cols = c(mse, cosine)) %>% 
  ggplot() +
  geom_boxplot(aes(x = spn, y = cosine, col = spn)) +
  scale_color_manual(values = col_spn) +
  facet_grid(coverage~purity) + 
  ggtitle(paste0('Simulated vs inferred exposures per sample')) + 
  theme_bw() +
  cosine_mse %>% 
  #pivot_longer(cols = c(mse, cosine)) %>% 
  ggplot() +
  geom_boxplot(aes(x = spn, y = mse, col = spn)) +
  scale_color_manual(values = col_spn) +
  facet_grid(coverage~purity) + 
  theme_bw()  +
  plot_layout(ncol = 2, guides = 'collect')
  
metric_plot <- metrics_sample %>% 
  pivot_longer(cols = c(precision, recall)) %>% 
  ggplot() +
  geom_col(aes(x = name, y = value, fill = spn),position=position_dodge()) +
  scale_fill_manual(values = col_spn) +
  facet_grid(coverage~purity) + 
  ggtitle(paste0('Simulated vs inferred SBS signature')) +
  theme_bw() 

wrap_plots(metric_plot, cosine_mse_plot, nrow = 2) 
