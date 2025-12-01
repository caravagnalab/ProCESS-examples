library(ProCESS)
library(ggplot2)
library(dplyr)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
# source("compute_FGA.R")
scout_dir <-"/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
SPNS <- c("SPN01","SPN02","SPN03","SPN04", "SPN06", 'SPN07')
COVERAGES <- c("50","100", "150")
PURITIES <- c("0.3","0.6","0.9")
SPN_colors <-c("SPN01"='steelblue', "SPN02"='seagreen', "SPN03"='goldenrod', 
               "SPN04"='coral', "SPN06"='palevioletred', "SPN07"='indianred3')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")
WGD_samples <- c("SPN01_1.1","SPN01_1.3","SPN06_3.1","SPN06_3.2")

validation_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION"
params_grid = expand.grid(COVERAGES, PURITIES)
colnames(params_grid) = c("coverage", "purity")
df_all_SPN <- list()
for (SPN in SPNS){
  spn <- SPN
  df = lapply(1:nrow(params_grid), function(i) {
    coverage = params_grid[i,]$coverage
    purity = params_grid[i,]$purity
    comb <- paste0(coverage,"x_",purity)
    samples <- get_sample_names(spn = spn)
    validation_dir_cna <- file.path(validation_dir,spn,"cna")
    all_metrics_comb <- list()
    for (sample in samples){
      metrics_filename <- paste0(validation_dir_cna,"/",comb,"/",sample,"/metrics.rds")
      metrics_bp_filename <- paste0(validation_dir_cna,"/",comb,"/",sample,"/metrics_bp.rds")
      if (!file.exists(metrics_filename)){
        message("File not found: ", metrics_filename)
        # all_metrics_comb[[sample]] <- NA
      } else{
        metrics_df <- readRDS(metrics_filename) %>% 
          mutate(delta_purity=as.numeric(true_purity)-as.numeric(purity)) %>% 
          mutate(delta_ploidy=as.numeric(true_ploidy)-as.numeric(ploidy)) %>% 
          filter(tool!="cnvkit")
        metrics_bp_df <- readRDS(metrics_bp_filename) %>% 
          filter(chr=="genome") %>% 
          filter(tool!="cnvkit")
        all_metrics_comb[[sample]] <- inner_join(metrics_df,metrics_bp_df,by=c("tool","sample","spn","coverage","fga","fgs","true_purity"))
      }
    }
    all_combinations_SPN <- do.call("rbind",all_metrics_comb)
  })
  df_all_SPN[[spn]] <- do.call("rbind",df)  
}
df_all_combs_SPN <- do.call("rbind",df_all_SPN)
df_all_combs_SPN <- df_all_combs_SPN %>% mutate(tool = case_when(
  tool == 'ascat' ~ 'ASCAT',
  tool == 'sequenza' ~ 'Sequenza',
  tool == 'battenberg' ~ 'Battenberg',
))

cna_corr <- df_all_combs_SPN %>%
  ggplot(aes(x = tool, y = correctness_clonal*100, col = tool)) +
  geom_violin(size = .3, 
              show.legend = F) + 
  geom_boxplot(size = .3, 
               width = 0.1, outlier.size = .6) +
  scale_color_manual('CN caller',values = c('deepskyblue4', 'maroon', 'sienna')) +
  theme_minimal() +
  xlab('Tool') + 
  ylim(70, 100) +
  ylab('% correct CN') + 
  my_ggplot_theme()
