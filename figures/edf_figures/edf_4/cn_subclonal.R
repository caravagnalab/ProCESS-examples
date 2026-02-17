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
SPNS <- c("SPN01","SPN02","SPN03","SPN04", "SPN05", "SPN06", 'SPN07')

COVERAGES <- c("50","100", "150")
PURITIES <- c("0.3","0.6","0.9")

source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
WGD_samples <- c("SPN01_1.1","SPN01_1.3","SPN06_3.1","SPN06_3.2")
fga_df <- readRDS("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/fga_df.rds")



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
)) %>% 
  left_join(fga_df %>% mutate(sample = spn), by = join_by(sample))

stat.test <- df_all_combs_SPN %>%
  filter(tool == 'Battenberg', fgs.y > 2) %>% 
  mutate(fgs_class = ifelse(fgs.y > 30, 'High FGS', 'Low FGS')) %>% 
  rstatix::wilcox_test(correctness_subclonal ~ fgs_class) %>%
  rstatix::adjust_pvalue(method = "BH") %>%       
  rstatix::add_significance("p.adj") %>%          
  rstatix::add_xy_position(x = "fgs_class") %>% 
  mutate(y.position = 100)

subclonal_cna <- df_all_combs_SPN %>%
  filter(tool == 'Battenberg', fgs.y > 2) %>% 
  mutate(fgs_class = ifelse(fgs.y > 30, 'High FGS', 'Low FGS')) %>% 
  ggplot(aes(x = fgs_class, y = correctness_subclonal*100, col = fgs_class)) +
  geom_boxplot(aes(col = fgs_class),size = .5, 
               width = 0.2, outlier.size = .6) +
  geom_jitter(size = 1, width = .1, alpha = .5) +
  #scale_color_manual('CN caller',values = c('ASCAT' = 'mediumpurple3', 'Sequenza' = 'darkorange2', 'Battenberg' = 'steelblue3', 'ProCESS' = 'gray')) +
  scale_color_manual('FGS Class', values = c("High FGS" = "indianred2", "Low FGS"  = "dodgerblue3")) + 
  stat_pvalue_manual(
    stat.test,
    label = "p.adj.signif",   # "***", "**", etc.
    tip.length = 0.01
  ) +
  xlab('') + 
  #ylim(75, 100) +
  ylab('% correct subclonal CN') + 
  my_ggplot_theme()
subclonal_cna
