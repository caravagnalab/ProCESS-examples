rm(list=ls())
library(ProCESS)
library(ggplot2)
library(dplyr)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(tidyr)
library(ComplexHeatmap)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/nf-validation/bin/somatic/utils/plot_utils.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")

# source("compute_FGA.R")
scout_dir <-"/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
SPNS <- c("SPN01","SPN02","SPN03","SPN04","SPN06","SPN07")
# SPNS <- cC("SPN01","SPN02","SPN03","SPN04")
COVERAGES <- c("50","100","150")
PURITIES <- c("0.3","0.6","0.9")
params_grid = expand.grid(COVERAGES, PURITIES)
colnames(params_grid) = c("coverage", "purity")
After_therapy_samples <- c("SPN04_2.1","SPN06_2.1","SPN06_3.1","SPN06_3.2","SPN07_2.1","SPN07_2.2")




####
cna_caller = "ascat"
vcf_caller ="mutect2"
params_grid = expand.grid(COVERAGES, PURITIES)
colnames(params_grid) = c("coverage", "purity")
validation_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION"

i = 1
df_all_SPN <- list()
for (SPN in SPNS){
  
  spn = SPN
  validation_dir_driver <- file.path(validation_dir,spn,"driver")
  df <- lapply(1:nrow(params_grid), function(i) {
    coverage <- params_grid[i, ]$coverage
    purity <- params_grid[i, ]$purity
    combination <- paste0(coverage, "x_", purity)
    combination1 <- paste0(coverage, "x_", purity,"p")
    callers_comb <-paste0(vcf_caller,"_",cna_caller)
    mutations_count_file <- file.path("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/mutations_counts/",spn,paste0("mutations_counts_",combination1,".rds"))
    mutations_count <- readRDS(mutations_count_file) %>% 
      mutate(mutation_count=log10(mutation_count)) %>% 
      pivot_wider(
        id_cols = c(sample, coverage, purity),
        names_from = type,
        values_from = mutation_count,
        values_fill = 0,   # replace missing values with 0
        names_glue = "{type}_count"
      ) %>% 
      select(!c("coverage","purity"))
    results_folder_path <- file.path(validation_dir_driver, combination,callers_comb)
    #confusion_matrix_file <- file.path(results_folder_path, "confusion_matrix.rds")
    all_drivers_comparison_file <- file.path(results_folder_path, "all_driver_comparsion.rds")
    if (file.exists(all_drivers_comparison_file)) {
      #confusion_matrix <- readRDS(confusion_matrix_file)
      metrics_df <- readRDS(all_drivers_comparison_file) %>%
        group_by(sample) %>%
        summarise(
          TP = sum(driver_class == "Process True - Tumourevo True"),
          FP = sum(driver_class == "Process False - Tumourevo True"),
          FN = sum(driver_class == "Process True - Tumourevo False"),
          .groups = "drop"
        ) %>%
        mutate(
          precision = ifelse(TP + FP > 0, TP / (TP + FP), NA),
          recall    = ifelse(TP + FN > 0, TP / (TP + FN), NA)
        ) %>% 
        mutate(coverage=coverage) %>% 
        mutate(purity=purity) %>% 
        mutate(spn=spn)
      all_true_drivers <- readRDS(all_drivers_comparison_file) %>% 
        filter(is_driver_process==TRUE) %>% 
        group_by(sample) %>% 
        summarize(driver_mutations = n())
      all_tumourevo_drivers <- readRDS(all_drivers_comparison_file) %>% 
        filter(is_driver_tumourevo==TRUE) %>% 
        group_by(sample) %>% 
        summarize(driver_mutations_tumourevo = n())
      all_combinations_SPN <- inner_join(metrics_df,mutations_count,by=c("sample")) %>% 
        inner_join(all_true_drivers,by="sample") %>% 
        inner_join(all_tumourevo_drivers,by="sample")
      return(all_combinations_SPN)
    } else {
      message("File not found: ", confusion_matrix_file)
      return(NULL)
    }
    
  }) %>% bind_rows()
  df_all_SPN[[spn]] <- df
}
df_all_combs_SPN <- do.call("rbind",df_all_SPN)
driver_counts <- bind_rows(
  lapply(SPNS, function(spn) {
    t <- readRDS(file.path(
      "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/drivers/",
      spn,
      "process_drivers.rds"
    ))
    
    tibble(
      sample = spn,
      spn = spn,
      driver_mutations = t %>%
        filter(type == "SID") %>%
        nrow()
    )
  })
)

df_all_combs_SPN<- df_all_combs_SPN %>% 
  arrange(spn) %>% 
  mutate(id=sample) %>% 
  mutate(comb=paste(coverage,purity,sep=":")) %>% 
  mutate(Therapy= case_when(sample%in%After_therapy_samples ~ "Therapy",
                            TRUE ~ "not Therapy")) 
all_ids <- df_all_combs_SPN %>% pull(id) %>% unique()

list_precision <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, precision) %>%
  tidyr::pivot_wider(names_from = id, values_from = precision) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()
list_recall <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, recall) %>%
  tidyr::pivot_wider(names_from = id, values_from = recall) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()
samples <- sapply(strsplit(colnames(list_precision), ":"), `[`, 1)

spn_ids <- sapply(strsplit(colnames(list_precision), "_"), `[`, 1)
coverages <- as.numeric(sapply(strsplit(rownames(list_precision), ":"), `[`, 1))
purities <- as.numeric(sapply(strsplit(rownames(list_precision), ":"), `[`, 2))
#te_muts <- as.numeric(sapply(strsplit(rownames(list_precision), ":"), `[`, 3))
therapy_map <- df_all_combs_SPN %>% distinct(sample, Therapy)
therapy_values <- therapy_map$Therapy[match(samples, therapy_map$sample)]

indel_count_map <- df_all_combs_SPN %>% distinct(sample, INDEL_count)
indel_count_values <- indel_count_map$INDEL_count[match(samples, indel_count_map$sample)]
snv_count_map <- df_all_combs_SPN %>% distinct(sample, SNV_count)
snv_count_values <- snv_count_map$SNV_count[match(samples, snv_count_map$sample)]
muts_count_matrix<- rbind(indel_count_values, snv_count_values) %>% t()

driver_count_map <- df_all_combs_SPN %>% distinct(sample, driver_mutations)
driver_count_values <- driver_count_map$driver_mutations[match(samples, driver_count_map$sample)]


graphics = list(
  "Therapy" = function(x, y, w, h) {
    grid.points(x, y, gp = gpar(col = "black"), pch = 8)
  },
  "not Therapy" = function(x, y, w, h) {
    grid.points(x, y, gp = gpar(col = "white"), pch = 8)
  }
)

column_ha <- HeatmapAnnotation(
  mut_counts = anno_points(muts_count_matrix,gp = gpar(col = c("grey","black")), add_points = TRUE, pt_gp = gpar(col = 5:6), pch = c(16, 16)),
  driver_count = anno_barplot(driver_count_values,gp = gpar(fill = "navajowhite1", col = NA),bar_width = 0.9),
  spn = spn_ids,
  Therapy = anno_customize(
    therapy_values, 
    graphics = graphics   # symbols drawn above the bar
  ),
  col = list(spn=SPN_colors),
  gap = unit(2, "mm"),
  annotation_label = c("Log Mutation Count", 
                       "N. driver mutations", 
                       "SPN ID","Therapy")
)


column_bottom_ha <- HeatmapAnnotation(
  Therapy = anno_customize(therapy_values, graphics = graphics)
)

row_ha <- rowAnnotation(
  
  coverage = coverages,
  purity = purities,
  col = list(coverage=coverage_colors, purity=purity_colors),
  show_annotation_name = F
)


####


col_fun = circlize::colorRamp2(c(min(list_precision,na.rm=TRUE), max(list_precision,na.rm=TRUE)), c("darkorange","#FCF8F5"))
h_prec = ComplexHeatmap::Heatmap(list_precision,cluster_rows = F,cluster_columns = F,
                                 top_annotation = column_ha,
                                col=col_fun, left_annotation = row_ha,show_column_names = F,
                                show_row_names = F,rect_gp = gpar(col = "grey", lwd = 0.5),
                                # right_annotation = right_ha,
                                column_gap = unit(4, "mm"),
                                row_gap = unit(2, "mm"),
                                name = "Caller precision"
)
col_fun = circlize::colorRamp2(c(min(list_precision,na.rm=TRUE), max(list_precision,na.rm=TRUE)), c("#6DA16A","#FCF8F5"))
h_recall = ComplexHeatmap::Heatmap(list_recall,cluster_rows = F,cluster_columns = F,
                                 col=col_fun, left_annotation = row_ha,show_column_names = F,
                                 show_row_names = F,rect_gp = gpar(col = "grey", lwd = 0.5),
                                 column_gap = unit(4, "mm"),
                                 row_gap = unit(2, "mm"),
                                 name = "Caller recall"
)
h_final_drivers <- h_prec %v% h_recall
pdf("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/Final_Driver_SCOUT_Validation.pdf",width = 10,height = 12)
draw(object = h_final_drivers,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     merge_legends = TRUE)
dev.off()



############ main text #########
counts_tot_mutations <- df_all_combs_SPN %>% 
  group_by(sample) %>% 
  mutate(mean_SNV_count=mean(SNV_count)) %>% 
  mutate(mean_INDEL_count=mean(INDEL_count)) %>%
  select(sample,mean_INDEL_count,mean_SNV_count) %>% 
  unique()
confusion_all <-df_all_combs_SPN %>% 
  group_by(sample) %>% 
  mutate(median_TP=median(TP)) %>% 
  mutate(median_FP=median(FP)) %>% 
  mutate(median_FN=median(FN)) %>% 
  select(sample,spn,median_FN,median_FP,median_TP) %>% 
  unique() %>% 
  ungroup() %>% 
  pivot_longer(
    cols = starts_with("median_"),
    names_to = "metric",
    values_to = "count"
  ) %>%
  mutate(
    metric = recode(metric,
                    median_TP = "TP",
                    median_FP = "FP",
                    median_FN = "FN")
  ) %>% 
  left_join(counts_tot_mutations,by = "sample")
ggplot(confusion_all, aes(x = sample, y = count, fill = metric)) +
  geom_col(position = "dodge") +
  labs(
    x = "Sample",
    y = "Count",
    fill = "Metric"
  ) +
  my_ggplot_theme()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
  

confusion_all %>% 
  ggplot(aes(x=sample,))

cm <- tibble(
  TP = sum(confusion_all$median_TP),
  FP = sum(confusion_all$median_FP),
  FN = sum(confusion_all$median_FN),
  TN = 0
)



df_all_SPN_driver <- list()
for (SPN in SPNS){
  
  spn = SPN
  validation_dir_driver <- file.path(validation_dir,spn,"driver")
  df <- lapply(1:nrow(params_grid), function(i) {
    coverage <- params_grid[i, ]$coverage
    purity <- params_grid[i, ]$purity
    combination <- paste0(coverage, "x_", purity)
    combination1 <- paste0(coverage, "x_", purity,"p")
    callers_comb <-paste0(vcf_caller,"_",cna_caller)
    results_folder_path <- file.path(validation_dir_driver, combination,callers_comb)
    all_drivers_comparison_file <- file.path(results_folder_path, "all_driver_comparsion.rds")
    if (file.exists(all_drivers_comparison_file)) {
      #confusion_matrix <- readRDS(confusion_matrix_file)
      metrics_df <- readRDS(all_drivers_comparison_file)
    }
  }) %>% bind_rows()
  df_all_SPN_driver[[spn]] <- df
}
df_all_combs_SPN_drivers <- do.call("rbind",df_all_SPN_driver)

FN_drivers <- df_all_combs_SPN_drivers %>% 
  filter(driver_class=="Process True - Tumourevo False") %>% 
  select(sample,driver_label) %>% unique()

FP_drivers <- df_all_combs_SPN_drivers %>% 
  filter(driver_class=="Process False - Tumourevo True") %>% 
  select(sample,driver_label) %>% unique()
