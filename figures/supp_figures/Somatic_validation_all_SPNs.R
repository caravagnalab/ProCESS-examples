rm(list=ls())
library(ProCESS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/nf-validation/bin/somatic/utils/plot_utils.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")
# source("compute_FGA.R")
scout_dir <-"/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
SPNS <- c("SPN01","SPN02","SPN03","SPN04","SPN05","SPN06","SPN07")
# SPNS <- cC("SPN01","SPN02","SPN03","SPN04")
COVERAGES <- c("50","100","150")
PURITIES <- c("0.3","0.6","0.9")
params_grid = expand.grid(COVERAGES, PURITIES)
colnames(params_grid) = c("coverage", "purity")




####
MUT_TYPES = c("INDEL", "SNV")

params_grid = expand.grid(COVERAGES, PURITIES, MUT_TYPES)
colnames(params_grid) = c("coverage", "purity", "mut")
validation_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION"

i = 1
df_all_SPN <- list()
for (SPN in SPNS){
  
  spn = SPN
  validation_dir_somatic <- file.path(validation_dir,spn,"somatic")
  df <- lapply(1:nrow(params_grid), function(i) {
    coverage <- params_grid[i, ]$coverage
    purity <- params_grid[i, ]$purity
    mut_type <- params_grid[i, ]$mut
    combination <- paste0(coverage, "x_", purity,"p")
    
    results_folder_path <- file.path(validation_dir_somatic, combination, "report", mut_type)
    file_name <- file.path(results_folder_path, "metrics.rds")
    
    mutations_count_file <- file.path("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/mutations_counts/",spn,paste0("mutations_counts_",combination,".rds"))
    mutations_count <- readRDS(mutations_count_file) %>% 
      pivot_wider(
        id_cols = c(sample, coverage, purity),
        names_from = type,
        values_from = mutation_count,
        values_fill = 0,   # replace missing values with 0
        names_glue = "{type}_count"
      ) %>% 
      select(!c("coverage","purity"))
    if (file.exists(file_name)) {
      metrics <- readRDS(file_name)

      parsed_metrics <- lapply(names(metrics), function(caller_name) {
        caller_metrics <- metrics[[caller_name]]
        
        lapply(names(caller_metrics), function(sample_id) {
          #overall <- caller_metrics[[sample_id]]$overall_metrics
          # fdr = FP / (TP + FP)
          detection_summary <- caller_metrics[[sample_id]]$detection_summary
          fdr = detection_summary["False Positive"]/(detection_summary["False Positive"]+detection_summary["True Positive"])
          tpr = detection_summary["True Positive"]/(detection_summary["False Negative"]+detection_summary["True Positive"])
          binned <- caller_metrics[[sample_id]]$performance_table
          binned <- binned %>%
            dplyr::mutate(CCF_bin_class = case_when(
              CCF_bin == "Clonal" ~ "Clonal",
              CCF_bin %in% c("10-25%", "25-50%", "50-99%") ~ "Subclonal High CCF",
              TRUE ~ "Subclonal Low CCF"
            )) %>%
            dplyr::group_by(CCF_bin_class) %>%
            dplyr::summarise(
              mean_sensitivity = mean(sensitivity),
              .groups = "drop"
            ) %>%
            dplyr::mutate(sample = sample_id,
                          caller = caller_name,
                          coverage = coverage,
                          purity = purity,
                          mut_type = mut_type,
                          spn = spn) %>% 
            dplyr::mutate(TPR = tpr,
                          FDR = fdr) %>% 
            inner_join(y = mutations_count,by = "sample")
            

        }) %>% bind_rows()
      }) %>% bind_rows()
      
      return(parsed_metrics)
    } else {
      message("File not found: ", file_name)
      return(NULL)
    }
  }) %>% bind_rows()
  df_all_SPN[[spn]] <- df
}
df_all_combs_SPN <- do.call("rbind",df_all_SPN)
df_all_combs_SPN_somatic <- df_all_combs_SPN


###### get annotation for heatmap ######


create_heatmap <- function(df_all_combs_SPN,name_measure,max_color,column_annotation,row_title,mutation_type,CCF_class){
  
  df_all_combs_SPN<- df_all_combs_SPN %>% 
    dplyr::filter(mut_type==mutation_type) %>% 
    dplyr::mutate(CCF_bin_class = factor(CCF_bin_class, 
                                         levels = c("Subclonal Low CCF", 
                                                    "Subclonal High CCF", 
                                                    "Clonal"))) %>% 
    arrange(CCF_bin_class,sample) %>% 
    mutate(id=paste(sample,caller,sep=":")) %>% 
    mutate(comb=paste(coverage,purity,sep=":")) %>% 
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  


  all_ids <- df_all_combs_SPN %>% pull(id) %>% unique()
  
  # Pivot wider with missing values filled as NA
  list_measure <- df_all_combs_SPN %>%
    filter(CCF_bin_class == CCF_class) %>%
    select(comb, id, !!sym(name_measure)) %>%
    pivot_wider(
      names_from = id,
      values_from = !!sym(name_measure),
      values_fill = setNames(list(NA), name_measure)  # <- correct tidyverse way
    ) %>%
    tibble::column_to_rownames("comb") %>%
    as.matrix()
  
  # Add missing columns if some IDs are absent in this class
  missing_ids <- setdiff(all_ids, colnames(list_measure))
  if(length(missing_ids) > 0){
    list_measure <- cbind(
      list_measure,
      matrix(NA, nrow = nrow(list_measure), ncol = length(missing_ids),
             dimnames = list(rownames(list_measure), missing_ids))
    )
  }
  
  # Reorder columns to consistent order
  list_measure <- list_measure[, all_ids, drop = FALSE]
  samples <- sapply(strsplit(colnames(list_measure), ":"), `[`, 1)
  tools <- sapply(strsplit(colnames(list_measure), ":"), `[`, 2)
  col_tools <- method_colors
  if (mutation_type=="SNV"){
    muts_count_map <- df_all_combs_SPN %>% distinct(sample, SNV_count)
    muts_count_values <- muts_count_map$SNV_count[match(samples, muts_count_map$sample)]
  } else {
    muts_count_map <- df_all_combs_SPN %>% distinct(sample, INDEL_count)
    muts_count_values <- muts_count_map$INDEL_count[match(samples, muts_count_map$sample)]
  }
  spn_ids <- sapply(strsplit(colnames(list_measure), "_"), `[`, 1)
  coverages <- as.numeric(sapply(strsplit(rownames(list_measure), ":"), `[`, 1))
  purities <- as.numeric(sapply(strsplit(rownames(list_measure), ":"), `[`, 2))
  #mut_types <- sapply(strsplit(rownames(list_mean_sensitivity), ":"), `[`, 3)
  # ccf_classes <- sapply(strsplit(rownames(list_mean_sensitivity), ":"), `[`, 3)
  
  column_ha <- HeatmapAnnotation(
    # ccf_classes = ccf_classes,
    muts_count = anno_barplot(muts_count_values,beside = T, border = F, bar_width = 1, gp = gpar(fill = '#DBD7D2', col = 'white'), height = unit(1, "cm")),
    spn = spn_ids, 
    col = list(spn=SPN_colors,ccf_classes=col_ccf_classes),
    annotation_label = c(muts_count=paste0(mutation_type," count"),spn="SPN ID")
  )
  
  
  row_ha <- rowAnnotation(
    coverage = coverages,
    purity = purities,
    # mutation_type =mut_types,
    #coverage = anno_simple(coverages, col= col_coverages),
    col = list(coverage=coverage_colors, purity=purity_colors),
    show_annotation_name = F
  )
  
  right_ha <- rowAnnotation(
    metric= rep("sensitivity",nrow(list_measure))
  )
  col_fun = circlize::colorRamp2(c(min(list_measure,na.rm=TRUE), max(list_measure,na.rm=TRUE)), c("white", max_color))
  if (column_annotation==TRUE){
    ht = ComplexHeatmap::Heatmap(list_measure,cluster_rows = F,cluster_columns = F,
                                 col=col_fun,
                                 top_annotation = column_ha,
                                 left_annotation = row_ha,show_column_names = F,
                                 show_row_names = F,rect_gp = gpar(col = "white", lwd = 1),
                                 column_split = tools,
                                 column_gap = unit(4, "mm"),
                                 row_title = row_title,
                                 row_title_side = "right",
                                 row_gap = unit(2, "mm"),
                                 name = row_title,row_title_rot = 0)
  } else {
    ht = ComplexHeatmap::Heatmap(list_measure,cluster_rows = F,cluster_columns = F,
                                 col=col_fun,
                                 left_annotation = row_ha,show_column_names = F,
                                 show_row_names = F,rect_gp = gpar(col = "white", lwd = 1),
                                 column_split = tools,
                                 column_gap = unit(4, "mm"),
                                 row_title = row_title,
                                 row_title_side = "right",
                                 row_gap = unit(2, "mm"),
                                 name = row_title,row_title_rot = 0)
  }
  return(ht)
  
  
}

h_mean_sensitivity_low_subclonal <- create_heatmap(df_all_combs_SPN=df_all_combs_SPN,mutation_type = "SNV",
                                                   name_measure = "mean_sensitivity" ,max_color = "goldenrod1",
                                                   column_annotation = T,row_title = "Subclonal CCF <10%",
                                                   CCF_class = "Subclonal Low CCF")
h_mean_sensitivity_high_subclonal <- create_heatmap(df_all_combs_SPN=df_all_combs_SPN,mutation_type = "SNV",
                                                   name_measure = "mean_sensitivity",max_color = "deepskyblue4",
                                                   column_annotation = F,row_title = "Subclonal CCF >10%",
                                                   CCF_class = "Subclonal High CCF")
h_mean_sensitivity_clonal <- create_heatmap(df_all_combs_SPN=df_all_combs_SPN,mutation_type = "SNV",
                                                    name_measure = "mean_sensitivity",max_color = "seagreen4",
                                                    column_annotation = F,row_title = "Clonal",CCF_class = "Subclonal High CCF")

h_fdr <- create_heatmap(df_all_combs_SPN=df_all_combs_SPN,mutation_type = "SNV",
                        name_measure = "FDR",max_color = "coral3",
                        column_annotation = F,row_title = "FDR",CCF_class = "Clonal")

h_final_somatic_SNV <- h_mean_sensitivity_low_subclonal %v% h_mean_sensitivity_high_subclonal %v% h_mean_sensitivity_clonal %v% h_fdr
########## INDEL 
h_mean_sensitivity_low_subclonal <- create_heatmap(df_all_combs_SPN=df_all_combs_SPN,mutation_type = "INDEL",
                                                   name_measure = "mean_sensitivity" ,max_color = "goldenrod1",
                                                   column_annotation = T,row_title = "Subclonal CCF <10%",
                                                   CCF_class = "Subclonal Low CCF")
h_mean_sensitivity_high_subclonal <- create_heatmap(df_all_combs_SPN=df_all_combs_SPN,mutation_type = "INDEL",
                                                    name_measure = "mean_sensitivity",max_color = "deepskyblue4",
                                                    column_annotation = F,row_title = "Subclonal CCF >10%",
                                                    CCF_class = "Subclonal High CCF")
h_mean_sensitivity_clonal <- create_heatmap(df_all_combs_SPN=df_all_combs_SPN,mutation_type = "INDEL",
                                            name_measure = "mean_sensitivity",max_color = "seagreen4",
                                            column_annotation = F,row_title = "Clonal",CCF_class = "Subclonal High CCF")

h_fdr <- create_heatmap(df_all_combs_SPN=df_all_combs_SPN,mutation_type = "INDEL",
                        name_measure = "FDR",max_color = "coral3",
                        column_annotation = F,row_title = "FDR",CCF_class = "Clonal")
h_final_somatic_INDEL <- h_mean_sensitivity_low_subclonal %v% h_mean_sensitivity_high_subclonal %v% h_mean_sensitivity_clonal %v% h_fdr
h_final_somatic_SNV
h_final_somatic <- h_final_somatic_SNV %v% h_final_somatic_INDEL

# pdf("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/supplementary/Final_Somatic_SCOUT_Validation_INDEL.pdf",width = 9,height = 4)
# draw(object = h_final_somatic_INDEL,
#      heatmap_legend_side = "bottom",
#      annotation_legend_side = "bottom",
#      merge_legends = TRUE )
# dev.off()
# pdf("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/supplementary/Final_Somatic_SCOUT_Validation_SNV.pdf",width = 9,height = 4)
# draw(object = h_final_somatic_SNV,
#      heatmap_legend_side = "bottom",
#      annotation_legend_side = "bottom",
#      merge_legends = TRUE )
# dev.off()
pdf("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/supp_figures/Final_Somatic_SCOUT_Validation.pdf",
    width = 12,height = 10)
draw(object = h_final_somatic,
     #heatmap_legend_side = "bottom",
     #annotation_legend_side = "bottom",
     merge_legends = TRUE )
dev.off()


