library(ProCESS)
library(ggplot2)
library(dplyr)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
library(tidyr)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/nf-validation/bin/somatic/utils/plot_utils.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

# source("compute_FGA.R")
scout_dir <-"/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
SPNS <- c("SPN01","SPN02","SPN03","SPN04","SPN06")
# SPNS <- c("SPN04")

COVERAGES <- c("50","100","150")
PURITIES <- c("0.3","0.6","0.9")
# After_therapy_samples <- c("SPN06_2.1","SPN07_2.1","SPN07_2.2")


cna_caller = "sequenza"
vcf_caller ="mutect2"
params_grid = expand.grid(COVERAGES, PURITIES)
colnames(params_grid) = c("coverage", "purity")

CONTEXTS = c("SBS96", "ID83")

params_grid = expand.grid(COVERAGES, PURITIES, CONTEXTS)
colnames(params_grid) = c("coverage", "purity", "context")
validation_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION"


i = 1
df_all_SPN <- list()
for (SPN in SPNS){
  message(paste0("Processing ",SPN))
  spn = SPN
  validation_dir_somatic <- file.path(validation_dir,spn,"signature")
  df <- lapply(1:nrow(params_grid), function(i) {
    coverage <- params_grid[i, ]$coverage
    purity <- params_grid[i, ]$purity
    ctx <- params_grid[i, ]$context
    combination <- paste0(coverage, "x_", purity)
    combination1 <- paste0(coverage, "x_", purity,"p")
    message(paste0("Processing combination ",combination))
    mutations_count_file <- file.path("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/mutations_counts/",spn,paste0("mutations_counts_",combination1,".rds"))
    mutations_count <- readRDS(mutations_count_file) %>% 
      pivot_wider(
        id_cols = c(sample, coverage, purity),
        names_from = type,
        values_from = mutation_count,
        values_fill = 0,   # replace missing values with 0
        names_glue = "{type}_count"
      ) %>% 
      select(!c("coverage","purity"))

    results_folder_path <- file.path(validation_dir_somatic, combination,paste0(vcf_caller,"_",cna_caller))
    file_name_cosine <- file.path(results_folder_path, paste0("cosine_mse_",ctx,".rds"))
    file_name_metrics <- file.path(results_folder_path, paste0("metrics_",ctx,"_sample.rds"))
    
    if (file.exists(file_name_cosine)) {
      cosine <- readRDS(file_name_cosine) %>% 
        dplyr::mutate(context=ctx)
      metrics <- readRDS(file_name_metrics) %>% unique() %>% 
        dplyr::mutate(context=ctx)
      all_combinations_SPN <- inner_join(cosine,metrics,by=c("sample","caller","spn","coverage","purity","context")) %>% 
        inner_join(y = mutations_count,by = "sample")
      
      return(all_combinations_SPN)
    } else {
      message("File not found: ", file_name_cosine)
      return(NULL)
    }
  }) %>% bind_rows()
  df_all_SPN[[spn]] <- df
}
df_all_combs_SPN <- do.call("rbind",df_all_SPN)


df_all_combs_SPN<- df_all_combs_SPN %>% 
  arrange(caller) %>% 
  mutate(id=paste(sample,caller,sep=":")) %>% 
  mutate(comb=paste(coverage,purity,context,sep=":"))
  # mutate(Therapy= case_when(sample%in%After_therapy_samples ~ "Therapy",
  #                       TRUE ~ "not Therapy"))


list_precision <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, precision) %>%
  tidyr::pivot_wider(names_from = id, values_from = precision) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()

list_sensitivity <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, recall) %>%
  tidyr::pivot_wider(names_from = id, values_from = recall) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()



list_cosine <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, cosine) %>%
  tidyr::pivot_wider(names_from = id, values_from = cosine) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()
###### get annotation for heatmap ######
coverages <- as.numeric(sapply(strsplit(rownames(list_precision), ":"), `[`, 1))

purities <- as.numeric(sapply(strsplit(rownames(list_precision), ":"), `[`, 2))
mut_types <- sapply(strsplit(rownames(list_precision), ":"), `[`, 3)

#therapy_map <- df_all_combs_SPN %>% distinct(sample, Therapy)
#therapy_values <- therapy_map$Therapy[match(samples, therapy_map$sample)]

samples <- sapply(strsplit(colnames(list_precision), ":"), `[`, 1)
tools <- sapply(strsplit(colnames(list_precision), ":"), `[`, 2)
col_tools <- method_colors
indel_count_map <- df_all_combs_SPN %>% distinct(sample, INDEL_count)
indel_count_values <- indel_count_map$INDEL_count[match(samples, indel_count_map$sample)]
snv_count_map <- df_all_combs_SPN %>% distinct(sample, SNV_count)
snv_count_values <- snv_count_map$SNV_count[match(samples, snv_count_map$sample)]
muts_count_matrix<- rbind(indel_count_values, snv_count_values) %>% t()

spn_ids <- sapply(strsplit(colnames(list_precision), "_"), `[`, 1)
#unique_spns <- unique(spn_ids)
# 
# graphics = list(
#   "Therapy" = function(x, y, w, h) {
#     grid.points(x, y, gp = gpar(col = "black"), pch = 16)
#   },
#   "not Therapy" = function(x, y, w, h) {
#     grid.points(x, y, gp = gpar(col = "white"), pch = 16)
#   }
# )

column_ha <- HeatmapAnnotation(
  mut_counts = anno_barplot(muts_count_matrix, 
                     beside = TRUE),
  spn = spn_ids, 
  col = list(spn=SPN_colors)
)

# column_bottom_ha <- HeatmapAnnotation(
#   Therapy = anno_customize(therapy_values, graphics = graphics)
# )

row_ha <- rowAnnotation(

  coverage = coverages,
  purity = purities,
  col = list(coverage=coverage_colors, purity=purity_colors),
  show_annotation_name = F
)

right_ha <- rowAnnotation(
  metric= rep("sensitivity",nrow(list_sensitivity))
)
####
col_fun = circlize::colorRamp2(c(min(list_precision), max(list_precision)), c("white", "#6DA16A"))
h_prec = ComplexHeatmap::Heatmap(list_precision,cluster_rows = F,cluster_columns = F,
                                top_annotation = column_ha,col=col_fun, left_annotation = row_ha,show_column_names = F,
                                show_row_names = F,rect_gp = gpar(col = "grey", lwd = 0.5),
                                #right_annotation = right_ha,
                                column_split = tools,
                                row_split = mut_types,
                                column_gap = unit(4, "mm"),
                                row_gap = unit(2, "mm"),
                                name = "Caller precision"
)

col_fun = circlize::colorRamp2(c(min(list_sensitivity), max(list_sensitivity)), c("white", "darkorange"))
h_sen = ComplexHeatmap::Heatmap(list_sensitivity,cluster_rows = F,cluster_columns = F,
                                 col=col_fun, left_annotation = row_ha,show_column_names = F,
                                 show_row_names = F,rect_gp = gpar(col = "grey", lwd = 0.5),
                                 #right_annotation = right_ha,
                                 column_split = tools,
                                 row_split = mut_types,
                                 column_gap = unit(4, "mm"),
                                 row_gap = unit(2, "mm"),
                                 name = "Caller sensitivity"
)

col_fun = circlize::colorRamp2(c(min(list_cosine), max(list_cosine)), c("white", "#B8799B"))
h_cosine = ComplexHeatmap::Heatmap(list_cosine,cluster_rows = F,cluster_columns = F,
                                col=col_fun, left_annotation = row_ha,show_column_names = F,
                                show_row_names = F,rect_gp = gpar(col = "grey", lwd = 0.5),
                                #right_annotation = right_ha,
                                column_split = tools,
                                row_split = mut_types,
                                column_gap = unit(4, "mm"),
                                row_gap = unit(2, "mm"),
                                name = "Caller cosine similarity"
)
h_final_signatures <- h_prec %v% h_sen %v% h_cosine
pdf("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/Final_Signature_SCOUT_Validation.pdf",width = 10,height = 10)
draw(object = h_final_signatures)
dev.off()
