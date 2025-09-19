rm(list=ls())
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
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/nf-validation/bin/signatures/utils_getters.R")

# source("compute_FGA.R")
scout_dir <-"/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
SPNS <- c("SPN01","SPN02","SPN03","SPN04","SPN06")
# SPNS <- c("SPN04","SPN03")

COVERAGES <- c("50","100","150")
PURITIES <- c("0.3","0.6","0.9")
After_therapy_samples <- c("SPN04_2.1","SPN06_2.1","SPN06_3.1","SPN06_3.2")


cna_caller = "ascat"
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
      mutate(mutation_count=log10(mutation_count)) %>% 
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


df_all_SPN_exposure_sbs <- list()
df_all_SPN_exposure_id <- list()
for (SPN in SPNS){
  message(paste0("Processing ",SPN))
  spn = SPN
  validation_dir_somatic <- file.path(validation_dir,spn,"signature")
  ctxs_base <- gsub("[0-9]", "", CONTEXTS)
  coverage <- 100
  purity <- 0.9
  process_exposures <- get_process_exposures(spn = spn,coverage = 100,purity = 0.9)
  process_exposures_sbs_m <- process_exposures$SBS %>% 
    mutate(across(everything(), ~replace_na(., 0))) %>% 
    select(!c("Sample_ID")) %>% 
    as.matrix() %>% t()
  m_sbs = cbind(process_exposures_sbs_m,process_exposures_sbs_m)
  process_exposures_id_m <- process_exposures$ID %>% 
    mutate(across(everything(), ~replace_na(., 0))) %>% 
    select(!c("Sample_ID")) %>% 
    as.matrix() %>% t()
  m_id = cbind(process_exposures_id_m,process_exposures_id_m)
  df_all_SPN_exposure_sbs[[spn]] <- process_exposures_sbs_m
  df_all_SPN_exposure_id[[spn]] <- process_exposures_id_m
}


all_signatures <- unique(unlist(lapply(df_all_SPN_exposure_sbs, rownames)))

df_all_SPN_exposure_sbs <- lapply(df_all_SPN_exposure_sbs, function(mat) {
  # Create a full matrix with all signatures
  full_mat <- matrix(0, 
                     nrow = length(all_signatures), 
                     ncol = ncol(mat),
                     dimnames = list(all_signatures, colnames(mat)))
  # Fill the existing rows back in
  full_mat[rownames(mat), ] <- mat
  return(full_mat)
})

all_signatures <- unique(unlist(lapply(df_all_SPN_exposure_id, rownames)))

df_all_SPN_exposure_id <- lapply(df_all_SPN_exposure_id, function(mat) {
  # Create a full matrix with all signatures
  full_mat <- matrix(0, 
                     nrow = length(all_signatures), 
                     ncol = ncol(mat),
                     dimnames = list(all_signatures, colnames(mat)))
  # Fill the existing rows back in
  full_mat[rownames(mat), ] <- mat
  return(full_mat)
})



m_sbs_all <- do.call("cbind",df_all_SPN_exposure_sbs)
m_id_all <- do.call("cbind",df_all_SPN_exposure_id)

df_all_combs_SPN<- df_all_combs_SPN %>% 
  arrange(caller) %>% 
  mutate(id=paste(sample,caller,sep=":")) %>% 
  mutate(comb=paste(coverage,purity,context,sep=":")) %>% 
  mutate(Therapy= case_when(sample%in%After_therapy_samples ~ "Therapy",
                        TRUE ~ "not Therapy")) 

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
samples <- sapply(strsplit(colnames(list_precision), ":"), `[`, 1)

coverages <- as.numeric(sapply(strsplit(rownames(list_precision), ":"), `[`, 1))
purities <- as.numeric(sapply(strsplit(rownames(list_precision), ":"), `[`, 2))
mut_types <- sapply(strsplit(rownames(list_precision), ":"), `[`, 3)


therapy_map <- df_all_combs_SPN %>% distinct(sample, Therapy)
therapy_values <- therapy_map$Therapy[match(samples, therapy_map$sample)]

tools <- sapply(strsplit(colnames(list_precision), ":"), `[`, 2)
col_tools <- method_colors
indel_count_map <- df_all_combs_SPN %>% distinct(sample, INDEL_count)
indel_count_values <- indel_count_map$INDEL_count[match(samples, indel_count_map$sample)]
snv_count_map <- df_all_combs_SPN %>% distinct(sample, SNV_count)
snv_count_values <- snv_count_map$SNV_count[match(samples, snv_count_map$sample)]
muts_count_matrix<- rbind(indel_count_values, snv_count_values) %>% t()

spn_ids <- sapply(strsplit(colnames(list_precision), "_"), `[`, 1)



graphics = list(
  "Therapy" = function(x, y, w, h) {
    grid.points(x, y, gp = gpar(col = "black"), pch = 8)
  },
  "not Therapy" = function(x, y, w, h) {
    grid.points(x, y, gp = gpar(col = "white"), pch = 8)
  }
)
column_ha <- HeatmapAnnotation(
  sbs = anno_barplot(t(cbind(m_sbs_all,m_sbs_all)), gp = gpar(fill = sbs_colors), 
                     bar_width = 1, height = unit(2, "cm")),
  id = anno_barplot(t(cbind(m_id_all,m_id_all)), gp = gpar(fill = id_colors), 
                     bar_width = 1, height = unit(2, "cm")),
  # id = anno_barplot(t(m_id_all), gp = gpar(fill = 6:10), 
  #                    bar_width = 1, height = unit(2, "cm")),
  mut_counts = anno_points(muts_count_matrix,gp = gpar(col = c("grey","black")), add_points = TRUE, pt_gp = gpar(col = 5:6), pch = c(16, 16)),
  spn = spn_ids,
  Therapy = anno_customize(
    therapy_values, 
    graphics = graphics   # symbols drawn above the bar
  ),
  col = list(spn=SPN_colors),
  gap = unit(2, "mm")
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

right_ha <- rowAnnotation(
  metric= rep("sensitivity",nrow(list_sensitivity))
)
####
col_fun = circlize::colorRamp2(c(min(list_precision), max(list_precision)), c("#6DA16A","#FCF8F5"))
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

col_fun = circlize::colorRamp2(c(min(list_sensitivity), max(list_sensitivity)), c("darkorange","#FCF8F5"))
h_sen = ComplexHeatmap::Heatmap(list_sensitivity,cluster_rows = F,cluster_columns = F,
                                 col=col_fun, left_annotation = row_ha,show_column_names = F,
                                 show_row_names = F,rect_gp = gpar(col = "grey", lwd = 0.5),
                                 #right_annotation = right_ha
                                 column_split = tools,
                                 row_split = mut_types,
                                 column_gap = unit(4, "mm"),
                                 row_gap = unit(2, "mm"),
                                 name = "Caller sensitivity"
)

col_fun = circlize::colorRamp2(c(min(list_cosine), max(list_cosine)), c("#B8799B","#FCF8F5"))
h_cosine = ComplexHeatmap::Heatmap(list_cosine,cluster_rows = F,cluster_columns = F,
                                col=col_fun, left_annotation = row_ha,show_column_names = F,
                                show_row_names = F,rect_gp = gpar(col = "grey", lwd = 0.5),
                                #right_annotation = right_ha,
                                column_split = tools,
                                # bottom_annotation = column_bottom_ha,
                                row_split = mut_types,
                                column_gap = unit(4, "mm"),
                                row_gap = unit(2, "mm"),
                                name = "Caller cosine similarity"
)
h_final_signatures <- h_prec %v% h_sen %v% h_cosine
pdf("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/Final_Signature_SCOUT_Validation.pdf",width = 10,height = 10)
draw(object = h_final_signatures)
dev.off()
