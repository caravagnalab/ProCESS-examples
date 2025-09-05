library(ProCESS)
library(ggplot2)
library(dplyr)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/nf-validation/bin/somatic/utils/plot_utils.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/colors.R")

# source("compute_FGA.R")
scout_dir <-"/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
SPNS <- c("SPN01","SPN02","SPN03","SPN04","SPN06","SPN07")
COVERAGES <- c("50","100")
PURITIES <- c("0.3","0.6","0.9")


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
  
  spn = SPN
  validation_dir_somatic <- file.path(validation_dir,spn,"signature")
  df <- lapply(1:nrow(params_grid), function(i) {
    coverage <- params_grid[i, ]$coverage
    purity <- params_grid[i, ]$purity
    ctx <- params_grid[i, ]$context
    combination <- paste0(coverage, "x_", purity)
    process_muts <- readRDS(get_mutations(spn = spn,coverage = coverage,purity = purity,type = "tumour")) %>% 
      dplyr::filter(classes!="germinal") %>% 
      select(contains("occurrences")) %>%         # keep only occurrence columns
      summarise(across(everything(), ~ sum(.x > 0)))
      

    results_folder_path <- file.path(validation_dir_somatic, combination,paste0(vcf_caller,"_",cna_caller))
    file_name_cosine <- file.path(results_folder_path, paste0("cosine_mse_",ctx,".rds"))
    file_name_metrics <- file.path(results_folder_path, paste0("metrics_",ctx,"_sample.rds"))
    
    if (file.exists(file_name_cosine)) {
      cosine <- readRDS(file_name_cosine) %>% 
        dplyr::mutate(context=ctx) %>% 
        dplyr::mutate(total_mutations=n_muts)
      metrics <- readRDS(file_name_metrics) %>% unique() %>% 
        dplyr::mutate(context=ctx)
      all_combinations_SPN <- inner_join(cosine,metrics,by=c("sample","caller","spn","coverage","purity","context"))
      
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



samples <- sapply(strsplit(colnames(list_precision), ":"), `[`, 1)
tools <- sapply(strsplit(colnames(list_precision), ":"), `[`, 2)
col_tools <- method_colors

spn_ids <- sapply(strsplit(colnames(list_precision), "_"), `[`, 1)
#unique_spns <- unique(spn_ids)


column_ha <- HeatmapAnnotation(
  spn = spn_ids, 
  col = list(spn=SPN_colors)
)


row_ha <- rowAnnotation(
  coverage = coverages,
  purity = purities,
  # mutation_type =mut_types,
  #coverage = anno_simple(coverages, col= col_coverages),
  col = list(coverage=col_coverages, purity=col_purities)
)

right_ha <- rowAnnotation(
  metric= rep("sensitivity",nrow(list_sensitivity))
)
####
col_fun = circlize::colorRamp2(c(0, 1), c("white", "#6DA16A"))
h_prec = ComplexHeatmap::Heatmap(list_precision,cluster_rows = F,cluster_columns = F,
                                top_annotation = column_ha,col=col_fun, left_annotation = row_ha,show_column_names = F,
                                show_row_names = F,rect_gp = gpar(col = "black", lwd = 1),
                                #right_annotation = right_ha,
                                column_split = tools,
                                row_split = mut_types,
                                column_gap = unit(4, "mm"),
                                row_gap = unit(2, "mm"),
                                name = "Caller precision"
)

col_fun = circlize::colorRamp2(c(0, 1), c("white", "darkorange"))
h_sen = ComplexHeatmap::Heatmap(list_sensitivity,cluster_rows = F,cluster_columns = F,
                                 col=col_fun, left_annotation = row_ha,show_column_names = F,
                                 show_row_names = F,rect_gp = gpar(col = "black", lwd = 1),
                                 #right_annotation = right_ha,
                                 column_split = tools,
                                 row_split = mut_types,
                                 column_gap = unit(4, "mm"),
                                 row_gap = unit(2, "mm"),
                                 name = "Caller sensitivity"
)

col_fun = circlize::colorRamp2(c(0, 1), c("white", "#B8799B"))
h_cosine = ComplexHeatmap::Heatmap(list_cosine,cluster_rows = F,cluster_columns = F,
                                col=col_fun, left_annotation = row_ha,show_column_names = F,
                                show_row_names = F,rect_gp = gpar(col = "black", lwd = 1),
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
