rm(list=ls())
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
SPNS <- c("SPN01","SPN02","SPN03","SPN04","SPN06","SPN07")
COVERAGES <- c("50","100")
PURITIES <- c("0.3","0.6","0.9")
SPN_colors <-c("SPN01"='steelblue', "SPN02"='seagreen', "SPN03"='goldenrod', 
               "SPN04"='coral', "SPN06"='palevioletred', "SPN07"='indianred3')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/colors.R")
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
          mutate(delta_ploidy=as.numeric(true_ploidy)-as.numeric(ploidy))
        metrics_bp_df <- readRDS(metrics_bp_filename) %>% 
          filter(chr=="genome")
        all_metrics_comb[[sample]] <- inner_join(metrics_df,metrics_bp_df,by=c("tool","sample","spn","coverage","fga","fgs","true_purity"))
      }
    }
    all_combinations_SPN <- do.call("rbind",all_metrics_comb)
  })
  df_all_SPN[[spn]] <- do.call("rbind",df)  
}
df_all_combs_SPN <- do.call("rbind",df_all_SPN)

######## heatmap purity ########

df_all_combs_SPN<- df_all_combs_SPN %>% 
  arrange(tool) %>% 
  mutate(id=paste(sample,tool,sep=":")) %>% 
  mutate(comb=paste(coverage,true_purity,sep=":")) %>% 
  
  mutate(delta_purity_class = case_when(delta_purity>=0.3  ~ "highly underestimated",
                                        delta_purity<=-0.3 ~ "highly overestimated",
                                        delta_purity>=0.1 & delta_purity<0.3~ "poorly underestimated",
                                        delta_purity<=-0.1 & delta_purity>-0.3~ "poorly overestimated",
                                        is.na(delta_ploidy) ~ "not estimated",
                                        TRUE ~ "correctly estimated",
                                        )) %>% 
  mutate(delta_ploidy_class = case_when(delta_ploidy>1  ~ "highly underestimated",
                                        delta_ploidy<=-1 ~ "highly overestimated",
                                        delta_ploidy>=0.5 & delta_ploidy<1~ "poorly underestimated",
                                        delta_ploidy<=-0.5 & delta_ploidy>-1~ "poorly overestimated",
                                        is.na(delta_ploidy) ~ "not estimated",
                                        TRUE ~ "correctly estimated")) %>% 
  mutate(correctness_clonal_class= case_when(correctness_clonal>=0.95  ~ "Excellent",
                                             correctness_clonal<0.95 & correctness_clonal>=0.80   ~ "High",
                                             correctness_clonal<0.80 & correctness_clonal>=0.60   ~ "Medium",
                                             is.na(delta_ploidy) ~ "not estimated",
                                             TRUE ~ "Low")) %>% 
  mutate(WGD= case_when(sample%in%WGD_samples ~ "WGD",
                        TRUE ~ "not WGD"))
  

## purity matrix

list_purities <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, delta_purity_class) %>%
  tidyr::pivot_wider(names_from = id, values_from = delta_purity_class) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()


## ploidy matrix


list_ploidy <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, delta_ploidy_class) %>%
  tidyr::pivot_wider(names_from = id, values_from = delta_ploidy_class) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()

## correctness clonal matrix

list_correctness_clonal <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, correctness_clonal_class) %>%
  tidyr::pivot_wider(names_from = id, values_from = correctness_clonal_class) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()


## correctness bp detection

list_bp_precision <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, precision) %>%
  tidyr::pivot_wider(names_from = id, values_from = precision) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()

list_bp_recall <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, recall) %>%
  tidyr::pivot_wider(names_from = id, values_from = recall) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()



###### get annotation for heatmap ######


samples <- sapply(strsplit(colnames(list_purities), ":"), `[`, 1)
tools_with_cnvkit <- sapply(strsplit(colnames(list_purities), ":"), `[`, 2)
tools <- sapply(strsplit(colnames(list_purities), ":"), `[`, 2)
col_tools <- c('ascat' = 'palegreen4','cnvkit' = 'orange', 'sequenza' = 'steelblue', 'battenberg' = 'palevioletred')

spn_ids <- sapply(strsplit(colnames(list_purities), "_"), `[`, 1)

#unique_spns <- unique(spn_ids)
#col_spns <-c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
#names(col_spns) <- SPNS

fga_map <- df_all_combs_SPN %>% distinct(sample, fga)
fga_values <- fga_map$fga[match(samples, fga_map$sample)]
fgs_map <- df_all_combs_SPN %>% distinct(sample, fgs)
fgs_values <- fgs_map$fgs[match(samples, fgs_map$sample)]
wgd_map <- df_all_combs_SPN %>% distinct(sample, WGD)
wgd_values <- wgd_map$WGD[match(samples, wgd_map$sample)]
true_ploidy_map <- df_all_combs_SPN %>% distinct(sample, true_ploidy)
true_ploidy_values <- true_ploidy_map$true_ploidy[match(samples, true_ploidy_map$sample)]
col_fun_true_ploidy = circlize::colorRamp2(c(2, 3, 4,5), c("khaki1", "khaki2", "khaki3","khaki4"))

coverages <- as.numeric(sapply(strsplit(rownames(list_purities), ":"), `[`, 1))
purities <- as.numeric(sapply(strsplit(rownames(list_purities), ":"), `[`, 2))



graphics = list(
  "WGD" = function(x, y, w, h) {
    grid.points(x, y, gp = gpar(col = "black"), pch = 16)
  },
  "not WGD" = function(x, y, w, h) {
    grid.points(x, y, gp = gpar(col = "white"), pch = 16)
  }
)

column_ha <- HeatmapAnnotation(
  fga = anno_barplot(fga_values),
  fgs = anno_barplot(fgs_values),
  #coverage = anno_simple(coverages, col= col_coverages),
  spn = spn_ids, 
  # tool = tools,
  col = list(spn=SPN_colors)
)

column_bottom_ha <- HeatmapAnnotation(
  true_ploidy = true_ploidy_values,
  WGD = anno_customize(wgd_values, graphics = graphics),
  #coverage = anno_simple(coverages, col= col_coverages),
  # tool = tools,
  col = list(true_ploidy=col_fun_true_ploidy)
)


row_ha <- rowAnnotation(
  coverage = coverages,
  purity = purities,
  #coverage = anno_simple(coverages, col= col_coverages),
  col = list(coverage=col_coverages, purity=col_purities)
)
####
#col_fun_purity = circlize::colorRamp2(c(-1,0, 1), c("forestgreen","white", "darkorange"))
# col_fun_classes <- c(
#   "correctly estimated"   = "#A5D6A7",  # soft green
#   "highly overestimated"  = "#EF9A9A",  # soft red
#   "poorly overestimated"  = "#FFEBEE",  # very light red/pink
#   "highly underestimated" = "#90CAF9",  # soft blue
#   "poorly underestimated" = "#E3F2FD",   # very light blue
#   "not estimated" = "grey"
# )

col_fun_classes <-  c(
  "correctly estimated"   = "snow",  # green
  "highly overestimated"  = "#D32F2F",  # dark red
  "poorly overestimated"  = "#FFCDD2",  # light red
  "highly underestimated" = "#1976D2",  # dark blue
  "poorly underestimated" = "#BBDEFB",   # light blue,
  "not estimated" = "grey"
)

col_fun_classes_correct <-  c(
  "Low"   = "goldenrod4",  # green
  "Excellent" ="forestgreen",
  "Medium"  = "goldenrod3",  # dark red
  "High"  = "goldenrod1",
  "not estimated" = "grey"
)

h_purity = ComplexHeatmap::Heatmap(list_purities,cluster_rows = F,cluster_columns = F,
                                top_annotation = column_ha,col=col_fun_classes,
                                left_annotation = row_ha,show_column_names = F,
                                show_row_names = F,rect_gp = gpar(col = "black", lwd = 1),
                                column_split = tools,
                                row_title = "Purity",
                                row_title_side = "right",
                                row_title_gp = gpar(fontsize = 12, lineheight = 0.8,fontface="bold"),
                                row_title_rot = 90,
                                name = "Purity/Ploidy estimation classes"
)
# col_fun_ploidy = circlize::colorRamp2(c(-3,0, 3), c("forestgreen","white", "darkorange"))
h_ploidy = ComplexHeatmap::Heatmap(list_ploidy,cluster_rows = F,cluster_columns = F,
                                   # top_annotation = column_ha,
                                   left_annotation = row_ha,
                                   col=col_fun_classes, show_column_names = F,
                                   show_row_names = F,rect_gp = gpar(col = "black", lwd = 1),
                                   column_split = tools,
                                   row_title = "Ploidy",
                                   row_title_side = "right",
                                   row_title_gp = gpar(fontsize = 12, lineheight = 0.8,fontface="bold"),
                                   row_title_rot = 90,
                                   show_heatmap_legend = F
)
#col_fun_correctness_clonal = circlize::colorRamp2(c(0, 1), c("white", "goldenrod2"))
h_correc = ComplexHeatmap::Heatmap(list_correctness_clonal,cluster_rows = F,cluster_columns = F,
                                   # top_annotation = column_ha,
                                   left_annotation = row_ha,
                                   col=col_fun_classes_correct, show_column_names = F,
                                   show_row_names = F,rect_gp = gpar(col = "black", lwd = 1),
                                   column_split = tools,
                                   name = "Correclty inferred clonal CNAs"
)

col_fun_correctness_clonal = circlize::colorRamp2(c(0, 1), c("white", "#D64F54"))
h_precision_bp = ComplexHeatmap::Heatmap(list_bp_precision,cluster_rows = F,cluster_columns = F,
                                   # top_annotation = column_ha,
                                   left_annotation = row_ha,
                                   col=col_fun_correctness_clonal,
                                   show_column_names = F,
                                   show_row_names = F,rect_gp = gpar(col = "black", lwd = 1),
                                   column_split = tools,
                                   name = "Precision break point"
)
h_recall_bp = ComplexHeatmap::Heatmap(list_bp_recall,cluster_rows = F,cluster_columns = F,
                                         # top_annotation = column_ha,
                                          bottom_annotation = column_bottom_ha,
                                         left_annotation = row_ha,
                                         col=col_fun_correctness_clonal,
                                         show_column_names = F,
                                         show_row_names = F,rect_gp = gpar(col = "black", lwd = 1),
                                         column_split = tools,
                                         name = "Recall break point"
)
h_final_cna <- h_purity %v% h_ploidy %v% h_correc %v% h_precision_bp %v% h_recall_bp
pdf("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/Final_CNA_SCOUT_Validation.pdf",width = 15,height = 10)
draw(object = h_final_cna)
dev.off()
