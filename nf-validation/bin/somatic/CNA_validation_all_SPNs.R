library(ProCESS)
library(ggplot2)
library(dplyr)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
source("../../getters/process_getters.R")
# source("compute_FGA.R")
scout_dir <-"/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
SPNS <- c("SPN01","SPN02","SPN03","SPN04","SPN06","SPN07")
COVERAGES <- c("50","100")
PURITIES <- c("0.3","0.6","0.9")

params_grid = expand.grid(COVERAGES, PURITIES)
colnames(params_grid) = c("coverage", "purity")

df_all_SPN <- list()
for (SPN in SPNS){
  spn = SPN
  df = lapply(1:nrow(params_grid), function(i) {
    coverage = params_grid[i,]$coverage
    purity = params_grid[i,]$purity
    comb <- paste0(coverage,"x_",purity,"p")
    samples <- get_sample_names(spn = spn)
    validation_dir_cna <- paste0(scout_dir,spn,"/validation/cna/",spn)
    all_metrics_comb <- list()
    for (sample in samples){
      metrics_filename <- paste0(validation_dir_cna,"/",comb,"/",sample,"/metrics.rds")
      if (!file.exists(metrics_filename)){
        message("File not found: ", metrics_filename)
        # all_metrics_comb[[sample]] <- NA
      } else{
        all_metrics_comb[[sample]] <- readRDS(metrics_filename) %>% 
          mutate(delta_purity=as.numeric(true_purity)-as.numeric(purity)) %>% 
          mutate(delta_ploidy=as.numeric(true_ploidy)-as.numeric(ploidy))
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
  mutate(comb=paste(coverage,true_purity,sep=":"))

## purity matrix

list_purities <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, delta_purity) %>%
  tidyr::pivot_wider(names_from = id, values_from = delta_purity) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()


## ploidy matrix


list_ploidy <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, delta_ploidy) %>%
  tidyr::pivot_wider(names_from = id, values_from = delta_ploidy) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()

## correctness clonal matrix

list_correctness_clonal <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, correctness_clonal) %>%
  tidyr::pivot_wider(names_from = id, values_from = correctness_clonal) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()


## correctness clonal matrix

list_bp_distance <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, bp_distance) %>%
  tidyr::pivot_wider(names_from = id, values_from = bp_distance) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()



###### get annotation for heatmap ######
coverages <- as.numeric(sapply(strsplit(rownames(list_purities), ":"), `[`, 1))
col_coverages <- circlize::colorRamp2(
  c(0, 50, 100),
  hcl.colors(3, "Purp", rev = TRUE)  # reverse = FALSE to go from light to dark
)
purities <- as.numeric(sapply(strsplit(rownames(list_purities), ":"), `[`, 2))

col_purities <- circlize::colorRamp2(
  c(0, 0.3, 0.6,0.9),
  hcl.colors(4, "Blues", rev = TRUE)  # reverse = FALSE to go from light to dark
)

samples <- sapply(strsplit(colnames(list_purities), ":"), `[`, 1)
tools_with_cnvkit <- sapply(strsplit(colnames(list_purities), ":"), `[`, 2)
tools <- sapply(strsplit(colnames(list_purities), ":"), `[`, 2)
col_tools <- c('ascat' = 'coral2', 'sequenza' = 'darkslategray4', 'cnvkit' = 'maroon')

spn_ids <- sapply(strsplit(colnames(list_purities), "_"), `[`, 1)
#unique_spns <- unique(spn_ids)
col_spns <-c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")
names(col_spns) <- SPNS

fga_map <- df_all_combs_SPN %>% distinct(sample, fga)
fga_values <- fga_map$fga[match(samples, fga_map$sample)]
fgs_map <- df_all_combs_SPN %>% distinct(sample, fgs)
fgs_values <- fgs_map$fgs[match(samples, fgs_map$sample)]


column_ha <- HeatmapAnnotation(
  fga = anno_barplot(fga_values),
  fgs = anno_barplot(fgs_values),
  #coverage = anno_simple(coverages, col= col_coverages),
  spn = spn_ids, 
  # tool = tools, 
  col = list(spn=col_spns, tool=col_tools)
)


row_ha <- rowAnnotation(
  coverage = coverages,
  purity = purities,
  #coverage = anno_simple(coverages, col= col_coverages),
  col = list(coverage=col_coverages, purity=col_purities)
)
####
col_fun_purity = circlize::colorRamp2(c(-1,0, 1), c("forestgreen","white", "darkorange"))
h_purity = ComplexHeatmap::Heatmap(list_purities,cluster_rows = F,cluster_columns = F,
                                top_annotation = column_ha,col=col_fun_purity, left_annotation = row_ha,show_column_names = F,
                                show_row_names = F,rect_gp = gpar(col = "black", lwd = 1),
                                column_split = tools,
                                name = "Delta purity"
)
col_fun_ploidy = circlize::colorRamp2(c(-3,0, 3), c("forestgreen","white", "darkorange"))
h_ploidy = ComplexHeatmap::Heatmap(list_ploidy,cluster_rows = F,cluster_columns = F,
                                   # top_annotation = column_ha,
                                   left_annotation = row_ha,
                                   col=col_fun_ploidy, show_column_names = F,
                                   show_row_names = F,rect_gp = gpar(col = "black", lwd = 1),
                                   column_split = tools,
                                   name = "Delta ploidy"
)
col_fun_correctness_clonal = circlize::colorRamp2(c(0, 1), c("white", "red"))
h_correc = ComplexHeatmap::Heatmap(list_correctness_clonal,cluster_rows = F,cluster_columns = F,
                                   # top_annotation = column_ha,
                                   left_annotation = row_ha,
                                   col=col_fun_correctness_clonal, show_column_names = F,
                                   show_row_names = F,rect_gp = gpar(col = "black", lwd = 1),
                                   column_split = tools,
                                   name = "Correclty inferred clonal CNAs"
)

#col_fun_correctness_clonal = circlize::colorRamp2(c(0, 1), c("white", "red"))
h_bp_distance = ComplexHeatmap::Heatmap(list_bp_distance,cluster_rows = F,cluster_columns = F,
                                   # top_annotation = column_ha,
                                   left_annotation = row_ha,
                                   show_column_names = F,
                                   show_row_names = F,rect_gp = gpar(col = "black", lwd = 1),
                                   column_split = tools,
                                   name = "Break point distance"
)
h_final_cna <- h_purity %v% h_ploidy %v% h_correc %v% h_bp_distance
pdf("Final_CNA_SCOUT_Validation.pdf",width = 20,height = 8)
draw(object = h_final_cna)
dev.off()
