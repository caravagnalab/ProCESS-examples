library(ProCESS)
library(ggplot2)
library(dplyr)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
source("../../getters/process_getters.R")
source("utils/plot_utils.R")
# source("compute_FGA.R")
scout_dir <-"/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
SPNS <- c("SPN01","SPN02","SPN03","SPN04","SPN06","SPN07")
COVERAGES <- c("50","100")
PURITIES <- c("0.3","0.6","0.9")

params_grid = expand.grid(COVERAGES, PURITIES)
colnames(params_grid) = c("coverage", "purity")




####
MUT_TYPES = c("INDEL", "SNV")

params_grid = expand.grid(COVERAGES, PURITIES, MUT_TYPES)
colnames(params_grid) = c("coverage", "purity", "mut")

input_dir <-  paste0(scout_dir,spn_id,"/validation/somatic/")

i = 1
df_all_SPN <- list()
for (SPN in SPNS){
  
  spn_id = SPN
  input_dir <-  paste0(scout_dir,spn_id,"/validation/somatic/")
  df <- lapply(1:nrow(params_grid), function(i) {
    coverage <- params_grid[i, ]$coverage
    purity <- params_grid[i, ]$purity
    mut_type <- params_grid[i, ]$mut
    combination <- paste0(coverage, "x_", purity, "p")
    
    results_folder_path <- file.path(input_dir, spn_id, combination, "allCaller", mut_type)
    file_name <- file.path(results_folder_path, "metrics.rds")
    
    if (file.exists(file_name)) {
      metrics <- readRDS(file_name)
      
      parsed_metrics <- lapply(names(metrics), function(caller_name) {
        caller_metrics <- metrics[[caller_name]]
        
        lapply(names(caller_metrics), function(sample_id) {
          overall <- caller_metrics[[sample_id]]$overall_metrics
          
          data.frame(
            sample = sample_id,
            precision = overall$precision,
            sensitivity = overall$sensitivity,
            specificity = overall$specificity,
            fpr = overall$fpr,
            caller = caller_name,
            coverage = coverage,
            purity = purity,
            mut_type = mut_type,
            spn = spn_id,
            stringsAsFactors = FALSE
          )
        }) %>% bind_rows()
      }) %>% bind_rows()
      
      return(parsed_metrics)
    } else {
      message("File not found: ", file_name)
      return(NULL)
    }
  }) %>% bind_rows()
  df_all_SPN[[spn_id]] <- df
}
df_all_combs_SPN <- do.call("rbind",df_all_SPN)







######## heatmap purity ########
#df_all_combs_SPN <- do.call("rbind",df_all_SPN)
df_all_combs_SPN<- df_all_combs_SPN %>% 
  arrange(caller) %>% 
  mutate(id=paste(sample,caller,sep=":")) %>% 
  mutate(comb=paste(coverage,purity,mut_type,sep=":"))

## precision matrix
list_precision <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, precision) %>%
  tidyr::pivot_wider(names_from = id, values_from = precision) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()
## sensitivity matrix


list_sensitivity <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, sensitivity) %>%
  tidyr::pivot_wider(names_from = id, values_from = sensitivity) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()

list_specificity <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, specificity) %>%
  tidyr::pivot_wider(names_from = id, values_from = specificity) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()

list_fpr <- df_all_combs_SPN %>% 
  dplyr::select(comb, id, fpr) %>%
  tidyr::pivot_wider(names_from = id, values_from = fpr) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()

###### get annotation for heatmap ######
coverages <- as.numeric(sapply(strsplit(rownames(list_sensitivity), ":"), `[`, 1))
col_coverages <- circlize::colorRamp2(
  c(0, 50, 100),
  hcl.colors(3, "Purp", rev = TRUE)  # reverse = FALSE to go from light to dark
)
purities <- as.numeric(sapply(strsplit(rownames(list_sensitivity), ":"), `[`, 2))

col_purities <- circlize::colorRamp2(
  c(0, 0.3, 0.6,0.9),
  hcl.colors(4, "Blues", rev = TRUE)  # reverse = FALSE to go from light to dark
)

mut_types <- sapply(strsplit(rownames(list_sensitivity), ":"), `[`, 3)
col_mut_types <-c("SNV"="black", "INDEL"="grey")


samples <- sapply(strsplit(colnames(list_sensitivity), ":"), `[`, 1)
tools <- sapply(strsplit(colnames(list_sensitivity), ":"), `[`, 2)
col_tools <- method_colors

spn_ids <- sapply(strsplit(colnames(list_sensitivity), "_"), `[`, 1)
#unique_spns <- unique(spn_ids)
col_spns <-c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")
names(col_spns) <- SPNS


column_ha <- HeatmapAnnotation(
  spn = spn_ids, 
  col = list(spn=col_spns)
)


row_ha <- rowAnnotation(
  coverage = coverages,
  purity = purities,
  mut_type =mut_types,
  #coverage = anno_simple(coverages, col= col_coverages),
  col = list(coverage=col_coverages, purity=col_purities,mut_type=col_mut_types)
)
####
col_fun = circlize::colorRamp2(c(0, 1), c("white", "darkorange"))
h_sen = ComplexHeatmap::Heatmap(list_sensitivity,cluster_rows = F,cluster_columns = F,
                                   top_annotation = column_ha,col=col_fun, left_annotation = row_ha,show_column_names = F,
                                   show_row_names = F,rect_gp = gpar(col = "black", lwd = 1),
                                    column_split = tools,
                                   name = "Caller sensitivity"
)
col_fun = circlize::colorRamp2(c(0, 1), c("white", "forestgreen"))
h_prec = ComplexHeatmap::Heatmap(list_precision,cluster_rows = F,cluster_columns = F,
                                # top_annotation = column_ha,
                                col=col_fun, 
                                left_annotation = row_ha,
                                show_column_names = F,
                                show_row_names = F,rect_gp = gpar(col = "black", lwd = 1),
                                column_split = tools,
                                name = "Caller precision"
)
col_fun = circlize::colorRamp2(c(0, 1), c("white", "goldenrod"))
h_spec = ComplexHeatmap::Heatmap(list_specificity,cluster_rows = F,cluster_columns = F,
                                 # top_annotation = column_ha,
                                 col=col_fun, 
                                 left_annotation = row_ha,
                                 show_column_names = F,
                                 show_row_names = F,rect_gp = gpar(col = "black", lwd = 1),
                                 column_split = tools,
                                 name = "Caller specificity"
)
col_fun = circlize::colorRamp2(c(0, 1), c("white", "purple"))
h_fpr = ComplexHeatmap::Heatmap(list_fpr,cluster_rows = F,cluster_columns = F,
                                 # top_annotation = column_ha,
                                col=col_fun, left_annotation = row_ha,show_column_names = F,
                                 show_row_names = F,rect_gp = gpar(col = "black", lwd = 1),
                                 column_split = tools,
                                 name = "Caller FPR"
)
h_final_somatic <- h_sen %v% h_prec %v% h_spec %v%h_fpr
pdf("Final_Somatic_SCOUT_Validation.pdf",width = 20,height = 10)
draw(object = h_final_somatic)
dev.off()


