library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)

purities <- c(0.3,0.6,0.9)
coverages <- c(50,100)
contexts <- c("SBS96","ID83")
all_combs <- list()
all_metrics <- list()
all_cosine  <- list()
SPNS <- c("SPN01","SPN02","SPN03","SPN04", "SPN06", 'SPN07')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
validation_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION/"

cna_caller <- "ascat"
variant_caller <- "mutect2"
color_caller <- list("ID83"=c('BASCULE' = '#E1BFF8','SigProfiler' = '#F6AF92'),
                     "SBS96"=c('BASCULE' = 'purple1','SigProfiler' = 'orange2'))
mut_counts_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/mutations_counts/"
params_grid = expand.grid(coverages, purities)
colnames(params_grid) = c("coverage", "purity")
df_all_SPN_muts_count <- list()
for (SPN in SPNS){
  spn <- SPN
  df_all_SPN_muts_count[[spn]] = lapply(1:nrow(params_grid), function(i) {
    coverage = params_grid[i,]$coverage
    purity = params_grid[i,]$purity
    comb <- paste0(coverage,"x_",purity,"p")
    mut_counts_dir_spn <- file.path(mut_counts_dir,spn)
    all_metrics_counts <- list()
    df <- readRDS(file.path(mut_counts_dir_spn,paste0("mutations_counts_",comb,".rds")))
  }) %>% bind_rows()
}
df_all_SPN_muts_count <- do.call("rbind",df_all_SPN_muts_count) %>% 
  mutate(context=case_when(type=="SNV"~"SBS96",
                           TRUE ~ "ID83"))
df_all_combs_SPN_signatures <- lapply(SPNS,function(spn){
  validation_dir_spn <- paste0(validation_dir,spn,"/signature/")
  all_combs <- list()
  all_metrics <- list()
  all_cosine  <- list()
  for (purity in purities){
    for (coverage in coverages){
      comb <- paste0(coverage,"x_",purity)
      metrics_spn <- list()
      cosine_spn  <- list()
      comb_caller <- paste0(variant_caller,"_",cna_caller)
      for (ctx in contexts){
        rds_metrics <- file.path(validation_dir_spn,comb,comb_caller,paste0("metrics_",ctx,"_spn.rds"))
        rds_cosine  <- file.path(validation_dir_spn,comb,comb_caller,paste0("cosine_mse_",ctx,".rds"))
        
        metrics_spn[[ctx]] <- readRDS(rds_metrics) %>% 
          dplyr::mutate(context=ctx,
                        purity=purity,
                        coverage=coverage)
        
        cosine_spn[[ctx]] <- readRDS(rds_cosine) %>% 
          dplyr::mutate(context=ctx,
                        purity=purity,
                        coverage=coverage)
      }
      
      # bind SBS96 + ID83 for this purity/coverage
      metrics_df <- dplyr::bind_rows(metrics_spn)
      cosine_df  <- dplyr::bind_rows(cosine_spn)
      
      # save into global list
      all_metrics[[comb]] <- metrics_df
      all_cosine[[comb]]  <- cosine_df
    }
  }
  cosine_metrics <- all_cosine %>% bind_rows()
  #metrics <- all_metrics %>% bind_rows()
  #all <- list("cosine"=cosine_metrics,"general_metrics"=metrics)
}) %>% bind_rows()

df_all_combs_SPN_signatures <-df_all_combs_SPN_signatures %>% 
  mutate(spn=case_when(spn=="SPN07_last"~"SPN07",
                       TRUE ~spn)) %>% 
  mutate(sample = str_replace_all(sample, "last_", "")) %>% 
  left_join(df_all_SPN_muts_count)

df_all_combs_SPN_signatures <- df_all_combs_SPN_signatures %>% 
  mutate(mutation_class=case_when((mutation_count>=1e4 & type=="INDEL") | (mutation_count>=1e5 & type=="SNV")~"high TMB",
                                  TRUE ~ "low TMB"))

df_all_combs_SPN_signatures<- df_all_combs_SPN_signatures %>% 
  # dplyr::filter(type=="SNV") %>% 
  mutate(id=paste(sample,caller,sep=":")) %>% 
  mutate(comb=paste(coverage,purity,context,sep=":")) %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0)))

list_cosine <- df_all_combs_SPN_signatures %>% 
  # dplyr::filter(type=="SNV") %>% 
  dplyr::select(comb, id, cosine) %>%
  tidyr::pivot_wider(names_from = id, values_from = cosine) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()

list_mse <- df_all_combs_SPN_signatures %>% 
  # dplyr::filter(type=="SNV") %>% 
  dplyr::select(comb, id, mse) %>%
  tidyr::pivot_wider(names_from = id, values_from = mse) %>% 
  tibble::column_to_rownames("comb") %>% 
  as.matrix()


###### get annotation for heatmap ######


samples <- sapply(strsplit(colnames(list_cosine), ":"), `[`, 1)
tools <- sapply(strsplit(colnames(list_cosine), ":"), `[`, 2)
col_tools <- c('BASCULE' = 'purple1','SigProfiler' = 'orange2')

spn_ids <- sapply(strsplit(colnames(list_cosine), "_"), `[`, 1)


muts_counts_map_snv <- df_all_combs_SPN_signatures %>% 
  dplyr::filter(type=="SNV") %>% 
  distinct(sample, mutation_count)
muts_counts_values_snv <- log10(muts_counts_map_snv$mutation_count[match(samples, muts_counts_map_snv$sample)])

muts_counts_map_indel <- df_all_combs_SPN_signatures %>% 
  dplyr::filter(type=="INDEL") %>% 
  distinct(sample, mutation_count)
muts_counts_values_indel <- log10(muts_counts_map_indel$mutation_count[match(samples, muts_counts_map_indel$sample)])

coverages <- as.numeric(sapply(strsplit(rownames(list_cosine), ":"), `[`, 1))
purities <- as.numeric(sapply(strsplit(rownames(list_cosine), ":"), `[`, 2))
contexts <- sapply(strsplit(rownames(list_cosine), ":"), `[`, 3)

column_ha <- HeatmapAnnotation(
  mutation_count_snv = anno_lines(muts_counts_values_snv,add_points = T),
  mutation_count_indel = anno_lines(muts_counts_values_indel,add_points = T),
  spn = spn_ids, 
  # tool = tools,
  col = list(spn=SPN_colors),
  annotation_label = c('log(TMB) SNV',
                       'log(TMB) INDEL',
                       "SPN")
)
row_ha <- rowAnnotation(
  coverage = coverages,
  purity = purities,
  context = contexts,
  #coverage = anno_simple(coverages, col= col_coverages),
  col = list(coverage=coverage_colors, purity=purity_colors,context=c("SBS96"="grey40","ID83"="grey91")),show_annotation_name = F
)
col_fun = circlize::colorRamp2(c(min(list_cosine,na.rm=TRUE), 
                                 max(list_cosine,na.rm=TRUE)), c("steelblue4","white"))
h_cosine = ComplexHeatmap::Heatmap(list_cosine,cluster_rows = F,cluster_columns = F,
                                   top_annotation = column_ha,
                                   col=col_fun,
                                   left_annotation = row_ha,show_column_names = F,
                                   column_split = tools,
                                   row_title = "Cosine",
                                   row_title_side = "right",
                                   row_title_gp = gpar(fontsize = 12, lineheight = 0.8),
                                   row_title_rot = 0,
                                   show_row_names = F,rect_gp = gpar(col = "white", lwd = 1),
                                   name = "Cosine similarity"
)

col_fun = circlize::colorRamp2(c(min(list_mse,na.rm=TRUE), 
                                 max(list_mse,na.rm=TRUE)), c("white","mediumpurple4"))
h_mse = ComplexHeatmap::Heatmap(list_mse,cluster_rows = F,cluster_columns = F,
                                   col=col_fun,
                                   left_annotation = row_ha,show_column_names = F,
                                   column_split = tools,
                                   row_title = "MSE",
                                   row_title_side = "right",
                                   row_title_gp = gpar(fontsize = 12, lineheight = 0.8),
                                   row_title_rot = 0,
                                   show_row_names = F,rect_gp = gpar(col = "white", lwd = 1),
                                   name = "MSE"
)

h_final_signatures <- h_cosine %v% h_mse

pdf("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/supplementary/Final_Signatures_SCOUT_Validation.pdf",width =9,height = 5)
draw(
  h_final_signatures,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE   # merges multiple legends into one row if possible
)
dev.off()

