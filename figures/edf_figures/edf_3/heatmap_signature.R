rm(list=ls())
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/oncoprint/oncoprint_function.R')
source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/oncoprint/get_fga_fgs.R')

spns_details <- readRDS("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/oncoprint/spns_all_info.rds")
spns_details_spn05 <- readRDS("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/oncoprint/spns_all_info_05.rds")
spns_details <- c(spns_details, list(spns_details_spn05))

order_sample <- c("SPN01_1.1","SPN01_1.2","SPN01_1.3","SPN02_1.1","SPN02_1.2","SPN03_1.1","SPN03_2.1","SPN03_3.1","SPN03_4.1",
                  "SPN04_1.1","SPN04_2.1","SPN05_1.1","SPN05_1.2","SPN05_1.3","SPN06_1.1","SPN06_1.2","SPN06_2.1","SPN06_3.1",
                  "SPN06_3.2","SPN07_1.1","SPN07_1.2","SPN07_1.3","SPN07_2.1","SPN07_2.2")
order_spn <- paste0('SPN0',1:7)

scout_exposure <- readRDS("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/signature_exposures/scout_exposure.rds")

sbs_df <- scout_exposure %>%
  filter(grepl("^SBS", signature)) %>%
  select(sample, signature, median_exposure, spn) %>%
  pivot_wider(
    names_from  = signature,
    values_from = median_exposure,
    values_fill = 0,
    values_fn   = median
  ) %>% 
  ungroup() %>% 
  dplyr::mutate(across(-c(sample, spn), ~ ifelse(.x == 0, 0, 1))) %>% 
  pivot_longer(
    cols = starts_with("SBS"),
    names_to = "signature",
    values_to = "exposure"
  ) %>% 
  group_by(spn,signature) %>% 
  summarize(
    exposure = as.integer(any(exposure > 0)),
    .groups = "drop"
  ) %>% mutate(type = case_when(
    signature %in% paste0('SBS', c(6,14,15,20,21,26,44)) ~ 'MMR deficency', 
    signature %in% paste0('SBS', c('10a','10b', '10c','10d',28)) ~ 'POL deficency',
    signature %in% paste0('SBS', c(3)) ~ 'HR deficency',
    signature %in% paste0('SBS', c(30, 36)) ~ 'BER deficency',
    signature %in% paste0('SBS', c(11, 25, 31, 32, 35, 86, 87, 90, 99)) ~ 'Treatment',
    signature %in% paste0('SBS', c(2, 13)) ~ 'APOBEC',
    signature %in% paste0('SBS', c(4, 29, 92)) ~ 'Tobacco',
    signature %in% paste0('SBS', c('7a', '7b', '7c', '7d', 38)) ~ 'UV',
    signature %in% paste0('SBS', c('22a', '22b')) ~ 'AA',
    signature %in% paste0('SBS', c(88)) ~ 'Colibactin',
    signature %in% paste0('SBS', c(27, 43, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 95)) ~ 'Artifact',
    signature %in% paste0('SBS', c(9, 84, 85)) ~ 'Lymphoid',
    signature %in% paste0('SBS', c(1,5)) ~ 'Clock-like',
    .default = 'Other'
  )) %>% 
  mutate(exposure = ifelse(signature == 'SBS11' & spn == 'SPN04', 1, exposure))



id_sample = scout_exposure %>%
  filter(grepl("^ID", signature)) %>%
  select(sample, signature, median_exposure, spn) %>%
  pivot_wider(
    names_from  = signature,
    values_from = median_exposure,
    values_fill = 0,
    values_fn   = median
  ) %>% 
  ungroup() %>% 
  dplyr::mutate(across(-c(sample, spn), ~ ifelse(.x == 0, 0, 1))) %>% 
  pivot_longer(
    cols = starts_with("ID"),
    names_to = "signature",
    values_to = "exposure"
  ) %>% 
  group_by(spn,signature) %>% 
  summarize(
    exposure = as.integer(any(exposure > 0)),
    .groups = "drop"
  ) %>% mutate(type = case_when(
    signature %in% paste0('SBS', c(6,14,15,20,21,26,44)) ~ 'MMR deficency', 
    signature %in% paste0('SBS', c('10a','10b', '10c','10d',28)) ~ 'POL deficency',
    signature %in% paste0('SBS', c(3)) ~ 'HR deficency',
    signature %in% paste0('SBS', c(30, 36)) ~ 'BER deficency',
    signature %in% paste0('SBS', c(11, 25, 31, 32, 35, 86, 87, 90, 99)) ~ 'Treatment',
    signature %in% paste0('SBS', c(2, 13)) ~ 'APOBEC',
    signature %in% paste0('SBS', c(4, 29, 92)) ~ 'Tobacco',
    signature %in% paste0('SBS', c('7a', '7b', '7c', '7d', 38)) ~ 'UV',
    signature %in% paste0('SBS', c('22a', '22b')) ~ 'AA',
    signature %in% paste0('SBS', c(88)) ~ 'Colibactin',
    signature %in% paste0('SBS', c(27, 43, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 95)) ~ 'Artifact',
    signature %in% paste0('SBS', c(9, 84, 85)) ~ 'Lymphoid',
    signature %in% paste0('SBS', c(1,5)) ~ 'Clock-like',
    .default = 'Other'
  )) %>% 
  mutate(type = case_when(
    signature %in% paste0('ID', c(1,2)) ~ 'Clock-like',
    signature %in% paste0('ID', c(7)) ~ 'MMR deficency',
    signature %in% paste0('ID', c(8)) ~ 'DSB repair',
    signature %in% paste0('ID', c(9,4, 5)) ~ 'Other',
    signature %in% paste0('ID', c(18)) ~ 'Colibactin',
    signature %in% paste0('ID', c(3)) ~ 'APOBEC',
    .default = 'Other'
  ))

df_spn <- sbs_df %>% 
  bind_rows(id_sample) %>% 
  mutate(exposure = ifelse(exposure =='NA', 0, exposure)) %>% 
  distinct() %>% 
  group_by(spn, signature) %>%
  filter(exposure == max(exposure)) %>%
  #slice(1) %>%
  ungroup()


mat <- df_spn %>%
  select(spn, signature, exposure) %>%
  tidyr::pivot_wider(names_from = spn, values_from = exposure) %>%
  tibble::column_to_rownames("signature") %>%
  as.matrix()

# Convert to numeric
mat <- apply(mat, 2, as.numeric)
rownames(mat) <- df_spn %>% 
  distinct(signature) %>% 
  pull(signature)

# Get signature type annotation
sig_types <- df_spn %>%
  distinct(signature, type) %>%
  tibble::column_to_rownames("signature")

# Define custom order for signature types
type_order <- c("Clock-like", "Treatment", "MMR deficency", 
                "POL deficency", "HR deficency", "DSB repair", 'APOBEC', "Tobacco", "Colibactin", 
                "Lymphoid",   "Other")

# Convert type to factor with custom order
sig_types$type <- factor(sig_types$type, levels = type_order)

# Define colors for signature types
type_colors <- c(
  "Clock-like" = "deepskyblue4",
  "Other" = "lightblue4",
  "Colibactin" = "palegreen4",
  "POL deficency" = "plum4",
  "MMR deficency" = "sienna3",
  "Lymphoid" = "goldenrod",
  "Treatment" = "darkorange",
  "Tobacco" = "palevioletred",
  "HR deficency" = "tomato1",
  "DSB repair" = 'steelblue2',
  "APOBEC" = 'lightpink'
)

# # Create row annotation
# row_ha <- topAnnotation(
#   Type = sig_types$type,
#   col = list(Type = type_colors),
#   annotation_name_side = "top",
#   show_annotation_name = FALSE,
#   annotation_legend_param = list(
#     Type = list(title = "Signature Type")
#   )
# )
sig_types <- sig_types %>% dplyr::rename(Type = type)
ann <- ComplexHeatmap::HeatmapAnnotation(df = sig_types, 
                                         col = list(Type = type_colors), 
                                         which = 'column',
                                         show_annotation_name = F,                                 
                                         na_col = '#EEEEEE',   annotation_legend_param = list(
                                           title_gp  = gpar(fontsize = 10),
                                           labels_gp = gpar(fontsize = 9),
                                           grid_height = unit(0.25, "cm"),
                                           grid_width  = unit(0.25, "cm")
                                         ))
draw(ann, test = '')


# Create heatmap
ht <- Heatmap(
  mat %>% t(),
  name = "Signature",
  col = c("0" = "white", "1" = "gainsboro"),
  show_heatmap_legend = F,
  
  #Row settings - no clustering, split by type
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9),
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  column_split = sig_types$Type,
  row_title_rot = 0,
  row_title = NULL,
  
  # Column settings - no clustering, names on top
  column_title = NULL,
  column_names_side = "top",
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  #column_names_centered = T,
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  row_split = colnames(mat),
  
  column_title_gp = gpar(fontsize = 9),
  row_title_gp = gpar(fontsize = 9),
  
  # Add row annotation
  top_annotation = ann,
  
  # Cell settings
  rect_gp = gpar(col = "gray90", lwd = 0.8),
  
  # Spacing settings
  row_gap = unit(1, "mm"),  # Space between cells within same group
  column_gap = unit(2, "mm"),
  border = F,  # Border around heatmap
  
  # Legend
  heatmap_legend_param = list(
    title = "Signature",
    at = c(0, 1),
    labels = c("Absent", "Present"),
    legend_direction = "vertical",
    title_gp  = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 9),
    grid_height = unit(0.3, "cm"),
    grid_width  = unit(0.3, "cm")
  )
)
#saveRDS(object = ht, file = '~/Desktop/signature.rds')

pdf('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_3/panela_signature.pdf', width = 8, height = 3)
# Draw the heatmap


draw(ht, heatmap_legend_side = "right", 
     annotation_legend_side = "right")
dev.off()
