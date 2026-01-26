rm(list=ls())
library(tidyverse)
library(ComplexHeatmap)

setwd("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/oncoprint")
source('oncoprint_function.R')
source('get_fga_fgs.R')

spns_details <- readRDS("spns_all_info.rds")
spns_details_spn05 <- readRDS("spns_all_info_05.rds")
spns_details <- c(spns_details, list(spns_details_spn05))

true_ploidy<- tibble(
  sample = c(
    "SPN01_1.1", "SPN01_1.2", "SPN01_1.3",
    "SPN02_1.1", "SPN02_1.2",
    "SPN03_1.1", "SPN03_2.1", "SPN03_3.1", "SPN03_4.1",
    "SPN04_1.1", "SPN04_2.1",
    "SPN05_1.1", "SPN05_1.2", "SPN05_1.3",
    "SPN06_1.1", "SPN06_1.2", "SPN06_2.1", "SPN06_3.1","SPN06_3.2",
    "SPN07_1.1", "SPN07_1.2", "SPN07_1.3", "SPN07_2.1", "SPN07_2.2"
  ),
  true_ploidy = c(
    2.97, 2.06, 4.12,
    1.99, 2.04,
    1.93, 1.92, 1.91, 1.92,
    1.96, 1.96,
    2.02, 2.04, 2.01, 
    1.96, 1.95, 1.94, 1.95, 3,
    2.53, 2.53, 2.53, 2.53, 2.53
  )
)

names(spns_details) = lapply(spns_details, function(s) {
  s$cna %>% 
    names %>% 
    lapply(., function(x) {str_split_1(x ,'_')[1]}) %>% 
    unlist %>% 
    unique
  })

snv_driver = lapply(spns_details %>% names, function(s) {
  spns_details[[s]]$forest_details$drivers %>%
    dplyr::filter(type == 'SID') %>% 
    dplyr::full_join(., spns_details[[s]]$drivers_snvs) %>% 
    dplyr::mutate(category = ifelse(str_length(ref) != 1 |  str_length(alt) != 1, 'SNV', 'SNV'))
})
names(snv_driver) = names(spns_details)

# map drivers on karyotypes
# 1. mutations in cnaqc format

snv_driver_pivoted = lapply(snv_driver, function(x) {
  vaf = x %>% 
    tidyr::pivot_longer(cols = ends_with('VAF'), names_to = 'sample', values_to = 'VAF') %>% 
    dplyr::mutate(sample = gsub('.VAF', '', sample)) %>% 
    dplyr::select(mutant, type, chr, start, end, ref, alt, code, VAF, sample, category)
  
  dp = x %>% 
    tidyr::pivot_longer(cols = ends_with('coverage'), names_to = 'sample', values_to = 'DP') %>% 
    dplyr::mutate(sample = gsub('.coverage', '', sample)) %>% 
    dplyr::select(mutant, type, chr, start, end, ref, alt, code, DP, sample, category)
  
  nv = x %>% 
    tidyr::pivot_longer(cols = ends_with('occurrences'), names_to = 'sample', values_to = 'NV') %>% 
    dplyr::mutate(sample = gsub('.occurrences', '', sample)) %>% 
    dplyr::select(mutant, type, chr, start, end, ref, alt, code, NV, sample, category)
  
  res = full_join(vaf, dp) %>% 
    full_join(., nv) %>% 
    group_by(sample) %>% 
    group_split()
  names(res) = lapply(res, function(s) {s$sample %>% unique}) %>% unlist
  return(res)
})

driver_karyotypes = lapply(names(snv_driver_pivoted), function(x) {

  spn_dr = snv_driver_pivoted[[x]]
  spn_cna = spns_details[[x]]$cna

  karyo_on_muts = lapply(names(spn_dr), function(s) {

    cna = spn_cna[[s]] %>%
      dplyr::rename(segment_start = begin) %>%
      dplyr::rename(segment_to = end) %>%
      dplyr::rename(Major = major)

    snv = spn_dr[[s]]
    
    res = apply(snv, 1, function(p) {
      cna %>%
        dplyr::filter(chr == p[3]) %>%
        dplyr::filter(segment_start < as.numeric(p[4]), segment_to > as.numeric(p[5])) %>% 
        dplyr::full_join(., as_tibble_row(p))
    }) %>% 
      bind_rows()
    
    return(res)

  })
  names(karyo_on_muts) = names(spn_dr)
  return(karyo_on_muts)
})
names(driver_karyotypes) = names(snv_driver_pivoted)

# copy number drivers
cna_driver_pos = lapply(spns_details, function(s) {
  s$forest_details$drivers %>% 
    dplyr::filter(type == 'CNA')
})

cna_driver = lapply(spns_details %>% names, function(s) {
  spns_details[[s]]$cna %>% 
    bind_rows() %>% 
    dplyr::right_join(., cna_driver_pos[[s]], by = join_by(
      'chr' == 'chr', 
      'begin' == 'start'
    ))
})
names(cna_driver) = names(spns_details)

cna_driver = cna_driver %>% 
  bind_rows() %>% 
  select(chr, begin, end.x,  major, minor, ratio, sample, mutant, order, CNA_type) %>% 
  dplyr::rename(end = end.x) %>% mutate(chr =as.integer(chr))
  

# information on the last driver for spn07 is missing, edit it manually -- pyou usare quello editato 
cna_v1 = read.delim('cna_driver.csv', sep = ',')
# cna_v2 = read.delim('cna_driver_v2.csv', sep = ',')
cna_v1 = cna_v1[,1:11] 
cna_v1 = cna_driver %>% 
  left_join(cna_v1, by = join_by(sample, mutant, major, minor), suffix = c('', '.X')) %>% 
  dplyr::select(-ends_with('X')) %>% 
  mutate(code = ifelse(chr == 13 & begin == 31315086, 'BRCA2', code)) %>% 
  mutate(code = ifelse(chr == 17 & begin == 6661779, 'TP53', code)) %>% 
  mutate(code = ifelse(chr == 9 & begin == 20967752, 'CDKN2A', code)) %>% 
  filter(ratio > 0.1)

driver_karyotypes = lapply(driver_karyotypes, bind_rows) %>% 
  bind_rows()

cna_v1 = cna_v1 %>% 
  dplyr::rename(segment_start = begin) %>%
  dplyr::rename(segment_to = end) %>%
  dplyr::rename(Major = major) %>%
  dplyr::rename(type = CNA_type) %>%
  dplyr::select(-order) %>%
  dplyr::mutate(start = NA, end = NA, ref = NA, alt = NA, VAF = NA, DP = NA, NV = NA) %>%
  dplyr::mutate(chr = as.character(chr)) %>% 
  dplyr::mutate(driver = 'CNA')

all_mut_info = driver_karyotypes %>% 
  dplyr::mutate(driver = category) %>% 
  bind_rows(., cna_v1)
all_mut_info = all_mut_info %>% 
  arrange(sample) %>% 
  as_tibble() %>% 
  dplyr::mutate(NV = as.numeric(NV)) %>% 
  mutate(mut_status = ifelse(NV != 0 & !is.na(NV), 'mutated', NA))

 all_mut_info = all_mut_info %>% 
  mutate(code = ifelse(is.na(code), 'MSH6 R361H', code)) 
 
all_mut_info = all_mut_info %>% 
  dplyr::mutate(karyotype = paste(Major, minor, sep = '-')) %>% 
  dplyr::select(sample, code, karyotype, mut_status, ratio, driver) %>% 
  dplyr::group_by(sample, code) %>% 
  dplyr::mutate(subclone = seq(1:length(sample))) %>% 
  dplyr::mutate(info = paste(karyotype, mut_status, ratio, driver, sep = ';')) %>% 
  tidyr::pivot_wider(names_from = c(sample, subclone), names_sep = '_', values_from = info) 


all_mut_info <- all_mut_info %>%
  dplyr::mutate(code = str_remove(code, " .*$"))  # Remove everything after first space

# Step 1: Get long format
df_long <- all_mut_info %>%
  pivot_longer(
    cols = starts_with("SPN"),
    names_to = "sample",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

# Step 2: Extract prefix (everything before the last underscore)
df_long <- df_long %>%
  mutate(sample = str_remove(sample, "_[^_]+$")) %>% 
  mutate(driver = ifelse(code == 'PIK3CA' & karyotype == '2-2', 'CNA', driver))

# Step 3: Collapse replicates by prefix
df_wide <- df_long %>%
  group_by(code, sample, driver) %>%
  summarise(all_info = paste(karyotype, collapse = ","), .groups = "drop") %>%
  pivot_wider(
    names_from = sample,
    values_from = all_info
  ) %>% 
  arrange(driver) 
 #mutate(driver = ifelse(grepl('\\*',code), 'Driver INDEL', driver))

karyo_col = CNAqc:::get_karyotypes_colors(karyotypes = all_mut_info$karyotype %>% unique)
names(karyo_col) = gsub(':', '-', karyo_col %>% names)  
karyo_col = karyo_col[all_mut_info$karyotype %>% unique]
karyo_col[6] = '#EBD9D1'
names(karyo_col)[6] = '4-2'


df_wide = df_wide %>% 
  dplyr::group_by(code) %>%
  dplyr::summarise(
    driver = paste(sort(unique(driver)), collapse = ","),
    dplyr::across(starts_with("SPN"), ~ {
      vals <- na.omit(.)
      if(length(vals) > 0) dplyr::first(vals) else NA_character_
    }),
    .groups = "drop"
  ) %>% 
  tibble::column_to_rownames('code') 


ann_data = read.delim('top_ann_info.csv', sep = ',') %>% 
  bind_rows(c(Sample = 'SPN05_1.1', SPN = 'SPN05', Tumour_type = 'BRCA', Treatment = 'No_treatment', Gender = 'XX', WGD = 'No', Hypermutant = 'No', Clonal_status = 'Monoclonal')) %>% 
  bind_rows(c(Sample = 'SPN05_1.2', SPN = 'SPN05', Tumour_type = 'BRCA', Treatment = 'No_treatment', Gender = 'XX', WGD = 'No', Hypermutant = 'No', Clonal_status = 'Polyclonal')) %>% 
  bind_rows(c(Sample = 'SPN05_1.3', SPN = 'SPN05', Tumour_type = 'BRCA', Treatment = 'No_treatment', Gender = 'XX', WGD = 'No', Hypermutant = 'No', Clonal_status = 'Monoclonal')) %>% 
  arrange(SPN) %>% 
  dplyr::rename(Samples = Clonal_status)

# p2 = c("#7F58A5", "#D87F27", "#DFBB24", "#67ADA9", "#B5CDE3")
# p1 = list("A"="#35978F", "B"="#EF6548",
#           "C"="#E6AB02", "D"="#045A8D", "jy"="#7570B3") %>% unlist()
# p3 = c('plum4', 'goldenrod2', 'darkseagreen4', 'coral3', 'cadetblue')
# 

sample_ann = ann_data %>%
  tidyr::separate(Sample, into = c('SPN', 'Sample'), sep = '_') %>%
  dplyr::select(Sample)
ha = HeatmapAnnotation(sample = anno_text(sample_ann$Sample, rot = -0, just = "center", gp = gpar(fontsize = 10), location = unit(2, 'mm')))

# top_ann_data = ann_data %>% 
#   dplyr::select(Tumour_type, Gender, Treatment) %>% 
#   dplyr::mutate(Treatment = gsub('_I|_II', '', Treatment)) %>% 
#   dplyr::mutate(Treatment = factor(Treatment, levels = c("No_treatment", "Before_treatment", "After_treatment")))
# 
# 
# top_ann_col = list(
#   Tumour_type = setNames(p3, nm = ann_data$Tumour_type %>% unique), 
#   Gender = setNames(nm = c('XX', 'XY'), 
#                     object = c('#a53860', '#748cab')),
#   Treatment = setNames(object = c("#f2d492", '#f29559', 'indianred3'), 
#                        #object = c("#FFFFB2", '#FD8D3C', '#BD0026'), 
#                        nm = top_ann_data$Treatment %>% unique)
# )
# 
# top_ann = create_annotation(x = top_ann_data, 
#                   ann_colors = top_ann_col, 
#                   position = 'column', 
#                   pos = 'left'
#                   )
# 
# #draw(top_ann, test = '')

# adding fga and fgs annotations
fga_fgs = lapply(spns_details, function(x) {
  lapply(x$cna, function(y) {
    get_fga_fgs(y)
  })
})

fga_fgs = lapply(fga_fgs, function(x) {
  
  tibble(
    spn = names(x), 
    fga = lapply(x, function(t) {t$fga}) %>% unlist, 
    fgs = lapply(x, function(t) {t$fgs}) %>% unlist
  )
  
}) %>% bind_rows() 

df <- fga_fgs %>% mutate(fga_class=case_when(fga>=15 ~ "High FGA",
                                         TRUE ~ "Low FGA")) 

names(true_ploidy$true_ploidy) = true_ploidy$sample
fga_fgs_ploidy <- fga_fgs %>% 
  left_join(true_ploidy %>% dplyr::rename(spn = sample)) %>% 
  arrange(spn) %>% 
  dplyr::select(-spn)


fga_ann = HeatmapAnnotation(FGA = anno_barplot(fga_fgs_ploidy$fga, beside = T, border = F, bar_width = 1, gp = gpar(fill = 'snow3', col = 'white'), height = unit(1, "cm")), which = 'column',
                            annotation_name_side = 'left', annotation_name_rot = 0)
fgs_ann = HeatmapAnnotation(FGS = anno_barplot(fga_fgs_ploidy$fgs, beside = T, border = F, bar_width = 1, gp = gpar(fill = 'snow3', col = 'white'), height = unit(1, "cm")), which = 'column',
                            annotation_name_side = 'left', annotation_name_rot = 0)
ploidy_ann = HeatmapAnnotation(Ploidy = anno_barplot(fga_fgs_ploidy$true_ploidy, beside = T, border = F, bar_width = 1, gp = gpar(fill = 'snow3', col = 'white'), height = unit(1, "cm")),
                               which = 'column',
                               annotation_name_side = 'left',
                               annotation_name_rot = 0)
# draw(ploidy_ann, test = '')

bottom_ann_data = ann_data %>% 
  dplyr::select(Hypermutant) %>% #, Samples) %>% WGD
  as_tibble() 

bottom_ann_col = list(
  # Samples = setNames(
  #   object = c('#856088', '#eec170'), 
  #   nm = c('Monoclonal', 'Polyclonal')
  # ), 
  # WGD = setNames(
  #   c('#555555', 'gray95'), 
  #   nm = c('Yes', 'No')
  # ),
  Hypermutant = setNames(
    c('#555555', 'gray95'), 
    nm = c('Yes', 'No')
  ))

bottom_ann = create_annotation(x = bottom_ann_data, 
                            ann_colors = bottom_ann_col, 
                            position = 'column', 
                            pos = 'left',
                            ann_name = T, 
                            show_legend = F
)
top_ann = c(ha, bottom_ann, fga_ann, fgs_ann, ploidy_ann, gap = unit(2, "mm")) #fga_ann, fgs_ann, ploidy_ann,

names <- sub(".*_", "", ann_data$Sample)
lf_spn = HeatmapAnnotation(gene = anno_text(names, rot = -0, just = "center", gp = gpar(fontsize = 10), location = unit(2, 'mm')),
                           which = 'row')
                           #pos = 'bottom')
left_ann_spn = create_annotation(x = bottom_ann_data, 
                             ann_colors = bottom_ann_col, 
                             position = 'row', 
                             pos = 'top',
                             ann_name = T, 
                             show_legend = F
)
left_ann_spn = c(lf_spn, left_ann_spn)


fga_ann_right = HeatmapAnnotation(FGA = anno_barplot(fga_fgs_ploidy$fga, beside = T, border = F, bar_width = 1, gp = gpar(fill = 'snow3', col = 'white'), 
                                  height = unit(1, "cm")), 
                                  which = 'row')
fgs_ann_right = HeatmapAnnotation(FGS = anno_barplot(fga_fgs_ploidy$fgs, beside = T, border = F, bar_width = 1, gp = gpar(fill = 'snow3', col = 'white'),
                                  height = unit(1, "cm")), 
                                  which = 'row')
ploidy_ann_right = HeatmapAnnotation(Ploidy = anno_barplot(fga_fgs_ploidy$true_ploidy, beside = T, border = F, bar_width = 1, gp = gpar(fill = 'snow3', col = 'white'), 
                                     height = unit(1, "cm")), 
                                     which = 'row')                           

rigth_ann_spn = c(fga_ann_right, fgs_ann_right, ploidy_ann_right, gap = unit(2, "mm")) 


df_wide = df_wide[,colnames(df_wide) %>% sort]
df_wide = df_wide %>% 
  tibble::rownames_to_column('genes')

# df_wide = df_wide %>% 
#   dplyr::filter(!genes %in% c('APC', 'KRAS', 'PTEN'))
df_wide = df_wide %>% 
  tibble::column_to_rownames('genes')


# signatures annotation 
# signatures = lapply(spns_details, function(x) {
#   sign = x$forest_details$signatures %>% 
#     group_by(type) %>% 
#     group_split()
#   names(sign) = lapply(sign, function(s) {s$type %>% unique}) %>% unlist
#   return(sign)
# }) 
# names(signatures) = names(spns_details)
# 
#  
# order <- colnames(df_wide)[2:ncol(df_wide)]
# 
# sbs = read.delim('sbs.csv', sep = ',')
# sbs = sbs %>% 
#   dplyr::select(-starts_with('X')) %>% 
#   tibble() %>% 
#   bind_rows(tibble(sample = 'SPN05_1.1', SBS1 =0.1, SBS17b=0.0, SBS18=0.0, SBS5=0.3, SBS88=0.0, SBS10b=0.0, SBS6=0.6, SBS9=0.0, SBS25=0.0, SBS4=0.0, SBS11=0.0, SBS3=0.0, SBS26=0.0,))  %>% 
#   bind_rows(tibble(sample = 'SPN05_1.2', SBS1 =0.1,SBS17b=0.0, SBS18=0.0, SBS5=0.3, SBS88=0.0, SBS10b=0.0, SBS6=0.6, SBS9=0.0, SBS25=0.0, SBS4=0.0, SBS11=0.0, SBS3=0.0, SBS26=0.0,))  %>% 
#   bind_rows(tibble(sample = 'SPN05_1.3', SBS1 =0.1,SBS17b=0.0, SBS18=0.0, SBS5=0.3, SBS88=0.0, SBS10b=0.0, SBS6=0.6, SBS9=0.0, SBS25=0.0, SBS4=0.0, SBS11=0.0, SBS3=0.0, SBS26=0.0,))  %>% 
#   tibble::column_to_rownames('sample') %>% 
#   t %>% 
#   as.data.frame() %>% 
#   dplyr::select(all_of(order))
# 
# 
# sbs_colors = setNames(
#   nm = c("SBS1",
#          "SBS17b",
#          "SBS18",
#          "SBS5",
#          "SBS88",
#          "SBS10b",
#          "SBS6",
#          "SBS9",
#          "SBS25",
#          "SBS4",
#          "SBS11",
#          "SBS3" ,
#          "SBS26"), 
#   object = c('#f1696bff', 
#              '#8fbd8cff', 
#              '#87c7d6ff', 
#              '#bac3deff', 
#              '#d7bfd9ff', 
#              '#a8a2a1ff', 
#              '#cfadb3ff', 
#              '#3c609aff', 
#              '#9a4564ff', 
#              '#fbcb5bff', 
#              '#c2b280ff', 
#              '#d47e2dff', 
#              '#5f8676ff')
# )
# 
# 
# sbs_ann = HeatmapAnnotation(SBS = anno_barplot(sbs %>% t, 
#                                                gp = gpar(fill = sbs_colors, col = sbs_colors), 
#                                                bar_width = 1,  
#                                                border = F),
#                             show_legend = T, 
#                             height = unit(1.5, 'cm'),
#                             annotation_name_side= 'left', annotation_name_rot = 0)
# 
# # ID annotation
# id = read.delim('ID.csv', sep = ',') %>% 
#   dplyr::select(-starts_with('X')) %>% 
#   tibble() %>% 
#   bind_rows(tibble(sample = 'SPN05_1.1', ID1 = 0.2, ID2 = 0.1, ID4 = 0.0, ID5 = 0.0, ID7 = 0.7, ID9 = 0.0, ID18 = 0.0, ID8 = 0.0))  %>% 
#   bind_rows(tibble(sample = 'SPN05_1.2', ID1 = 0.2, ID2 = 0.1, ID4 = 0.0, ID5 = 0.0, ID7 = 0.7, ID9 = 0.0, ID18 = 0.0, ID8 = 0.0))  %>% 
#   bind_rows(tibble(sample = 'SPN05_1.3', ID1 = 0.2, ID2 = 0.1, ID4 = 0.0, ID5 = 0.0, ID7 = 0.7, ID9 = 0.0, ID18 = 0.0, ID8 = 0.0))  %>% 
#   tibble::column_to_rownames('sample') %>% 
#   t %>% 
#   as.data.frame() %>% 
#   dplyr::select(all_of(order))
# 
# id_colors = setNames(
#   nm = c('ID1', 
#          'ID2', 
#          "ID4",  
#          "ID5",  
#          "ID7",  
#          "ID8",   
#          "ID9",  
#          "ID18"), 
#   object = c('#0c8281ff', 
#              '#f5a55fff', 
#              '#7d287eff', 
#              '#2e4f4fff', 
#              '#c4ddbcff', 
#              '#996869ff', 
#              '#daa627ff', 
#              '#bc8f8fff')
# )
# 
# id_ann = HeatmapAnnotation(ID = anno_barplot(id %>% t, 
#                                                gp = gpar(fill = id_colors, col = id_colors), 
#                                                bar_width = 1,  
#                                                border = F),
#                             show_legend = T, 
#                             height = unit(1.5, 'cm'),
#                             annotation_name_side = 'left', annotation_name_rot = 0)
# 
# 
# SBS_lgd = Legend(at = c("SBS1",
#                        "SBS17b",
#                        "SBS18",
#                        "SBS5",
#                        "SBS88",
#                        "SBS10b",
#                        "SBS6",
#                        "SBS9",
#                        "SBS25",
#                        "SBS4",
#                        "SBS11",
#                        "SBS3" ,
#                        "SBS26"), title = "SBS", legend_gp = gpar(fill = c('#f1696bff', 
#                                                                           '#8fbd8cff', 
#                                                                           '#87c7d6ff', 
#                                                                           '#bac3deff', 
#                                                                           '#d7bfd9ff', 
#                                                                           '#a8a2a1ff', 
#                                                                           '#cfadb3ff', 
#                                                                           '#3c609aff', 
#                                                                           '#9a4564ff', 
#                                                                           '#fbcb5bff', 
#                                                                           '#c2b280ff', 
#                                                                           '#d47e2dff', 
#                                                                           '#5f8676ff')), 
#                  nrow = 2)
# 
# ID_lgd = Legend(at = c('ID1', 
#                        'ID2', 
#                        "ID4",  
#                        "ID5",  
#                        "ID7",  
#                        "ID8",   
#                        "ID9",  
#                        "ID18"), 
#                 title = "ID", 
#                 legend_gp = gpar(fill = c('#0c8281ff', 
#                                           '#f5a55fff', 
#                                           '#7d287eff', 
#                                           '#2e4f4fff', 
#                                           '#c4ddbcff', 
#                                           '#996869ff', 
#                                           '#daa627ff', 
#                                           '#bc8f8fff')), 
#                 nrow =2)
# 
# #sign_legends = list(SBS_lgd, ID_lgd)

# order_sample <- true_ploidy$sample
# muts_count_rds <- readRDS('mut_counts.rds') %>%
#   mutate(mutation_count=log10(mutation_count)) %>%
#   pivot_wider(
#     id_cols = c(sample, coverage, purity),
#     names_from = type,
#     values_from = mutation_count,
#     values_fill = 0,   # replace missing values with 0
#     names_glue = "{type}_count"
#   ) %>%
#   select(!c("coverage","purity")) %>% 
#   bind_rows(tibble(sample = c('SPN05_1.1', 'SPN05_1.2', 'SPN05_1.3'), INDEL_count=c(0,0,0), SNV_count=c(0,0,0))) %>% 
#   arrange(factor(sample, levels = order_sample))
# 
#   
# samples <- muts_count_rds$sample
# indel_count_map <- muts_count_rds %>% distinct(sample, INDEL_count)
# indel_count_values <- indel_count_map$INDEL_count[match(samples, indel_count_map$sample)]
# snv_count_map <- muts_count_rds %>% distinct(sample, SNV_count)
# snv_count_values <- snv_count_map$SNV_count[match(samples, snv_count_map$sample)]
# muts_count_matrix<- rbind(indel_count_values, snv_count_values) %>% t()
# 
# mut_counts = HeatmapAnnotation(Counts =anno_points(muts_count_matrix, gp = gpar(col = c("seagreen","palevioletred3")), add_points = TRUE, pt_gp = gpar(col = 5:6), pch = c(16, 16), beside = T, border = F, bar_width = 1),
#                                annotation_name_side = 'left', annotation_name_rot = 0)
# #bottom_ann = c(mut_counts, sbs_ann, id_ann,  gap = unit(1.5, "mm"))
# #draw(bottom_ann, test = '')

oncogene <- readRDS('oncogenes.rds')
supp <- readRDS('suppressor.rds')

drivers_ann_data = all_mut_info %>%
  dplyr::select(driver) %>%
  distinct() %>%
  ungroup() %>%
  mutate(type = 'NA') %>% 
  mutate(type = ifelse(code %in% oncogene$gene, 'Oncogene', type)) %>% 
  mutate(type = ifelse(code %in% supp$gene, 'TSG', type)) %>% 
  mutate(type = ifelse(code == 'ATRX',  'TSG', type)) %>% 
  select(-driver) %>% 
  distinct() 

pathways <- readRDS('genes_pathways.rds') %>% 
  dplyr::filter(significant == TRUE) %>% 
  dplyr::select(p_value, source, term_name, evidence_codes, intersection) %>% 
  dplyr::rename(code = intersection)

drivers_ann_data <- left_join(drivers_ann_data, pathways) %>% #filter(code == '13q14.2')
   filter(source == "REAC" | is.na(source)) %>% 
   filter(term_name != 'Disease' | is.na(term_name)) %>% 
   group_by(code, type) %>% #filter(code == '13q14.2')
   filter(p_value == min(p_value) |  is.na(p_value)) %>% 
  dplyr::rename(pathway = term_name) %>% 
   group_by(code) %>% 
   slice_head()

order = df_wide %>% mutate(driver = ifelse(driver == 'SNV', 'SNV,INDEL', driver)) %>% pull(driver)
genes_order = rownames(df_wide)

df_driver = tibble(gene = rownames(df_wide), df_wide$driver)
#saveRDS(object = df_driver, 'driver_ann.rds')

tmp <- drivers_ann_data[which(drivers_ann_data$code %in% genes_order), ] %>%   
  arrange(factor(code, levels = genes_order)) %>% 
  dplyr::select(type, pathway)  %>% 
  mutate(pathway = case_when(
    pathway == "Diseases of signal transduction by growth factor receptors and second messengers" ~ "GFS",
    pathway == "Nuclear events mediated by NFE2L2" ~ "NFE2L2/NRF2",
    pathway == "Diseases of DNA repair" ~ "DNA Repair",
    pathway == "Signal Transduction" ~ "Cell Signaling",
    pathway == "Alternative Lengthening of Telomeres (ALT)" ~ "ALT"
  )) %>% 
  mutate(pathway = ifelse(is.na(pathway), 'NA', pathway))  %>% 
  ungroup() %>% 
  bind_cols(order = order) %>% 
  select(-type) %>% 
  arrange(factor(pathway, levels = c("GFS", "NFE2L2/NRF2", "Cell Signaling", "DNA Repair", "ALT", "NA")), factor(order, levels = c('CNA', 'CNA,SNV', 'SNV', 'INDEL')))

drivers_ann_data = drivers_ann_data[which(drivers_ann_data$code %in% genes_order), ] %>%   
  arrange(factor(code, levels = genes_order)) %>% 
  dplyr::select(type, pathway)  %>% 
  ungroup() %>% 
  mutate(pathway = case_when(
    pathway == "Diseases of signal transduction by growth factor receptors and second messengers" ~ "GFS",
    pathway == "Nuclear events mediated by NFE2L2" ~ "NFE2L2/NRF2",
    pathway == "Diseases of DNA repair" ~ "DNA Repair",
    pathway == "Signal Transduction" ~"Cell Signaling",
    pathway == "Alternative Lengthening of Telomeres (ALT)" ~ "ALT"
  )) %>% 
  mutate(pathway = ifelse(is.na(pathway), 'NA', pathway)) %>% 
  dplyr::rename(Gene = type, Pathway = pathway)


tmp_drivers_ann_data = drivers_ann_data %>% select(-code)
right_ann = create_annotation(x = tmp_drivers_ann_data, 
                             position = 'row',
                             ann_colors = list(
                               Gene = setNames(c( 'gainsboro','#3e5c76', '#bc4b51'),
                                               nm = tmp_drivers_ann_data$Gene %>% unique),
                               Pathway = setNames(c('gainsboro',
                                                    '#68b0ab',
                                                    '#F0AB00',
                                                    '#586ba4',
                                                    '#C50084',
                                                    '#fa7921'),
                                                  nm = tmp_drivers_ann_data$Pathway %>% unique)
                             ),
                             ann_name = T, 
                             pos = 'bottom',
                             annotation_legend_param = list(ncol = 1))


ha_spn = HeatmapAnnotation(gene = anno_text(drivers_ann_data$code, rot = 45, just = "center", gp = gpar(fontsize = 10), location = unit(5.5, 'mm')))
top_ann_spn = create_annotation(x = tmp_drivers_ann_data, 
                              position = 'column',
                              ann_colors = list(
                                Gene = setNames(c( 'gainsboro','palegreen4', 'hotpink4'),
                                                nm = tmp_drivers_ann_data$Gene %>% unique),
                                Pathway = setNames(c('gainsboro',
                                                     '#68b0ab',
                                                     '#F0AB00',
                                                     '#586ba4',
                                                     '#C50084',
                                                     '#fa7921'),
                                                   nm = tmp_drivers_ann_data$Pathway %>% unique)
                              ),
                              ann_name = T, 
                              pos = 'right',
                              annotation_legend_param = list(ncol = 1))

top_ann_spn = c(ha_spn, top_ann_spn)

names <- sub(".*_", "", ann_data$Sample)
ht_oncoprint = oncoPrint(t(df_wide %>% select(-driver)),  
               alter_fun = 
                 function(x, y, w, h, v) {
                   n = sum(v)
                   h = h*0.9
                   grid.rect(x, y, w * 0.9, h, gp = gpar(fill = "gray95", col = NA))
                   
                   if(n) {grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.9, 1/n*h,
                                    gp = gpar(fill = karyo_col[names(which(v))], 
                                              karyo_col = NA,
                                              col = NA
                                    ), 
                                    just = "top")
                   }
                 }, 
               na_col = "gray95",
               col = karyo_col, 
               show_column_names = F, 
               show_row_names = F, 
               row_order = ann_data$Sample, 
               #top_annotation = top_ann,
               row_split = factor(ann_data$SPN),
               show_pct = F, 
               #right_annotation = right_ann, 
               #bottom_annotation = bottom_ann, 
               right_annotation = rigth_ann_spn,
               left_annotation = left_ann_spn,
               top_annotation = top_ann_spn,
               #row_names_side = 'left', 
               heatmap_legend_param = list(direction = 'vertical', ncol = 1), 
               column_split = order,
               column_order = tmp$code,
               name = 'CN',
               #column_names_gp = 7, 
               row_names_gp = 7,
               row_title_gp = gpar(fontsize = 10),
               row_gap = unit(3, "mm"))

pdf('scout_oncoprint_ploidy.pdf', width = 7.5, height = 5.5)
draw(ht_oncoprint, heatmap_legend_side = "right", annotation_legend_side = 'right') #annotation_legend_list = sign_legends)
dev.off()

