rm(list=ls())
library(tidyverse)
library(ComplexHeatmap)
# setwd('oncoprint/')
source('oncoprint_function.R')

spns_details <- readRDS("spns_all_info.rds")

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
    dplyr::full_join(., spns_details[[s]]$drivers_snvs)
})
names(snv_driver) = names(spns_details)

# map drivers on karyotypes
# 1. mutations in cnaqc format

snv_driver_pivoted = lapply(snv_driver, function(x) {
  vaf = x %>% 
    tidyr::pivot_longer(cols = ends_with('VAF'), names_to = 'sample', values_to = 'VAF') %>% 
    dplyr::mutate(sample = gsub('.VAF', '', sample)) %>% 
    dplyr::select(mutant, type, chr, start, end, ref, alt, code, VAF, sample)
  
  dp = x %>% 
    tidyr::pivot_longer(cols = ends_with('coverage'), names_to = 'sample', values_to = 'DP') %>% 
    dplyr::mutate(sample = gsub('.coverage', '', sample)) %>% 
    dplyr::select(mutant, type, chr, start, end, ref, alt, code, DP, sample)
  
  nv = x %>% 
    tidyr::pivot_longer(cols = ends_with('occurrences'), names_to = 'sample', values_to = 'NV') %>% 
    dplyr::mutate(sample = gsub('.occurrences', '', sample)) %>% 
    dplyr::select(mutant, type, chr, start, end, ref, alt, code, NV, sample)
  
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
  bind_rows()

# write.table(cna_driver, 'cna_driver_v2.csv', sep = ',', quote = F, row.names = F)

# information on the last driver for spn07 is missing, edit it manually -- pyou usare quello editato 

cna_v1 = read.delim('cna_driver.csv', sep = ',')
# cna_v2 = read.delim('cna_driver_v2.csv', sep = ',')
cna_v1 = cna_v1[,1:11]

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
  dplyr::mutate(driver = 'Driver CNA')

all_mut_info = driver_karyotypes %>% 
  dplyr::mutate(driver = 'Driver SNV') %>% 
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

# prefixes you want to merge by (everything before last underscore)
# prefixes <- c("SPN01_1.1", "SPN01_1.2", "SPN01_1.3")

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
  mutate(sample = str_remove(sample, "_[^_]+$"))

# Step 3: Collapse replicates by prefix
df_wide <- df_long %>%
  group_by(code, sample) %>%
  # dplyr::mutate(karyotype = case_when((karyotype == '1-1',  mut_status != 'mutated') ~ NA, .default = karyotype))
  # mutate(info = paste(karyotype, mut_status, sep = ',')) %>% 
  summarise(all_info = paste(karyotype, collapse = ","), .groups = "drop") %>%
  pivot_wider(
    names_from = sample,
    values_from = all_info
  )

karyo_col = CNAqc:::get_karyotypes_colors(karyotypes = all_mut_info$karyotype %>% unique)
names(karyo_col) = gsub(':', '-', karyo_col %>% names)  
karyo_col = karyo_col[all_mut_info$karyotype %>% unique]
karyo_col[6] = '#EBD9D1'
names(karyo_col)[6] = '4-2'

# now merge the subclones together!

# vars_fun = c(
#   'mutated' = alter_graphic("rect", width = 0.9, height = 0.9, fill = 'steelblue'), 
#   'background' = alter_graphic("rect", width = 0.9, height = 0.9, fill = "#CCCCCC")
# )

# vars_fun = list(
#   'background' = alter_graphic("rect", width = 0.9, height = 0.9, fill = "#CCCCCC"), 
#   # 'ratio' = alter_graphic('rect', width = 0.9, height = 0.9, fill = "#CCCCCC", )
#   'karyotype' = function(x, y, w, h, v) {
#     n = sum(v)  # how many alterations for current gene in current sample
#     h = h*0.9
#     # use `names(which(v))` to correctly map between `v` and `col`
#     if(n) grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.9, 1/n*h, 
#                     gp = gpar(fill = karyo_col[names(which(v))], col = NA), just = "top")}
# )

df_wide = df_wide %>% 
  tibble::column_to_rownames('code')

# this would be nice but not working atm
# all_ratios = all_mut_info$ratio %>% unique

# oncoPrint(df_wide,
#           alter_fun = list(
#             'karyotype' = function(x, y, w, h, v, t) {
#             n = sum(v)  # how many alterations for current gene in current sample
#             h = h*0.9
#             # use `names(which(v))` to correctly map between `v` and `col`
#             if(n) grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.9, 1/n*h, 
#                             gp = gpar(fill = adjustcolor(karyo_col[names(which(v))], alpha.f = t[which(t == all_ratios)]), 
#                                       karyo_col = NA), 
#                             just = "top")
#           }, 
#           'background' = alter_graphic("rect", width = 0.9, height = 0.9, fill = "#CCCCCC")
#           ), col = karyo_col)


# creating annotations
# top_ann_data = lapply(spns_details %>% names, function(s) {
#   tibble(
#     Sample = spns_details[[s]]$cna %>% names
#   ) %>% 
#     mutate(SPN = s)
#   }) %>% 
#   bind_rows
# write.table(top_ann_data, file = 'top_ann_info.csv', sep = ',', quote = F, row.names = F)

ann_data = read.delim('top_ann_info.csv', sep = ',')

# top_ann_data = top_ann_data %>% 
#   dplyr::select(SPN)
# 

top_ann_col = list(
  # SPN = setNames(c('steelblue', 'seagreen', 'goldenrod', 'coral', 'palevioletred', 'indianred3'),
  #                nm = ann_data$SPN %>% unique),
  Tumour_type = setNames(c(awtools::a_palette[1:length(ann_data$Tumour_type %>% unique)]), 
                         nm = ann_data$Tumour_type %>% unique)
)

sample_ann = ann_data %>% 
  tidyr::separate(Sample, into = c('SPN', 'Sample'), sep = '_') %>% 
  dplyr::select(Sample)
ha = HeatmapAnnotation(sample = anno_text(sample_ann$Sample, rot = -45, location = unit(2, 'mm')))

top_ann_data = ann_data %>% 
  dplyr::select(Tumour_type, Treatment, Gender) %>% 
  dplyr::mutate(Treatment = gsub('_I|_II', '', Treatment)) %>% 
  dplyr::mutate(Treatment = factor(Treatment, levels = c("No_treatment", "Before_treatment", "After_treatment")))

top_ann_col = list(
  # SPN = setNames(c('steelblue', 'seagreen', 'goldenrod', 'coral', 'palevioletred', 'indianred3'),
  #                nm = ann_data$SPN %>% unique),
  Tumour_type = setNames(c(awtools::a_palette[1:length(ann_data$Tumour_type %>% unique)]), 
                         nm = ann_data$Tumour_type %>% unique), 
  # Treatment = setNames(object = RColorBrewer::brewer.pal(length(top_ann_data$Treatment %>% unique), 'YlOrRd') , 
  #                      nm = top_ann_data$Treatment %>% unique), 
  Treatment = setNames(object = c("#FFFFB2", '#FD8D3C', '#BD0026'), 
                       nm = top_ann_data$Treatment %>% unique), 
  Gender = setNames(nm = c('XX', 'XY'), 
                    object = c('#239BA7', '#8FA31E'))
)

top_ann = create_annotation(x = top_ann_data, 
                  ann_colors = top_ann_col, 
                  position = 'column', 
                  pos = 'left'
                  )

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
  
}) %>% bind_rows() %>% 
  dplyr::select(-spn)

fga_ann = HeatmapAnnotation(FGA = anno_barplot(fga_fgs$fga, beside = T, border = F), which = 'column')
fgs_ann = HeatmapAnnotation(FGS = anno_barplot(fga_fgs$fgs, beside = T, border = F), which = 'column')

top_ann = c(ha, top_ann, fga_ann, fgs_ann, gap = unit(1.5, "mm"))

# creating bottom annotations 
#   - mut rate
#   - cna rate
#   - wgd
#   - signatures

bottom_ann_data = ann_data %>% 
  dplyr::select(WGD, Hypermutant, Clonal_status) %>% 
  dplyr::mutate(WGD = ifelse(WGD == 'No', NA, 'WGD'),
                Hypermutant = ifelse(Hypermutant == 'No', NA, 'Hypermutant')) %>%
  as_tibble() %>%
  tidyr::unite(col = 'Status', WGD, Hypermutant, na.rm = T, sep = '')

bottom_ann_col = list(
  Clonal_status = setNames(
    object = c('#643c7bff', '#daa627ff'), 
    nm = c('Monoclonal', 'Polyclonal')
  ), 
  Status = setNames(
    c('#064232', '#660B05'), 
    nm = c('WGD', 'Hypermutant')
  ))

#   WGD = setNames(c('#064232', '#DEE8CE'),
#                  nm = c('WGD', 'No')), 
#   Hypermutant = setNames(c('#', '#CCC'), 
#                          nm = c('Hypermutant', 'No'))
# )

bottom_ann = create_annotation(x = bottom_ann_data, 
                            ann_colors = bottom_ann_col, 
                            position = 'column', 
                            pos = 'left'
)

# ht = oncoPrint(mat, 
#           alter_fun = vars_fun, 
#           top_annotation = top_ann, 
#           column_split = factor(top_ann_data$SPN), 
#           show_pct = F, 
#           bottom_annotation = bottom_ann, 
#           row_names_side = 'left'
#           
#           )
# draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

muts = c('mutated' = 16)

df_wide = df_wide[,colnames(df_wide) %>% sort]

df_wide = df_wide %>% 
  tibble::rownames_to_column('genes')

df_wide = df_wide %>% 
  dplyr::filter(!genes %in% c('APC', 'KRAS', 'PTEN'))
df_wide = df_wide %>% 
  tibble::column_to_rownames('genes')


# signatures annotation 
signatures = lapply(spns_details, function(x) {
  sign = x$forest_details$signatures %>% 
    group_by(type) %>% 
    group_split()
  names(sign) = lapply(sign, function(s) {s$type %>% unique}) %>% unlist
  return(sign)
}) 
names(signatures) = names(spns_details)

sbs = read.delim('sbs.csv', sep = ',')
sbs = sbs %>% 
  tibble::column_to_rownames('sample') %>% 
  t

sbs_colors = setNames(
  nm = c("SBS1",
         "SBS17b",
         "SBS18",
         "SBS5",
         "SBS88",
         "SBS10b",
         "SBS6",
         "SBS9",
         "SBS25",
         "SBS4",
         "SBS11",
         "SBS3" ,
         "SBS26"), 
  object = c('#f1696bff', 
             '#8fbd8cff', 
             '#87c7d6ff', 
             '#bac3deff', 
             '#d7bfd9ff', 
             '#a8a2a1ff', 
             '#cfadb3ff', 
             '#3c609aff', 
             '#9a4564ff', 
             '#fbcb5bff', 
             '#c2b280ff', 
             '#d47e2dff', 
             '#5f8676ff')
)

sbs_ann = HeatmapAnnotation(SBS = anno_barplot(sbs %>% t,
                                          gp = gpar(fill = sbs_colors)
                                          ), show_legend = T, height = unit(1.5, 'cm'))

# ID annotation
id = read.delim('ID.csv', sep = ',')
id = id %>% 
  tibble::column_to_rownames('sample') %>% 
  dplyr::select(-starts_with('X')) %>% 
  t

id_colors = setNames(
  nm = c('ID1', 
         'ID2', 
         "ID4",  
         "ID5",  
         "ID7",  
         "ID8",   
         "ID9",  
         "ID18"), 
  object = c('#0c8281ff', 
             '#f5a55fff', 
             '#7d287eff', 
             '#2e4f4fff', 
             '#c4ddbcff', 
             '#996869ff', 
             '#daa627ff', 
             '#bc8f8fff')
)

id_ann = HeatmapAnnotation(ID = anno_barplot(id %>% t,
                                               gp = gpar(fill = id_colors)
), show_legend = T, height = unit(1.5, 'cm'))

SBS_lgd = Legend(at = c("SBS1",
                       "SBS17b",
                       "SBS18",
                       "SBS5",
                       "SBS88",
                       "SBS10b",
                       "SBS6",
                       "SBS9",
                       "SBS25",
                       "SBS4",
                       "SBS11",
                       "SBS3" ,
                       "SBS26"), title = "SBS", legend_gp = gpar(fill = c('#f1696bff', 
                                                                          '#8fbd8cff', 
                                                                          '#87c7d6ff', 
                                                                          '#bac3deff', 
                                                                          '#d7bfd9ff', 
                                                                          '#a8a2a1ff', 
                                                                          '#cfadb3ff', 
                                                                          '#3c609aff', 
                                                                          '#9a4564ff', 
                                                                          '#fbcb5bff', 
                                                                          '#c2b280ff', 
                                                                          '#d47e2dff', 
                                                                          '#5f8676ff')), 
                 nrow = 2)

ID_lgd = Legend(at = c('ID1', 
                       'ID2', 
                       "ID4",  
                       "ID5",  
                       "ID7",  
                       "ID8",   
                       "ID9",  
                       "ID18"), 
                title = "ID", 
                legend_gp = gpar(fill = c('#0c8281ff', 
                                          '#f5a55fff', 
                                          '#7d287eff', 
                                          '#2e4f4fff', 
                                          '#c4ddbcff', 
                                          '#996869ff', 
                                          '#daa627ff', 
                                          '#bc8f8fff')), 
                nrow =2)

sign_legends = list(SBS_lgd, ID_lgd)

bottom_ann = c(bottom_ann, sbs_ann, id_ann, gap = unit(1.5, "mm"))
# write.table(df_wide, 'df_wide.csv', sep = ',', quote = F,row.names = F, col.names = T)

ht = oncoPrint(df_wide, 
          alter_fun = 
            function(x, y, w, h, v) {
              n = sum(v)  # how many alterations for current gene in current sample
              h = h*0.9
              grid.rect(x, y, w * 0.9, h, gp = gpar(fill = "#CCCCCC", col = NA))
              
              # use `names(which(v))` to correctly map between `v` and `col`
              if(n) {grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.9, 1/n*h,
                               gp = gpar(fill = karyo_col[names(which(v))], 
                                         karyo_col = NA,
                                         col = NA
                                         ), 
                               just = "top")
                
                # grid.points(x,y,h,w, pch = muts[names(which(v))])
                # grid.points(x,y,h,w, pch = muts[names(which(v))])
              } #else {
              # grid.points(x,y,h,w, pch = muts[names(which(v))])
                # grid.rect(x, y, w, h, gp = gpar(fill = "#CCCCCC", col = NA))
                # }
            }, 
          na_col = "#CCCCCC",
          
          #   background = function(x, y, w, h, v) {
          #     grid.rect(x, y, w, h, gp = gpar(fill = "#CCCCCC", col = NA))
          #   }
          # ),
          col = karyo_col, 
          show_column_names = F, 
          show_row_names = T, 
          column_order = as.character(ann_data$Sample), 
          top_annotation = top_ann,
          column_split = factor(ann_data$SPN), 
          show_pct = F, 
          right_annotation = NULL,
          bottom_annotation = bottom_ann, 
          # right_annotation = left_ann,
          row_names_side = 'left', 
          heatmap_legend_param = list(direction = 'horizontal', nrow = 2), 
          name = 'SCOUT')

# creating rigth annotation --> driver type 

drivers_ann_data = all_mut_info %>% 
  # ungroup() %>% 
  dplyr::select(driver) %>% 
  distinct() %>% 
  ungroup() %>% 
  dplyr::mutate(driver = case_when(
    code == "APC R1450*" ~ "Driver SNV/CNA",
    code == "PTEN R130G" ~ "Driver SNV/CNA",
    code == "KRAS G12D"  ~ "Driver SNV/CNA",
    TRUE ~ driver   # <- replaces .default
  ))

genes_order = df_wide[row_order(ht), ] %>% rownames()

drivers_ann_data = drivers_ann_data[which(drivers_ann_data$code %in% genes_order), ] %>% 
  dplyr::select(driver)

left_ann = create_annotation(drivers_ann_data, position = 'row', 
                             ann_colors = list(
                               driver = setNames(c('#640D5F', '#F78D60', '#77BEF0'), 
                                                 nm = drivers_ann_data$driver %>% unique)
                             ), 
                             pos = 'bottom')

pdf('scout_oncoprint_v2.pdf', width = 13, height = 12)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", annotation_legend_list = sign_legends)
dev.off()


