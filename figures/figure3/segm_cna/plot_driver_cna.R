library(tidyverse)
library(ggplot2)
library(ggrepel)
library(gmodels)
set.seed(123)


source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils.R")
setwd("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/segm_cna")
spns_details <- readRDS('spn_details.rds')
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

cna_v1 = read.delim('cna_driver.csv', sep = ',')
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


parse_Sequenza = function(segments_file, extra_file){
  # Extract the segments information
  segments = readr::read_tsv(segments_file, col_types = readr::cols()) %>%
    dplyr::rename(
      chr = chromosome,
      from = start.pos,
      to = end.pos,
      Major = A,
      minor = B) %>%
    dplyr::select(chr, from, to, Major, minor, dplyr::everything())

  solutions = readr::read_tsv(extra_file, col_types = readr::cols())
  purity = solutions[["cellularity"]][2]
  ploidy = solutions[["ploidy.estimate"]][2]
  return(list(segments = segments, purity = purity, ploidy = ploidy))
}

## genes ####
all_gene <-  all_mut_info %>%
   dplyr::mutate(code = str_remove(code, " .*$")) %>%
   pull(code) %>% unique()

gene_pos <- CNAqc::gene_coordinates_GRCh38 %>%
  filter(gene %in% all_gene)


source('../../../getters/process_getters.R')
source('../../../getters/tumourevo_getters.R')
source('../../../getters/sarek_getters.R')
SPNS <- c("SPN02")#","SPN02","SPN03","SPN04", "SPN06",'SPN07'
all_df <- tibble()

for (tool in c('ascat', 'sequenza')){
  print(tool)
  for (p in c(0.3,0.6, 0.9)){
    print(p)
    for (c in c(50, 100, 150)){
      print(c)
      for (spn in SPNS){
        print(spn)
        samples <- get_sample_names(spn = spn)
        for (s in samples){
          if (tool == 'ascat'){
            cna_data <- readRDS(get_tumourevo_qc(spn = spn,
                                                 coverage = c,
                                                 purity = p,
                                                 tool = 'cnaqc',
                                                 vcf_caller = 'mutect2',
                                                 cna_caller = tool,
                                                 sample = s)$qc_rds)

            pred_cna_gene <- CNAqc::CNA_gene(cna_data, genes = all_gene) %>%
              mutate(CN = Major + minor) %>%
              dplyr::select(gene, karyotype, CN) %>%
              filter(karyotype != 'NA:NA')

          } else if (tool == 'sequenza'){
            sequenza <- get_sarek_cna_file(spn = spn, sampleID = s, coverage = 100, purity = 0.9, caller = 'sequenza')
            cna <- parse_Sequenza(segments_file = sequenza$segments, extra_file = sequenza$confints_CP)
            muts <- readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/', spn, '/tumourevo/100x_0.3p_mutect2_ascat/formatter/vcf2cnaqc/SCOUT/', spn, '/', spn, '_', s, '/SCOUT_', spn, '_', spn, '_', s, '_snv.rds'))
            cna_data <- CNAqc::init(mutations = muts[[paste0(spn, '_', s)]]$mutations, cna = cna$segments, purity = cna$purity, ref = 'GRCh38')
            pred_cna_gene <- CNAqc::CNA_gene(cna_data, genes = all_gene) %>%
              mutate(CN = Major + minor) %>%
              select(gene, karyotype, CN) %>%
              filter(karyotype != 'NA:NA')
          }
          true_cna <- readRDS(get_process_cna(spn = spn,
                                              sample = s)) %>% mutate(chr = paste0('chr', chr))
          driver <- all_mut_info %>%
            filter(sample == s) %>%
            dplyr::mutate(code = str_remove(code, " .*$")) %>%
            select(sample, code, driver) %>%
            distinct()

          true_cna_gene <- gene_pos %>%
            left_join(true_cna, by = join_by(chr), relationship = "many-to-many") %>%
            filter(from >= begin, to <= end) %>%
            filter(ratio>.1) %>%
            group_by(gene) %>%
            slice_max(ratio, n = 1, with_ties = FALSE) %>%
            ungroup() %>%
            mutate(karyotype = paste(major, minor, sep = ':')) %>%
            mutate(CN = major + minor) %>%
            select(gene, karyotype, CN, ratio)


          df <- true_cna_gene %>%
            select(gene, true = karyotype, CN) %>%
            inner_join(pred_cna_gene %>% dplyr::rename(pred = karyotype, pred_CN = CN), by = "gene") %>%
            mutate(delta = CN - pred_CN) %>%
            mutate(spn = spn, sample = s, tool = tool) %>%
            left_join(driver %>% select(code = code, driver),
                      by = c("gene" = "code")) %>%
            mutate(is_driver = ifelse(is.na(driver), F, T)) %>%
            select(-driver) %>%
            distinct() %>%
            mutate(purity = p, coverage = c)
         all_df <- bind_rows(all_df, df)
        }
      }
    }
  }
}

all_df = readRDS( '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/segm_cna/driver_df.rds')
all_df_spn05 = readRDS( '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/segm_cna/driver_df_spn05.rds')
all_df = all_df %>% bind_rows(all_df_spn05)

fga <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure3/fga_df.rds') %>% 
  distinct() %>% 
  dplyr::rename(sample = spn)

oncogene <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/ref/oncogenes.rds')
supp <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/ref/suppressor.rds')
info <- readRDS('driver_ann.rds')
colnames(info) <- c('gene', 'class')

all_df = all_df %>%
  mutate(type = 'NA') %>% 
  mutate(type = ifelse(gene %in% oncogene$gene, 'Oncogene', type)) %>% 
  mutate(type = ifelse(gene %in% supp$gene, 'TSG', type)) %>% 
  mutate(type = ifelse(gene == 'ATRX',  'TSG', type)) %>% 
  left_join(info)

# do plot
all_df2 <- all_df %>%
  left_join(fga) %>% 
  filter(is_driver == TRUE)


# Exact string match: true == pred_CN (as strings)
# Integer match: CN == pred_CN (as integers)
accuracy <- all_df2 %>%
  summarise(
    n = n(),
    
    # Exact string match
    string_correct = sum(true == pred_CN, na.rm = TRUE),
    string_accuracy = string_correct / n,
    
    # Integer CN match
    int_correct = sum(CN == pred_CN, na.rm = TRUE),
    int_accuracy = int_correct / n
  )



gene_class_df <- all_df2 %>%
  distinct(gene, class) %>%
  mutate(
    gene = factor(gene, levels = levels(
      all_df2 %>% pull(gene) %>% factor()
    ))
  )

all_df2 <- all_df2 %>% 
  mutate(delta = abs(delta)) %>% 
  group_by(gene, fga_class, type) %>% 
  summarise(
    mean  = ci(delta)[[1]],
    variance = var(delta),
    sd = sd(delta),
    lower = mean - sd, #ci(delta)[[2]],
    upper = mean + sd, #ci(delta)[[3]],
    n=n(),
    .groups = "drop"
  ) %>% 
  mutate(gene = factor(gene, levels = unique(gene))) %>% 
  mutate(fga_class = factor(fga_class, levels = c("High FGA", "Low FGA"))) #%>% 
  #mutate(tool = ifelse(tool == 'ascat', 'ASCAT', 'Sequenza'))

gene_order = all_df2 %>% 
  mutate(gene = as.character(gene)) %>% 
  group_by(gene) %>% 
  reframe(max.m = max(abs(mean)), across(everything())) %>% 
  arrange(desc(max.m)) %>% pull(gene) %>% unique()

all_df2 = all_df2 %>% mutate(gene = factor(gene, levels = gene_order %>% rev()))
gene_class_df = gene_class_df %>% mutate(gene = factor(gene, levels = gene_order %>% rev()))


driver_gene_plot <- ggplot() + 
  geom_tile(data = gene_class_df, aes(x = gene, y = -.5, fill = class), height = .3, col = 'white') +
  geom_hline(aes(yintercept =0), linetype = 2, linewidth = .2, col ='gray30') +
  scale_fill_manual(
    name = "Driver class",
    values = c(
      "SNV" = "lightblue",
      "CNA" = "khaki",
      'CNA,SNV' = "darkseagreen"
    )
  ) +
  #geom_point(data = all_df2, aes(x = gene, y = mean, color = fga_class, group = fga_class, shape = as.factor(type)), position = position_dodge(width = 0.6), size = 2) + 
  geom_pointrange(data = all_df2, 
    aes(x = gene, 
        y = mean, 
        color = fga_class, 
        group = fga_class,
        shape = as.factor(type),
        ymin = lower, ymax = upper),
    position = position_dodge(width = 0.6), 
    size = .5
  ) +
  scale_color_manual('FGA class',
                     values = c("High FGA" = "indianred3", "Low FGA" = "dodgerblue3")) +
  scale_shape_manual('Gene type', values = c(16, 17)) + 
  my_ggplot_theme() +
  #facet_wrap(.~tool, ncol = 1) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x  = element_blank()
  ) + 
  ylab('MAE total CN') + xlab('Gene')  +
  theme(
    #axis.text.x = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0.00001, "cm")
  )

driver_gene_plot
#ggsave(filename = 'driver.pdf', dpi = 300, width = 6, height = 3)
