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
  rename(end = end.x) %>% mutate(chr =as.integer(chr))


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

parse_Mutect = function(vcf, tumour_id, normal_id){
  # Transform vcf to tidy
  tb = vcfR::vcfR2tidy(vcf)
  
  # Extract gt field and obtain coverage (DP) and variant allele frequency (VAF) fields
  gt_field = tb[["gt"]] %>%
    tidyr::separate(gt_AD, sep = ",", into = c("NR", "NV")) %>%
    dplyr::mutate(
      NR = as.numeric(NR),
      NV = as.numeric(NV),
      DP = NV + NR,
      VAF = NV/DP) %>%
    dplyr::rename(sample = Indiv)
  
  fix_field = tb[["fix"]] %>%
    dplyr::rename(
      chr = CHROM,
      from = POS,
      ref = REF,
      alt = ALT) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      from = as.numeric(from),
      to = from + nchar(alt)) %>%
    dplyr::ungroup() %>%
    dplyr::select(chr, from, to, ref, alt, dplyr::everything(), -ChromKey, -DP)
  
  # Extract sample names
  samples_list = gt_field[["sample"]] %>% unique
  
  calls = lapply(
    samples_list,
    function(s){
      gt_field_s = gt_field %>% dplyr::filter(sample == s)
      
      if(nrow(fix_field) != nrow(gt_field_s))
        stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")
      
      fits = list()
      fits[["sample"]] = s
      fits[["mutations"]] = dplyr::bind_cols(fix_field, gt_field_s) %>%
        dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, dplyr::everything())
      fits
    })
  
  names(calls) = samples_list
  samples = c(tumour_id, normal_id)
  calls = calls[samples]
  
  # check if VCF is annotated with VEP
  if ("CSQ" %in% tb[["meta"]][["ID"]]){
    # VEP specific field extraction
    # Take CSQ field names and split by |
    
    vep_field = tb[['meta']] %>%
      dplyr::filter(ID == "CSQ") %>%
      dplyr::select(Description) %>%
      dplyr::pull()
    
    tmp_vep_field = strsplit(vep_field, split = "|", fixed = TRUE) %>% unlist()
    vep_field = tmp_vep_field[1:length(tmp_vep_field)-1]
    
    # Tranform the fix field by splittig the CSQ and select the columns needed
    calls[[tumour_id]][['mutations']] = calls[[tumour_id]][['mutations']] %>%
      dplyr::mutate(CSQ = strsplit(CSQ, ",")) %>%
      tidyr::unnest(CSQ) %>%
      tidyr::separate(CSQ, vep_field, sep = "\\\\|") %>%
      dplyr::select(chr, from, to, ref, alt, IMPACT, SYMBOL, Gene, dplyr::everything())  #can add other thing, CSQ, HGSP
    
    calls[[normal_id]][['mutations']] = calls[[normal_id]][['mutations']] %>% dplyr::select(-CSQ) %>% dplyr::distinct()
  }
  return(calls)
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
SPNS <- c("SPN01","SPN02","SPN03","SPN04", "SPN06",'SPN07')
all_df <- tibble()

for (tool in c('ascat', 'sequenza')){
  print(tool)
  for (spn in SPNS){
    print(spn)
    samples <- get_sample_names(spn = spn)
    for (s in samples){
      
      if (tool == 'ascat'){
        cna_data <- readRDS(get_tumourevo_qc(spn = spn, 
                                             coverage = 100, 
                                             purity = 0.9,
                                             tool = 'CNAqc',
                                             vcf_caller = 'mutect2', 
                                             cna_caller = tool, 
                                             sample = s)$qc_rds)
        pred_cna_gene <- CNAqc::CNA_gene(cna_data, genes = all_gene) %>% 
          mutate(CN = Major + minor) %>% 
          select(gene, karyotype, CN) %>% 
          filter(karyotype != 'NA:NA') 
      } else if (tool == 'sequenza'){
        sequenza <- get_sarek_cna_file(spn = spn, sampleID = s, coverage = 100, purity = 0.9, caller = 'sequenza')
        cna <- parse_Sequenza(segments_file = sequenza$segments, extra_file = sequenza$confints_CP)
        vcf <- vcfR::read.vcfR(get_sarek_vcf_file(spn = spn, sampleID = s, coverage = 100, purity = 0.9, caller = 'mutect2', type = 'tumour')$vcf)
        muts <- parse_Mutect(vcf = vcf, 
                             tumour_id = paste0(spn, '_', s),
                             normal_id = paste0(spn, '_normal_sample'))
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
        inner_join(pred_cna_gene %>% rename(pred = karyotype, pred_CN = CN), by = "gene") %>% 
        mutate(delta = CN - pred_CN) %>% 
        mutate(spn = spn, sample = s, tool = tool) %>%
        left_join(driver %>% select(code = code, driver), 
                  by = c("gene" = "code")) %>% 
        mutate(is_driver = ifelse(is.na(driver), F, T)) %>% 
        select(-driver) %>% 
        distinct()
      all_df <- bind_rows(all_df, df)
      
    }
  }
}


all <- all_df %>% 
  ggplot() +
  geom_jitter(aes(x = tool, y = delta, shape = is_driver, col = spn), size = 2, height = 0.2) +
  theme_minimal()


library(dplyr)
library(ggrepel)

all_df <- all_df %>% 
  filter(is_driver == T)

set.seed(1)  # reproducible jitter
all_df2 <- all_df %>% 
  mutate(
    x_jit = as.numeric(factor(tool)) + runif(n(), -0.2, 0.2),
    y_jit = delta + runif(n(), -0.2, 0.2)
  )

df_lab <- subset(all_df2, is_driver == TRUE & delta != 0)

ggplot(all_df2) +
  geom_point(
    aes(x = x_jit, y = y_jit, shape = is_driver, col = spn),
    size = 2
  ) +
  geom_label_repel(
    data = subset(all_df2, is_driver == TRUE & delta != 0),
    aes(x = x_jit, y = y_jit, label = gene, col = spn),
    size = 3,
    arrow = arrow(length = unit(0.015, "npc"))
  ) +
  scale_x_continuous(
    breaks = unique(as.numeric(factor(all_df2$tool))),
    labels = unique(all_df2$tool)
  ) +
  theme_minimal() + 
  ylab('Delta CN') + 
  xlab('Tool')


only_driver <- all_df %>% 
  filter(is_driver == T) %>% 
  ggplot() +
  geom_jitter(aes(x = tool, y = delta, shape = is_driver, col = spn), size = 2, height = 0.2) +
  theme_minimal()

# all_df %>% 
#   ggplot() +
#   geom_jitter(aes(x = CN, y = pred_CN, col = is_driver))



# # list of all karyotype classes
# classes <- sort(unique(c(df$true, df$pred)))
# 
# # function to compute confusion stats for a single class
# compute_stats <- function(class) {
#   TP <- sum(df$true == class & df$pred == class)
#   FP <- sum(df$true != class & df$pred == class)
#   FN <- sum(df$true == class & df$pred != class)
#   TN <- sum(df$true != class & df$pred != class)
#   
#   tibble(
#     class = class,
#     TP = TP,
#     FP = FP,
#     FN = FN,
#     TN = TN,
#     sensitivity = ifelse(TP + FN == 0, NA, TP / (TP + FN)),
#     specificity = ifelse(TN + FP == 0, NA, TN / (TN + FP)),
#     precision   = ifelse(TP + FP == 0, NA, TP / (TP + FP)),
#     F1          = ifelse(TP == 0, 0, 2 * TP / (2*TP + FP + FN))
#   )
# }
# 
# # compute for all classes
# metrics <- map_df(classes, compute_stats)
# overall_accuracy <- mean(df$true == df$pred)
# overall_accuracy
