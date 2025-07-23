parse_HaplotypeCaller = function(file, out_file = '', save = TRUE){
  vcf = vcfR::read.vcfR(file)
  tb = vcfR::vcfR2tidy(vcf)
  
  gt_field = tb$gt %>% 
    tidyr::separate(gt_AD, sep = ",", into = c("NR", "NV")) %>%
    dplyr::mutate(
      NR = as.numeric(NR),
      NV = as.numeric(NV),
      DP = NV + NR,
      BAF = NV/DP) %>%  
    dplyr::rename(sample = Indiv)
  
  samples_list = gt_field$sample %>% unique
  
  fix_field = tb$fix %>%
    dplyr::rename(
      chr = CHROM,
      from = POS,
      ref = REF,
      alt = ALT
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      from = as.numeric(from),
      to = from + nchar(alt)) %>%
    dplyr::ungroup() %>%
    dplyr::select(chr, from, to, ref, alt, dplyr::everything(), -ChromKey, -DP)
  
  calls = lapply(
    samples_list,
    function(s){
      gt_field_s = gt_field %>% dplyr::filter(sample == s)
      
      if(nrow(fix_field) != nrow(gt_field_s))
        stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")
      
      fits = list()
      fits$sample = s
      fits$mutations = dplyr::bind_cols(fix_field, gt_field_s) %>%
        dplyr::select(chr, from, to, ref, alt, NV, DP, BAF, dplyr::everything())
      fits
    }) 
  
  names(calls) = samples_list
  
  if (save){
    saveRDS(object = calls, file = out_file)
  }
  
  return(calls)
}


parse_freebayes = function(file, out_file = '', save = TRUE, cutoff = 20){
  vcf = vcfR::read.vcfR(file)
  tb = vcfR::vcfR2tidy(vcf)
  
  gt_field = tb$gt %>% 
    tidyr::separate(gt_AD, sep = ",", into = c("NR", "NV")) %>%
    dplyr::mutate(
      NR = as.numeric(NR),
      NV = as.numeric(NV),
      DP = NV + NR,
      BAF = NV/DP) %>%  
    dplyr::rename(sample = Indiv)
  
  samples_list = gt_field$sample %>% unique
  
  fix_field = tb$fix %>%
    dplyr::rename(
      chr = CHROM,
      from = POS,
      ref = REF,
      alt = ALT
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      from = as.numeric(from),
      to = from + nchar(alt)) %>%
    dplyr::ungroup() %>%
    dplyr::select(chr, from, to, ref, alt, dplyr::everything(), -ChromKey, -DP) %>% 
    mutate(FILTER = ifelse(QUAL > cutoff, 'PASS', 'LOW QUALITY'))
    
  
  calls = lapply(
    samples_list,
    function(s){
      gt_field_s = gt_field %>% dplyr::filter(sample == s)
      
      if(nrow(fix_field) != nrow(gt_field_s))
        stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")
      
      fits = list()
      fits$sample = s
      fits$mutations = dplyr::bind_cols(fix_field, gt_field_s) %>%
        dplyr::select(chr, from, to, ref, alt, NV, DP, BAF, dplyr::everything())
      fits
    }) 
  
  names(calls) = samples_list
  
  if (save){
    saveRDS(object = calls, file = out_file)
  }
  
  return(calls)
}

parse_strelka = function(file, out_file = '', save = TRUE){
  vcf = vcfR::read.vcfR(file)
  tb = vcfR::vcfR2tidy(vcf)
  
  gt_field = tb$gt %>% 
    tidyr::separate(gt_AD, sep = ",", into = c("NR", "NV")) %>%
    dplyr::mutate(
      NR = as.numeric(NR),
      NV = as.numeric(NV),
      DP = NV + NR,
      BAF = NV/DP) %>%  
    dplyr::rename(sample = Indiv)
  
  samples_list = gt_field$sample %>% unique
  
  fix_field = tb$fix %>%
    dplyr::rename(
      chr = CHROM,
      from = POS,
      ref = REF,
      alt = ALT
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      from = as.numeric(from),
      to = from + nchar(alt)) %>%
    dplyr::ungroup() %>%
    dplyr::select(chr, from, to, ref, alt, dplyr::everything(), -ChromKey)
  
  
  calls = lapply(
    samples_list,
    function(s){
      gt_field_s = gt_field %>% dplyr::filter(sample == s)
      
      if(nrow(fix_field) != nrow(gt_field_s))
        stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")
      
      fits = list()
      fits$sample = s
      fits$mutations = dplyr::bind_cols(fix_field, gt_field_s) %>%
        dplyr::select(chr, from, to, ref, alt, NV, DP, BAF, dplyr::everything())
      fits
    }) 
  
  names(calls) = samples_list
  
  if (save){
    saveRDS(object = calls, file = out_file)
  }
  
  return(calls)
}

