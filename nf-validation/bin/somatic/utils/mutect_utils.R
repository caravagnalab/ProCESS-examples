library(tidyverse)
library(vcfR)
COLS_TO_KEEP <- c("chr", "from", "to", "ref","alt","mutationID","NV","DP","VAF","FILTER")

parse_mutect2 = function(vcf, sample_id, filter_mutations = FALSE, chromosome = NULL) {
  # Filter by chromosome if specified
  if (!is.null(chromosome)) {
    vcf = vcf[vcfR::getCHROM(vcf) == chromosome, ]
  }
  
  # Transform vcf to tidy 
  tb = vcfR::vcfR2tidy(vcf)
  
  # Extract gt field and obtain coverage (DP) and variant allele frequency (VAF) fields
  gt_field = tb$gt %>% 
    tidyr::separate(gt_AD, sep = ",", into = c("NR", "NV")) %>%
    dplyr::mutate(
      NR = as.numeric(NR),
      NV = as.numeric(NV),
      DP = NV + NR,
      VAF = NV/DP) %>%  
    dplyr::rename(sample = Indiv)
  
  # Extract sample names
  samples_list = gt_field$sample %>% unique
  
  # check if VCF is annotated with VEP
  if ("CSQ" %in% tb$meta$ID){
    # VEP specific field extraction
    # Take CSQ field names and split by |
    
    vep_field = tb$meta %>% 
      dplyr::filter(ID == "CSQ") %>% 
      dplyr::select(Description) %>% 
      dplyr::pull() 
    
    vep_field = strsplit(vep_field, split = "|", fixed = TRUE)[[1]]
    vep_field = vep_field[2:length(vep_field)-1]
    
    # Tranform the fix field by splittig the CSQ and select the columns needed
    fix_field = tb$fix %>%
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
      dplyr::select(chr, from, to, ref, alt, CSQ, dplyr::everything()) %>% 
      tidyr::separate(CSQ, vep_field, sep = "\\|") %>% 
      dplyr::select(chr, from, to, ref, alt, IMPACT, SYMBOL, Gene, dplyr::everything()) #can add other thing, CSQ, HGSP
    
  } else {
    # Take from fix field some columns
    fix_field = tb$fix %>%
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
      dplyr::select(chr, from, to, ref, alt, dplyr::everything(), -ChromKey, -DP) #-DP
  }
  
  # For each sample create the table of mutations 
  calls = lapply(
    samples_list,
    function(s){
      gt_field_s = gt_field %>% dplyr::filter(sample == s)
      
      if(nrow(fix_field) != nrow(gt_field_s))
        stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")
      
      fits = list()
      fits$sample = s
      fits$mutations = dplyr::bind_cols(fix_field, gt_field_s) %>%
        dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, dplyr::everything())
      # dplyr::filter(FILTER %in% filter)
      fits
    })
  
  names(calls) = samples_list
  # normal = names(calls)[length(calls)] 
  # calls = calls[c(sample_id, normal)]
  
  sample_id_vcf = names(calls)[grepl(sample_id, names(calls))]
  
  calls[[sample_id_vcf]]$mutations <- calls[[sample_id_vcf]]$mutations %>%
    #dplyr::filter(chr%in%chromsomes) %>%
    dplyr::mutate(mutationID=paste0(chr,":",from,":",ref,":",alt))
  
  calls[[sample_id_vcf]]$mutations[,c(COLS_TO_KEEP)]
}

process_mutect2_results = function(spn, purity, coverage, chromosome, outdir, vcf_path, sample_id) {

  combination = paste0(coverage, "x_", purity, "p")
  print(combination)


  message(paste0("Parsing sample ", sample_id, "..."))
  
  vcf = vcfR::read.vcfR(vcf_path)
  caller_res = parse_mutect2(vcf = vcf, sample_id = sample_id, chromosome = paste0("chr", chromosome))

  mut_data_snv = caller_res %>% dplyr::filter(nchar(ref) == 1 & nchar(alt) == 1)
  mut_data_indel = caller_res %>% dplyr::filter(nchar(ref) != 1 | nchar(alt) != 1)
  mut_data = list('SNV' = mut_data_snv, 'INDEL' = mut_data_indel)

  file_name <- file.path(paste0(paste0("chr", chromosome), ".rds"))
  saveRDS(mut_data, file_name)
}

#process_mutect2_results = function(spn, purity, coverage, chromosome, base_path, outdir,vcf_path) {
#  # Extract purity and coverage values from the file path
#  # spn <- gsub(".*SCOUT/(SPN[0-9]+).*", "\\1", gt_path)
#  # purity <- gsub(".*purity_([0-9.]+).*", "\\1", gt_path)
#  # coverage <- gsub(".*coverage_([0-9]+).*", "\\1", gt_path)
#  combination = paste0(coverage, "x_", purity, "p")
#  print(combination)
#  sample_names = get_sample_names(spn=spn,base_path=base_path) #list.files(file.path(outdir, spn, combination, "process"), full.names = F)
#  print(sample_names) 
#  folder_path <- file.path(outdir, "somatic", combination, "mutect2")
#  dir.create(folder_path, recursive = T, showWarnings = T)
#  
#  for (sample in sample_names) {
#    message(paste0("Parsing sample ", sample, "..."))
#    
#    sample_path <- file.path(folder_path, sample)
#    dir.create(sample_path, recursive = TRUE, showWarnings = FALSE)
#    
#    #vcf_path = get_sarek_vcf_file(spn, sample, coverage, purity, caller = "mutect2", type = "tumour", basedir = base_path)$vcf
#    
#    #vcf_path = file.path(mutect_vcfs_dir, paste0(spn, ".mutect2.filtered.vcf.gz"))
#    vcf = vcfR::read.vcfR(vcf_path)
#    
#    message(paste0("Parsing ", chromosome, "..."))
#    caller_res = parse_mutect2(vcf = vcf, sample_id = sample, chromosome = paste0("chr", chromosome))
#    
#    for (mutation in c("SNV", "INDEL")) {
#      message(paste0("Parsing ", mutation, " mutations..."))
#      mut_path <- file.path(sample_path, mutation)
#      dir.create(mut_path, recursive = TRUE, showWarnings = FALSE)
#      
#      if (mutation == "SNV") {
#        mut_data = caller_res %>% 
#          dplyr::filter(nchar(ref) == 1 & nchar(alt) == 1)
#      } else if (mutation == "INDEL") {
#        mut_data = caller_res %>% 
#          dplyr::filter(nchar(ref) != 1 | nchar(alt) != 1)
#      }
#      
#      # Save the processed mutation data
#      file_name <- file.path(mut_path, paste0(paste0("chr", chromosome), ".rds"))
#      saveRDS(mut_data, file_name)
#    }
#    
#    # chromosomes = unique(vcfR::getCHROM(vcf))
#    # for (chromosome in chromosomes) {
#    #   message(paste0("Parsing ", chromosome, "..."))
#    #   caller_res = parse_mutect2(vcf, sample_id = sample, chromosome = chromosome)
#    #   
#    #   for (mutation in c("SNV", "INDEL")) {
#    #     message(paste0("Parsing ", mutation, " mutations..."))
#    #     mut_path <- file.path(sample_path, mutation)
#    #     dir.create(mut_path, recursive = TRUE, showWarnings = FALSE)
#    #     
#    #     if (mutation == "SNV") {
#    #       mut_data = caller_res %>% 
#    #         dplyr::filter(nchar(ref) == 1 & nchar(alt) == 1)
#    #     } else if (mutation == "INDEL") {
#    #       mut_data = caller_res %>% 
#    #         dplyr::filter(nchar(ref) != 1 | nchar(alt) != 1)
#    #     }
#    #     
#    #     # Save the processed mutation data
#    #     file_name <- file.path(mut_path, paste0(chromosome, ".rds"))
#    #     saveRDS(mut_data, file_name)
#    #   }
#    # }
#  }
#}
