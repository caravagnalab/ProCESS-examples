library(tidyverse)
library(vcfR)

COLS_TO_KEEP <- c("chr", "from", "to", "ref","alt","mutationID","NV","DP","VAF","FILTER")


process_strelka_results = function(gt_path, spn, purity, coverage, chromosome, base_path, outdir) {
  # Extract purity and coverage values from the file path
  # spn <- gsub(".*SCOUT/(SPN[0-9]+).*", "\\1", gt_path)
  # purity <- gsub(".*purity_([0-9.]+).*", "\\1", gt_path)
  # coverage <- gsub(".*coverage_([0-9]+).*", "\\1", gt_path)
  combination = paste0(coverage, "x_", purity, "p")
  
  sample_names = list.files(file.path(outdir, spn, combination, "process"), full.names = F)
  folder_path <- file.path(outdir, spn, combination, "strelka")
  dir.create(folder_path, recursive = T, showWarnings = T)
  
  for (sample in sample_names) {
    message(paste0("Parsing sample ", sample, "..."))
    
    sample_path <- file.path(folder_path, sample)
    dir.create(sample_path, recursive = TRUE, showWarnings = FALSE)
    
    #vcf_folder = file.path(strelka_vcfs_dir, paste0(sample, "_vs_normal_sample"))
    #vcf_folder = paste0("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/",spn,"/sarek/",coverage,"x_",purity,"p/variant_calling/strelka/", sample, "_vs_normal_sample/")
    # vcf_files = list.files(vcf_folder, full.names = T)
    # vcf_files = vcf_files[!grepl(".tbi", vcf_files)]
    
    for (mutation in c("SNV", "INDEL")) {
      message(paste0("Parsing ", mutation, " mutations..."))
      mut_path <- file.path(sample_path, mutation)
      dir.create(mut_path, recursive = TRUE, showWarnings = FALSE)
      
      if (mutation == "SNV") {
        vcf_path = get_sarek_vcf_file(spn, sample, coverage, purity, caller = "strelka", type = "tumour", basedir = base_path)$snvs_vcf
      } else {
        vcf_path = get_sarek_vcf_file(spn, sample, coverage, purity, caller = "strelka", type = "tumour", basedir = base_path)$indels_vcf
      }
      
      vcf = vcfR::read.vcfR(vcf_path)
      
      message(paste0("Parsing ", chromosome, "..."))
      strelka_res = get_strelka_res(vcf, sample_id = sample, filter_mutations = F, chromosome = paste0("chr", chromosome), mut_type = mutation)
      
      # Save the processed mutation data
      file_name <- file.path(mut_path, paste0("chr", chromosome, ".rds"))
      saveRDS(strelka_res, file_name)
      
      # chromosomes = unique(vcfR::getCHROM(vcf))
      # 
      # # Process each chromosome in parallel
      # parallel::mclapply(chromosomes, function(chromosome) {
      #   message(paste0("Parsing ", chromosome, "..."))
      #   strelka_res = get_strelka_res(vcf, sample_id = sample, filter_mutations = F, chromosome = chromosome, mut_type = mutation)
      #   
      #   # Save the processed mutation data
      #   file_name <- file.path(mut_path, paste0(chromosome, ".rds"))
      #   saveRDS(strelka_res, file_name)
      # }, mc.cores = 1)  # Utilize 4 CPU cores
    }
  }
}

retrieve_ref_alt = function(row, mut_type) {
  ref = row$ref
  alt = row$alt
  
  if (mut_type == "SNV") {
    # Process SNVs using AU, CU, GU, TU fields
    if (ref == 'A') {NR = as.integer(strsplit(row$gt_AU, split=',')[[1]][1])}
    if (ref == 'T') {NR = as.integer(strsplit(row$gt_TU, split=',')[[1]][1])}
    if (ref == 'G') {NR = as.integer(strsplit(row$gt_GU, split=',')[[1]][1])}
    if (ref == 'C') {NR = as.integer(strsplit(row$gt_CU, split=',')[[1]][1])}
    
    if (alt == 'A') {NV = as.integer(strsplit(row$gt_AU, split=',')[[1]][1])}
    if (alt == 'T') {NV = as.integer(strsplit(row$gt_TU, split=',')[[1]][1])}
    if (alt == 'G') {NV = as.integer(strsplit(row$gt_GU, split=',')[[1]][1])}
    if (alt == 'C') {NV = as.integer(strsplit(row$gt_CU, split=',')[[1]][1])}
  } 
  else if (mut_type == "INDEL") {
    # Process INDELs using TAR (reference) and TIR (indel) fields
    # TAR = tier counts for reference alleles
    # TIR = tier counts for indel alleles
    NR = as.integer(strsplit(row$gt_TAR, split=',')[[1]][1])
    NV = as.integer(strsplit(row$gt_TIR, split=',')[[1]][1])
  }
  else {
    stop("Invalid mut_type. Must be 'SNV' or 'INDEL'")
  }
  
  ref_alt = paste0(NR, ',', NV)
  ref_alt
}

parse_Strelka = function(vcf, sample_id, filter_mutations = FALSE, chromosome = NULL, mut_type = NULL) {
  # Validate mut_type parameter
  if (!(mut_type %in% c("SNV", "INDEL"))) {
    stop("mut_type must be either 'SNV' or 'INDEL'")
  }
  
  # Filter by chromosome if specified
  if (!is.null(chromosome)) {
    vcf = vcf[vcfR::getCHROM(vcf) == chromosome, ]
  }
  
  # Parse the VCF file into a tidy format
  tb = vcfR::vcfR2tidy(vcf)
  
  # Extract genotype fields and identify samples
  gt_field = tb$gt %>% dplyr::rename(sample = Indiv)
  samples_list = gt_field$sample %>% unique
  
  # Process VEP annotations if present
  if ("CSQ" %in% tb$meta$ID) {
    # VEP specific field extraction
    # Take CSQ field names and split by |
    vep_field = tb$meta %>% 
      dplyr::filter(ID == "CSQ") %>% 
      dplyr::select(Description) %>% 
      dplyr::pull() 
    
    vep_field = strsplit(vep_field, split = "|", fixed = TRUE)[[1]]
    vep_field = vep_field[2:length(vep_field)-1]
    
    # Transform the fix field by splitting the CSQ and select the columns needed
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
      dplyr::select(chr, from, to, ref, alt, CSQ, dplyr::everything(), -ChromKey) %>% 
      tidyr::separate(CSQ, vep_field, sep = "\\|") %>% 
      dplyr::select(chr, from, to, ref, alt, IMPACT, SYMBOL, Gene, dplyr::everything())
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
      dplyr::select(chr, from, to, ref, alt, dplyr::everything(), -ChromKey)
  }
  
  # Filter mutations if requested
  if (filter_mutations) { 
    filter = c('PASS')
  } else {
    filter = fix_field$FILTER %>% unique()
  }
  
  # Determine variant type based on reference and alternate allele length
  if (mut_type == "SNV") {
    # For SNVs: ref and alt should have the same length (1 bp)
    fix_field = fix_field %>% 
      dplyr::filter(nchar(ref) == 1 & nchar(alt) == 1)
  } else if (mut_type == "INDEL") {
    # For INDELs: ref and alt should have different lengths
    fix_field = fix_field %>% 
      dplyr::filter(nchar(ref) != nchar(alt))
  }
  
  # Process each sample
  calls = lapply(
    samples_list,
    function(s) {
      gt_field_s = gt_field %>% dplyr::filter(sample == s)
      
      fits = list()
      fits$sample = s
      
      # Join fix fields and genotype fields
      mutations = dplyr::bind_cols(fix_field, gt_field_s) 
      
      # Skip if no mutations after filtering
      if (nrow(mutations) == 0) {
        fits$mutations = data.frame()
        return(fits)
      }
      
      # Extract ref and alt read counts
      ref_alt = lapply(1:nrow(mutations), function(r) {
        retrieve_ref_alt(mutations[r,], mut_type)
      })
      
      ref_alt = ref_alt %>% unlist()
      mutations$ref_alt = ref_alt
      
      # Process the data
      mutations = mutations %>% 
        tidyr::separate(ref_alt, into = c('NR', 'NV')) %>%
        dplyr::mutate(NR = as.integer(NR), 
                      NV = as.integer(NV)) %>%
        dplyr::mutate(DP = NR+NV) %>% 
        dplyr::mutate(VAF = NV/DP) %>% 
        dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, everything()) %>% 
        dplyr::filter(FILTER %in% filter)
      
      fits$mutations = mutations
      fits
    })
  
  names(calls) = samples_list
  return(calls)
}

get_strelka_res = function(vcf, sample_id, filter_mutations = FALSE, chromosome = NULL, mut_type = NULL) {
  muts <- parse_Strelka(vcf = vcf, filter_mutations = FALSE, chromosome = chromosome, mut_type = mut_type)
  
  muts$TUMOR$mutations <- muts$TUMOR$mutations %>%
    # dplyr::filter(chr == str_replace(chromosome, "chr", "")) %>% 
    dplyr::mutate(mutationID = paste0(chr,":",from,":",ref,":",alt))
  
  strelka_res = muts$TUMOR$mutations[,c(COLS_TO_KEEP)]
  strelka_res
}