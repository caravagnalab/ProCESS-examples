library(tidyverse)
library(vcfR)

COLS_TO_KEEP <- c("chr", "from", "to", "ref","alt","mutationID","NV","DP","VAF","FILTER")
# parse_FreeBayes = function(vcf, sample_id, filter_mutations = FALSE, chromosome = NULL, mut_type = NULL){
#   # Validate mut_type parameter
#   if (!is.null(mut_type) && !(mut_type %in% c("SNV", "INDEL"))) {
#     stop("mut_type must be either 'SNV' or 'INDEL'")
#   }
#   
#   # Filter by chromosome if specified
#   if (!is.null(chromosome)) {
#     vcf = vcf[vcfR::getCHROM(vcf) == chromosome, ]
#   }
#   
#   # Parse the VCF file into a tidy format
#   tb = vcfR::vcfR2tidy(vcf)
#   
#   gt_field = tb$gt %>% 
#     tidyr::separate(gt_AD, sep = ",", into = c("NR", "NV")) %>%
#     dplyr::mutate(
#       NR = as.numeric(NR),
#       NV = as.numeric(NV),
#       DP = NV + NR,
#       VAF = NV/DP) %>% 
#     dplyr::rename(sample = Indiv)
#   
#   samples_list = gt_field$sample %>% unique
#   
#   if ("CSQ" %in% tb$meta$ID){
#     # VEP specific field extraction
#     # Take CSQ field names and split by |
#     
#     vep_field = tb$meta %>%
#       dplyr::filter(ID == "CSQ") %>%
#       dplyr::select(Description) %>%
#       dplyr::pull()
#     
#     vep_field = strsplit(vep_field, split = "|", fixed = TRUE)[[1]]
#     vep_field = vep_field[2:length(vep_field)-1]
#     
#     # Transform the fix field by splitting the CSQ and select the columns needed
#     fix_field = tb$fix %>%
#       dplyr::rename(
#         chr = CHROM,
#         from = POS,
#         ref = REF,
#         alt = ALT,
#         qual = QUAL  # Explicitly rename the QUAL field
#       ) %>%
#       dplyr::rowwise() %>%
#       dplyr::mutate(
#         from = as.numeric(from),
#         to = from + nchar(alt),
#         qual = as.numeric(qual)  # Ensure quality is numeric
#       ) %>%
#       dplyr::ungroup() %>%
#       dplyr::select(chr, from, to, ref, alt, qual, CSQ, dplyr::everything(),  -ChromKey, -DP) %>%  # Include qual in selection
#       tidyr::separate(CSQ, vep_field, sep = "\\|") %>%
#       dplyr::select(chr, from, to, ref, alt, qual, IMPACT, SYMBOL, Gene, dplyr::everything())  # Include qual in final selection
#     
#   } else {
#     fix_field = tb$fix %>%
#       dplyr::rename(
#         chr = CHROM,
#         from = POS,
#         ref = REF,
#         alt = ALT,
#         qual = QUAL  # Explicitly rename the QUAL field
#       ) %>%
#       dplyr::rowwise() %>%
#       dplyr::mutate(
#         from = as.numeric(from),
#         to = from + nchar(alt),
#         qual = as.numeric(qual)  # Ensure quality is numeric
#       ) %>%
#       dplyr::ungroup() %>%
#       dplyr::select(chr, from, to, ref, alt, qual, dplyr::everything(), -ChromKey, -DP)  # Include qual in selection
#   }
#   
#   # if have to filter mutations
#   if (filter_mutations){
#     filter = c('PASS')
#   } else {
#     filter = fix_field$FILTER %>% unique()
#   }
#   
#   calls = lapply(
#     samples_list,
#     function(s){
#       gt_field_s = gt_field %>% dplyr::filter(sample == s)
#       
#       if(nrow(fix_field) != nrow(gt_field_s))
#         stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")
#       
#       fits = list()
#       fits$sample = s
#       fits$mutations = dplyr::bind_cols(fix_field, gt_field_s) %>%
#         dplyr::select(chr, from, to, ref, alt, qual, NV, DP, VAF, dplyr::everything()) %>%  # Include qual in final mutation data
#         dplyr::filter(FILTER %in% filter)
#       
#       # Apply mutation type filtering if specified
#       if (!is.null(mut_type)) {
#         if (mut_type == "SNV") {
#           fits$mutations = fits$mutations %>% 
#             dplyr::filter(nchar(ref) == 1 & nchar(alt) == 1)
#         } else if (mut_type == "INDEL") {
#           fits$mutations = fits$mutations %>% 
#             dplyr::filter(nchar(ref) != 1 | nchar(alt) != 1)
#         }
#       }
#       
#       fits
#     })
#   
#   names(calls) = samples_list
#   
#   # Filter by sample_id if provided
#   if (!is.null(sample_id) && sample_id %in% names(calls)) {
#     calls = calls[sample_id]
#   }
# 
#   return(calls)
# }

find_inserted_substring <- function(s1, s2) {
  
  if (nchar(s2) <= nchar(s1)) {
    return(NULL)
  }
  
  # Try each possible insertion position in s1
  for (i in 0:nchar(s1)) {
    # Split s1 at position i
    if (i == 0) {
      prefix <- ""
      suffix <- s1
    } else if (i == nchar(s1)) {
      prefix <- s1
      suffix <- ""
    } else {
      prefix <- substr(s1, 1, i)
      suffix <- substr(s1, i + 1, nchar(s1))
    }
    
    # Check if s2 starts with prefix and ends with suffix
    starts_with_prefix <- (prefix == "" || substr(s2, 1, nchar(prefix)) == prefix)
    ends_with_suffix <- (suffix == "" || substr(s2, nchar(s2) - nchar(suffix) + 1, nchar(s2)) == suffix)
    
    if (starts_with_prefix && ends_with_suffix) {
      # Calculate the length of the inserted part
      inserted_length <- nchar(s2) - nchar(prefix) - nchar(suffix)
      
      # Extract the inserted substring
      if (inserted_length > 0) {
        inserted <- substr(s2, nchar(prefix) + 1, nchar(prefix) + inserted_length)
        return(inserted)
      }
    }
  }
  
  return(NULL)
}

parse_FreeBayes = function(vcf, filter_mutations = FALSE, chromosome = NULL, mut_type = NULL, min_vaf = 0.01, max_normal_vaf = 0.02){
  # Validate mut_type parameter
  if (!is.null(mut_type) && !(mut_type %in% c("SNV", "INDEL"))) {
    stop("mut_type must be either 'SNV' or 'INDEL'")
  }
  
  # Filter by chromosome if specified
  if (!is.null(chromosome)) {
    vcf = vcf[vcfR::getCHROM(vcf) == chromosome, ]
  }
  
  # Parse the VCF file into a tidy format
  tb = vcfR::vcfR2tidy(vcf)
  
  gt_field = tb$gt %>% 
    tidyr::separate(gt_AD, sep = ",", into = c("NR", "NV")) %>%
    dplyr::mutate(
      NR = as.numeric(NR),
      NV = as.numeric(NV),
      DP = NV + NR,
      VAF = NV/DP) %>% 
    dplyr::rename(sample = Indiv)
  
  samples_list = gt_field$sample %>% unique
  
  # Get tumour and normal id
  tumor_id = samples_list[!grepl("normal", samples_list)]
  normal_id = samples_list[grepl("normal", samples_list)]
  
  if ("CSQ" %in% tb$meta$ID){
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
        alt = ALT,
        qual = QUAL  # Explicitly rename the QUAL field
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        from = as.numeric(from),
        to = from + nchar(alt),
        qual = as.numeric(qual)  # Ensure quality is numeric
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(chr, from, to, ref, alt, qual, CSQ, dplyr::everything(),  -ChromKey, -DP) %>%  # Include qual in selection
      tidyr::separate(CSQ, vep_field, sep = "\\|") %>%
      dplyr::select(chr, from, to, ref, alt, qual, IMPACT, SYMBOL, Gene, dplyr::everything())  # Include qual in final selection
    
  } else {
    fix_field = tb$fix %>%
      dplyr::rename(
        chr = CHROM,
        from = POS,
        ref = REF,
        alt = ALT,
        qual = QUAL  # Explicitly rename the QUAL field
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        from = as.numeric(from),
        to = from + nchar(alt),
        qual = as.numeric(qual)  # Ensure quality is numeric
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(chr, from, to, ref, alt, qual, dplyr::everything(), -ChromKey, -DP)  # Include qual in selection
  }
  
  # if have to filter mutations
  if (filter_mutations){
    filter = c('PASS')
  } else {
    filter = fix_field$FILTER %>% unique()
  }
  
  # Extract genotype data for tumor and normal samples
  gt_field_tumor = gt_field %>% dplyr::filter(sample == tumor_id)
  gt_field_normal = gt_field %>% dplyr::filter(sample == normal_id)
  
  # Check for mismatches
  if (nrow(fix_field) != nrow(gt_field_tumor) || nrow(fix_field) != nrow(gt_field_normal)) {
    stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")
  }
  
  # Filter for somatic mutations by using both tumor and normal data
  # First create a temporary combined dataset to apply filtering
  temp_combined = dplyr::bind_cols(
    fix_field,
    gt_field_tumor %>% dplyr::select(NR, NV, DP, VAF),
    gt_field_normal %>% dplyr::select(VAF) %>% dplyr::rename(N_VAF = VAF)
  )
  
  # Get the indices of somatic mutations
  somatic_indices = which(
    temp_combined$FILTER %in% filter &
      temp_combined$VAF >= min_vaf &
      temp_combined$N_VAF <= max_normal_vaf &
      temp_combined$DP >= 10
  )
  
  # Now create the final tumor-only mutation data
  somatic_mutations = dplyr::bind_cols(fix_field, gt_field_tumor) %>%
    dplyr::slice(somatic_indices) %>%
    dplyr::select(chr, from, to, ref, alt, qual, NV, DP, VAF, dplyr::everything())
  
  # Apply mutation type filtering if specified
  if (!is.null(mut_type)) {
    if (mut_type == "SNV") {
      somatic_mutations = somatic_mutations %>% 
        dplyr::filter(nchar(ref) == 1 & nchar(alt) == 1)
    } else if (mut_type == "INDEL") {
      somatic_mutations = somatic_mutations %>% 
        dplyr::filter(nchar(ref) != 1 | nchar(alt) != 1)
      
      # Parse correctly
      parsed_ref_alt_df = lapply(1:nrow(somatic_mutations), function(i) {
        ref = somatic_mutations[i,]$ref
        alt = somatic_mutations[i,]$alt
        
        alts = unlist(strsplit(alt, ","))
        if (length(alts) != 1) {
          df_alts = lapply(alts, function(alt) {
            insertion = find_inserted_substring(ref, alt)
            if (!is.null(insertion)) {
              new_ref = str_split(ref, insertion)[[1]][1]
              new_alt = paste0(new_ref, insertion)
              dplyr::tibble(new_ref, new_alt)
            }    
          }) %>% do.call("bind_rows", .)
          if (nrow(df_alts) != 0) {
            df_alts[1,]
          } else {
            dplyr::tibble(new_ref = ref, new_alt = alt)
          }
        } else {
          insertion = find_inserted_substring(ref, alt)
          if (is.null(insertion)) {
            dplyr::tibble(new_ref = ref, new_alt = alt)
          } else {
            new_ref = str_split(ref, insertion)[[1]][1]
            new_alt = paste0(new_ref, insertion)
            dplyr::tibble(new_ref, new_alt)
          }  
        } 
      }) %>% do.call(bind_rows, .)
      
      somatic_mutations = cbind(somatic_mutations, parsed_ref_alt_df)
      somatic_mutations = somatic_mutations %>% 
        dplyr::mutate(alt = new_alt, ref=new_ref)
    }
  }
  
  # Create the result structure in the original format but with somatic mutations only
  calls = list()
  calls[[tumor_id]] = list(
    sample = tumor_id,
    mutations = somatic_mutations
  )
  
  return(calls)
}

get_freeBayes_res = function(vcf, sample_id, filter_mutations = FALSE, chromosome = NULL, mut_type = NULL, pass_quality = 20, min_vaf = 0.01, max_normal_vaf = 0.02) {
  muts <- parse_FreeBayes(vcf = vcf, filter_mutations = FALSE, chromosome = chromosome, mut_type = mut_type, min_vaf = min_vaf, max_normal_vaf = max_normal_vaf)
  
  tum_name = names(muts)[!grepl("normal", names(muts))]
  muts[[tum_name]]$mutations <- muts[[tum_name]]$mutations %>%
    # dplyr::filter(chr == str_replace(chromosome, "chr", "")) %>% 
    dplyr::mutate(mutationID = paste0(chr,":",from,":",ref,":",alt)) %>% 
    dplyr::mutate(FILTER = if_else(qual >= pass_quality, "PASS", "LOW QUALITY"))
  
  muts[[tum_name]]$mutations[,c(COLS_TO_KEEP)]
}


process_freebayes_results = function(gt_path, 
                                     spn, 
                                     purity, 
                                     coverage, 
                                     chromosome, 
                                     base_path, 
                                     outdir, 
                                     pass_quality, 
                                     min_vaf, 
                                     max_normal_vaf) {
  
  # Extract purity and coverage values from the file path
  combination = paste0(coverage, "x_", purity, "p")
  
  sample_names = list.files(file.path(outdir, spn, combination, "process"), full.names = F)
  folder_path <- file.path(outdir, spn, combination, "freebayes")
  dir.create(folder_path, recursive = T, showWarnings = T)
  
  for (sample in sample_names) {
    message(paste0("Parsing sample ", sample, "..."))
    
    sample_path <- file.path(folder_path, sample)
    dir.create(sample_path, recursive = TRUE, showWarnings = FALSE)
    
    #vcf_folder = file.path(freebayes_vcfs_dir, paste0(sample, "_vs_normal_sample"))
    #vcf_folder = paste0("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/",spn,"/sarek/",coverage,"x_",purity,"p/variant_calling/freebayes/", ")
    # vcf_files = list.files(vcf_folder, full.names = T)
    # vcf_path = vcf_files[!grepl(".tbi", vcf_files)]
    
    vcf_path = get_sarek_vcf_file(spn, sample, coverage, purity, caller = "freebayes", type = "tumour", basedir = base_path)$vcf
    vcf = vcfR::read.vcfR(vcf_path)
    
    message(paste0("Parsing ", chromosome, "..."))
    
    mutation = "INDEL"
    for (mutation in c("SNV", "INDEL")) {
      message(paste0("Parsing ", mutation, " mutations..."))
      mut_path <- file.path(sample_path, mutation)
      dir.create(mut_path, recursive = TRUE, showWarnings = FALSE)
      
      mut_data = get_freeBayes_res(
        vcf, 
        sample, 
        filter_mutations = FALSE, 
        chromosome = paste0("chr", chromosome), 
        mut_type = mutation, 
        pass_quality = pass_quality, 
        min_vaf = min_vaf, 
        max_normal_vaf = max_normal_vaf
      )
      
      # if (mutation == "SNV") {
      #   mut_data = caller_res %>% 
      #     dplyr::filter(nchar(ref) == 1 & nchar(alt) == 1)
      # } else if (mutation == "INDEL") {
      #   mut_data = caller_res %>% 
      #     dplyr::filter(nchar(ref) != 1 | nchar(alt) != 1)
      # }
      
      # Save the processed mutation data
      file_name <- file.path(mut_path, paste0("chr", chromosome, ".rds"))
      saveRDS(mut_data, file_name)  
    }
    
    # chromosomes = unique(vcfR::getCHROM(vcf))
    # 
    # # Process each chromosome in parallel
    # parallel::mclapply(chromosomes, function(chromosome) {
    #   message(paste0("Parsing ", chromosome, "..."))
    #   caller_res = get_freeBayes_res(
    #     vcf, 
    #     sample, 
    #     filter_mutations = FALSE, 
    #     chromosome = chromosome, 
    #     mut_type = NULL, 
    #     pass_quality = pass_quality, 
    #     min_vaf = min_vaf, 
    #     max_normal_vaf = max_normal_vaf
    #   )
    #   
    #   for (mutation in c("SNV", "INDEL")) {
    #     message(paste0("Parsing ", mutation, " mutations..."))
    #     mut_path <- file.path(sample_path, mutation)
    #     dir.create(mut_path, recursive = TRUE, showWarnings = FALSE)
    #     
    #     if (mutation == "SNV") {
    #       mut_data = caller_res %>% 
    #         dplyr::filter(nchar(ref) == 1 & nchar(alt) == 1)
    #     } else if (mutation == "INDEL") {
    #       mut_data = caller_res %>% 
    #         dplyr::filter(nchar(ref) != 1 | nchar(alt) != 1)
    #     }
    #     
    #     # Save the processed mutation data
    #     file_name <- file.path(mut_path, paste0(chromosome, ".rds"))
    #     saveRDS(mut_data, file_name)  
    #   }
    # }, mc.cores = 1)  # Utilize 4 CPU cores
  }
}
