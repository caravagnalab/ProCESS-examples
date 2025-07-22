library(IRanges)
library(GenomicRanges)

color_by_state = c("INFERRED_Major1"='steelblue',
                   "INFERRED_minor1"='indianred',
                   "INFERRED_CN"='seagreen')

fill_by_match = c('match'= alpha('gainsboro', .03),
                  'no match subclone'= alpha('cadetblue3', .08),
                  'match subclone'= alpha('lightsalmon', .08),
                  'no match' = alpha('indianred', .08),
                  'subclonal' = alpha('lightsalmon', .08),
                  'clonal' = 'gainsboro')

color_type <- c('clonal' = 'gainsboro', 'sub-clonal' = 'lightsalmon')
color_process <- c('Major1'='steelblue', 'minor1'='indianred', 'Major2'='steelblue1', 'minor2'="orangered")
color_caller <- c('ascat' = 'palegreen4','cnvkit' = 'orange', 'sequenza' = 'steelblue')


compute_CNA_metrics <- function(confusion_df) {
  metrics_df <- confusion_df %>%
    mutate(
      precision = TP / (TP + FP),
      recall = TP / (TP + FN),
      f1 = 2 * precision * recall / (precision + recall)
    )
  
  # Handle NaNs
  metrics_df <- metrics_df %>%
    mutate(
      precision = ifelse(is.nan(precision), 0, precision),
      recall = ifelse(is.nan(recall), 0, recall),
      f1 = ifelse(is.nan(f1), 0, f1)
    )
  
  # Genome-wide metrics
  genome_total <- confusion_df %>%
    summarise(
      TP = sum(TP),
      FP = sum(FP),
      FN = sum(FN),
      chr = "genome"
    ) %>%
    mutate(
      precision = TP / (TP + FP),
      recall = TP / (TP + FN),
      f1 = 2 * precision * recall / (precision + recall)
    ) %>%
    mutate(
      precision = ifelse(is.nan(precision), 0, precision),
      recall = ifelse(is.nan(recall), 0, recall),
      f1 = ifelse(is.nan(f1), 0, f1)
    )
  
  bind_rows(metrics_df, genome_total)
}


# Select a chromosome
compute_confusion_matrix <- function(predicted, truth, delta = 10000, delta_merge = NULL) {
  # Step 0: Sort inputs
  predicted <- sort(predicted)
  truth <- sort(truth)
  
  # Step 1: Recursively merge predicted breakpoints
  if (!is.null(delta_merge) && length(predicted) > 1) {
    changed <- TRUE
    while (changed) {
      diffs <- diff(predicted)
      close_pairs <- which(diffs <= delta_merge)
      if (length(close_pairs) == 0) {
        changed <- FALSE
      } else {
        # Merge first close pair (can be done recursively or all at once, we'll go greedy)
        i <- close_pairs[1]
        new_val <- mean(predicted[i:(i+1)])
        predicted <- c(predicted[-c(i, i+1)], new_val)
        predicted <- sort(predicted)
      }
    }
  }
  
  # Step 2: Matching logic
  matched_truth <- logical(length(truth))
  matched_pred <- logical(length(predicted))
  
  for (i in seq_along(predicted)) {
    diffs <- abs(truth - predicted[i])
    # Find unmatched truth within delta
    candidates <- which(diffs <= delta & !matched_truth)
    if (length(candidates) > 0) {
      best_match <- candidates[which.min(diffs[candidates])]
      matched_truth[best_match] <- TRUE
      matched_pred[i] <- TRUE
    }
  }
  
  TP <- sum(matched_pred)
  FP <- sum(!matched_pred)
  FN <- sum(!matched_truth)
  
  dplyr::tibble(TP = TP, FP = FP, FN = FN)
}

absolute_to_relative_coordinates <- function(muts, reference = CNAqc::chr_coordinates_GRCh38, centromere = F){
  vfrom = reference$from
  names(vfrom) = reference$chr
  if (!centromere){
    muts %>%
      mutate(
        start = start + vfrom[chr],
        end = end + vfrom[chr])
  } else if (centromere){
    muts %>%
      mutate(
        start = start + vfrom[chr],
        end = end + vfrom[chr],
        centromere = centromere + vfrom[chr])
  }
}

relative_to_absolute_coordinates <- function(muts, reference = CNAqc::chr_coordinates_GRCh38, centromere = F){
  vfrom = reference$from
  names(vfrom) = reference$chr
  if (!centromere){
    muts %>%
      mutate(
        start = start - vfrom[chr],
        end = end - vfrom[chr])
  } else if (centromere){
    muts %>%
      mutate(
        start = start - vfrom[chr],
        end = end - vfrom[chr],
        centromere = centromere - vfrom[chr])
  }
}


absolute_to_relative_coordinates_muts <- function(muts, reference = CNAqc::chr_coordinates_GRCh38){
  vfrom = reference$from
  names(vfrom) = reference$chr
  muts <- muts %>% mutate(from = from + vfrom[chr])
  return(muts) 
}

## Read data
read_ProCESS = function(spn_id,sample_id,coverage,purity){
  
  # Read CNA data
  CNA_races = readRDS(get_process_cna(spn = spn_id,sample = sample_id)) %>%
    mutate(chr = paste0('chr',chr)) %>% dplyr::rename(from=begin,to=end) %>%
    mutate(seg_id = paste0(chr,':',from,':',to)) %>% as_tibble()
  
  return(list('CNA'=CNA_races))
}

read_ASCAT = function(spn_id,sample_id,coverage,purity){
  caller = 'ascat'
  ascat_results <- get_sarek_cna_file(spn = spn_id,
                                      coverage = coverage,
                                      purity = purity,
                                      caller = caller,
                                      type = "tumour",
                                      sampleID = sample_id)
  
  CNA_ascat = read.csv(ascat_results$cnvs, sep='\t') %>%
    mutate(chr = paste0('chr',chr)) %>% dplyr::rename(from=startpos,to=endpos,major=nMajor,minor=nMinor) %>% as_tibble()
  
  purity_ploidy = read.csv(ascat_results$purityploidy, sep='\t')
  return(list('CNA'=CNA_ascat, 'purity_ploidy'= purity_ploidy))
}

read_Sequenza = function(spn_id,sample_id,coverage,purity){
  caller = 'sequenza'
  sequenza_results <- get_sarek_cna_file(spn = spn_id,
                                         coverage = coverage,
                                         purity = purity,
                                         caller = caller,
                                         type = "tumour",
                                         sampleID = sample_id)
  
  CNA_sequenza = read.csv(sequenza_results$segments, sep='\t')
  CNA_sequenza = CNA_sequenza %>% select(!c(N.BAF,sd.BAF,N.ratio,sd.ratio,CNt,LPP))
  colnames(CNA_sequenza) = c('chr', 'from', 'to','BAF','DR','major','minor')
  purity_ploidy_sequenza = read.csv(sequenza_results$confints_CP, sep='\t')
  return(list('CNA'=CNA_sequenza, 'purity_ploidy'= purity_ploidy_sequenza))
}

read_CNVkit = function(spn_id,sample_id,coverage,purity){
  caller='cnvkit'
  cnvkit_results <- get_sarek_cna_file(spn = spn_id,
                                       coverage = coverage,
                                       purity = purity,
                                       caller = caller,
                                       type = "tumour",
                                       sampleID = sample_id)
  
  CNA_cnvkit = read.csv(cnvkit_results$somatic.call,sep='\t')
  CNA_cnvkit = CNA_cnvkit %>% select(!c(gene,ci_hi,ci_lo,depth,probes,weight))
  colnames(CNA_cnvkit) = c('chr', 'from','to','logR','CN')

  return(list('CNA'=CNA_cnvkit))
}

## Create joint segmentation
create_joint_segmentation = function(CNA_ProCESS, CNA_target, caller, chromosomes){
  
  joint_segmentation = lapply(chromosomes, function(c){
    
    CNA_ProCESS_chr = CNA_ProCESS %>% 
      filter(chr==c) %>% 
      mutate(group = cumsum(
        lag(major, default = dplyr::first(major)) != major |
          lag(minor, default = dplyr::first(minor)) != minor |
          lag(ratio, default = dplyr::first(ratio)) != ratio
      )) %>%
      group_by(chr, major, minor, ratio, group, from ,to) %>%
      select(chr, from, to, major, minor, ratio) %>% 
      arrange(from, to)
    
    CNA_target_chr = CNA_target %>% filter(chr==c)
    froms = c(CNA_ProCESS_chr$from %>% unique(), CNA_target_chr$from %>% unique())
    tos = c(CNA_ProCESS_chr$to %>% unique(), CNA_target_chr$to %>% unique())
    breakpoints = c(froms, tos) %>% sort() %>% unique()
    df = data.frame()
    
    for (i in 1:(length(breakpoints)-1)){
      f = breakpoints[i]
      t = breakpoints[i+1]
      true_cna= CNA_ProCESS_chr %>% filter(from <= f, from < t, to>=t) 
      inferred_cna= CNA_target_chr %>% filter(from <= f, from < t, to>=t) 
      if (nrow(true_cna)==0){
        true_M1=NA
        true_m1=NA
        true_M2=NA
        true_m2=NA
        true_cn=NA
      }
      if (nrow(true_cna)==1){
        true_M1=true_cna %>% pull(major)
        true_m1=true_cna %>% pull(minor)
        true_M2=NA
        true_m2=NA
        true_cn=true_M1+true_m1
      }
      if (nrow(true_cna)>1){
        true_cna = true_cna %>% arrange(desc(ratio))
        true_M1=true_cna[1,] %>% pull(major)
        true_m1=true_cna[1,] %>% pull(minor)
        true_M2=true_cna[2,] %>% pull(major)
        true_m2=true_cna[2,] %>% pull(minor)
        ratio=true_cna[1,] %>% pull(ratio)
        true_cn = ratio*(true_M1+true_m1) + (1-ratio)*(true_M2+true_m2) 
      }
      if (nrow(inferred_cna)==0){
        inferred_M=NA
        inferred_m=NA
        inferred_cn=NA
      }
      if (nrow(inferred_cna)==1){
        if (caller!='cnvkit'){
          inferred_M=inferred_cna %>% pull(major)
          inferred_m=inferred_cna %>% pull(minor)
          inferred_cn=inferred_M+inferred_m
        }else{
          inferred_M=NA
          inferred_m=NA
          inferred_cn=inferred_cna %>% pull(CN)
        }
      }
      df = rbind(df, data.frame(
        'from' =f,
        'to'=t,
        'chromosome'=c,
        'TRUE_Major1'= true_M1,
        'TRUE_minor1'= true_m1,
        'TRUE_Major2'= true_M2,
        'TRUE_minor2'= true_m2,
        'TRUE_CN' = true_cn,
        'INFERRED_Major1' = inferred_M,
        'INFERRED_minor1' = inferred_m,
        'INFERRED_CN' = inferred_cn,
        'caller'=caller
      ))
    }
    df
    
  }) #%>% Reduce(rbind)
  joint_segmentation = Reduce(rbind, joint_segmentation)
  
  joint_segmentation_shifted = lapply(1:nrow(joint_segmentation), function(r){
    chromosome = joint_segmentation[r,]$chromosome
    from_chromosome = CNAqc:::get_reference('hg38') %>% filter(chr==chromosome) %>% pull(from)
    joint_segmentation[r,] %>% mutate(from = from + from_chromosome, to = to + from_chromosome)
  })
  joint_segmentation_shifted = Reduce(rbind, joint_segmentation_shifted)
  
  if (caller %in% c('ascat','sequenza')){
    joint_segmentation_shifted_longer = joint_segmentation_shifted %>% 
      mutate(is_match = 
               case_when(
                 TRUE_Major1==INFERRED_Major1 & TRUE_minor1==INFERRED_minor1 & is.na(TRUE_Major2) ~ 'match',
                 TRUE_Major1==INFERRED_Major1 & TRUE_minor1==INFERRED_minor1 | 
                 TRUE_Major1==INFERRED_Major1 & TRUE_minor1==INFERRED_minor1 & !is.na(TRUE_Major2) ~ 'match subclone',
                 TRUE_Major1!=INFERRED_Major1 & TRUE_minor1!=INFERRED_minor1 & !is.na(TRUE_Major2) ~ 'no match subclone',
                 .default = 'no match'
               )) %>%
      mutate(type = ifelse(!is.na(TRUE_Major2) & !is.na(TRUE_Major2), 'subclonal', 'clonal')) %>% 
      pivot_longer(
        cols = c('TRUE_Major1', 'TRUE_minor1', 'TRUE_Major2', 'TRUE_minor2', 'INFERRED_Major1', 'INFERRED_minor1'),  
        names_to = "Type",  # New column for variable names
        values_to = "Value"  # New column for values
      )
    joint_segmentation_shifted = joint_segmentation_shifted %>% filter(!is.na(TRUE_CN))
    joint_segmentation_shifted_longer = joint_segmentation_shifted_longer %>% filter(!is.na(TRUE_CN))
  } else{
    joint_segmentation_shifted = joint_segmentation_shifted %>% filter(!is.na(TRUE_CN))
    joint_segmentation_shifted_longer = joint_segmentation_shifted %>% 
      mutate(is_match = 
               case_when(
                 abs(TRUE_CN - INFERRED_CN)<.1 ~ 'match',
                 .default = 'no match'
               )) %>% 
      mutate(type = ifelse(!is.na(TRUE_Major2) & !is.na(TRUE_Major2), 'subclonal', 'clonal')) %>% 
      select(from, to, chromosome, TRUE_CN, INFERRED_CN,is_match, type) %>%
      mutate(is_match = ifelse(type == 'subclonal', 'subclonal', is_match)) %>% 
      pivot_longer(
        cols = c('TRUE_CN', 'INFERRED_CN'),  
        names_to = "Type",  # New column for variable names
        values_to = "Value"  # New column for values
      )
  }
  
  return(list('joint_segmentation'=joint_segmentation_shifted,
              'joint_segmentation_long'=joint_segmentation_shifted_longer))
}

## Compute correctness 
compute_correctness = function(df, caller) {
  if (caller %in% c('ascat', 'sequenza')){
    clonal = 1-( (df %>% filter(is_match == 'no match') %>% filter(type == 'clonal') %>% mutate(len=to-from) %>%
           pull(len) %>% unique() %>% sum()) / (df %>% filter(type == 'clonal') %>% mutate(len=to-from) %>%
                                                  pull(len) %>% unique() %>% sum()))
    
    all = 1-( (df %>% filter(is_match %in% c('no match', 'no match subclone')) %>% mutate(len=to-from) %>%
                    pull(len) %>% unique() %>% sum()) / (df %>% mutate(len=to-from) %>%
                                                           pull(len) %>% unique() %>% sum()))
  } else{
    clonal = 1-( (df %>% filter(is_match == 'no match') %>% filter(type == 'clonal') %>% mutate(len=to-from) %>%
                    pull(len) %>% unique() %>% sum()) / (df %>% filter(type == 'clonal') %>% mutate(len=to-from) %>%
                                                           pull(len) %>% unique() %>% sum()))
    all <- NA
  }
  
  return(list('all'=all, 'clonal'=clonal))
}


## True Ploidy
compute_true_ploidy = function(CNA_ProCESS){
  # Compute real ploidy
  CNA_ProCESS = CNA_ProCESS %>% mutate(len=to-from)
  genome_len =CNA_ProCESS$len %>% unique() %>% sum()
  CNA_ProCESS = CNA_ProCESS %>% mutate(seg_id = paste0(chr, ':', from,':',to))
  segments = CNA_ProCESS$seg_id %>% unique()
  ploidy = 0
  # purity_number = as.double(purity)
  for (s in segments){
    seg = CNA_ProCESS %>% filter(seg_id == s)
    nA1 = seg[1,]$major
    nB1 = seg[1,]$minor
    genome_fraction = seg[1,]$len/genome_len
    if (nrow(seg)>1){
      nA2 = seg[2,]$major
      nB2 = seg[2,]$minor
      ccf = seg[1,]$ratio
    }else{
      nA2 = 0
      nB2 = 0
      ccf = 1
    }
    
    # ploidy = ploidy + ((nA1+nB1)*ccf + (nA2+nB2)*(1-ccf))*purity_number*genome_fraction
    ploidy = ploidy + ((nA1+nB1)*ccf + (nA2+nB2)*(1-ccf))*genome_fraction
  }
  # ploidy = ploidy + 2*(1-purity_number)
  ploidy
}

## Compute segmentation correctness
covered_genome = function(CNA_target, chromosome){
  chr_len = CNAqc:::get_reference('hg38') %>% filter(chr==chromosome) %>% pull(length)
  covered = CNA_target %>% filter(chr == chromosome) %>% mutate(l=to-from) %>% pull(l) %>% sum()
  (covered / chr_len)*100
}

remove_centromeric_signal_gr <- function(signal_df, centromere_df, callable_regions_df = NULL) {
  # 1. Convert signal to GRanges with metadata
  signal_gr <- GRanges(
    seqnames = signal_df$chr,
    ranges = IRanges(start = signal_df$from, end = signal_df$to),
    major = signal_df$major,
    minor = signal_df$minor,
    ratio = signal_df$ratio,
    type = signal_df$type,
  )
  
  # 2. Convert centromeres to GRanges
  centromere_gr <- GRanges(
    seqnames = centromere_df$chr,
    ranges = IRanges(start = centromere_df$centromerStart,
                     end = centromere_df$centromerEnd)
  )
  
  # 3. Subtract centromere regions from signal
  split_gr <- GenomicRanges::setdiff(signal_gr, centromere_gr)
  
  # 4. (Optional) Intersect with callable regions
  if (!is.null(callable_regions_df)) {
    callable_gr <- GRanges(
      seqnames = callable_regions_df$chr,
      ranges = IRanges(start = callable_regions_df$start,
                       end = callable_regions_df$end)
    )
    split_gr <- GenomicRanges::intersect(split_gr, callable_gr, ignore.strand = TRUE)
  }
  
  # 5. Restore metadata if any segments remain
  if (length(split_gr) != 0) {
    hits <- findOverlaps(split_gr, signal_gr, type = "within")
    
    if (length(hits) != 0) {
      split_gr$major <- signal_gr$major[subjectHits(hits)]
      split_gr$minor <- signal_gr$minor[subjectHits(hits)]
      split_gr$ratio <- signal_gr$ratio[subjectHits(hits)]
      split_gr$type  <- signal_gr$type[subjectHits(hits)]
      
      result_df <- as.data.frame(split_gr) %>%
        transmute(
          chr = as.character(seqnames),
          from = start,
          to = end,
          major,
          minor,
          ratio,
          type,
          seg_id = paste0(chr, ":", from, ":", to)
        )
      
      return(result_df)  
    }
    
    return(signal_df)
  } else {
    return(signal_df)
  }
}

convert_centromeres_to_relative <- function(centromere_df) {
  centromere_df %>%
    mutate(
      centromerStart = centromerStart - from + 1,
      centromerEnd = centromerEnd - from + 1
    ) %>%
    select(chr, centromerStart, centromerEnd)
}

filter_genomic_signal <- function(signal_df, centromere_df, callable_regions_df = NULL) {
  
  # Convert signal dataframe to GRanges
  signal_gr <- GRanges(
    seqnames = signal_df$chr,
    ranges = IRanges(start = signal_df$from, end = signal_df$to),
    major = signal_df$major,
    ratio = signal_df$ratio,
    type = signal_df$type,
    minor = signal_df$minor
  )
  
  # Convert centromere dataframe to GRanges
  centromere_gr <- GRanges(
    seqnames = centromere_df$chr,
    ranges = IRanges(start = centromere_df$centromerStart, end = centromere_df$centromerEnd)
  )
  
  # Step 1: Remove centromeric regions from signal
  # Find overlaps between signal and centromeres
  overlaps <- findOverlaps(signal_gr, centromere_gr)
  
  if (length(overlaps) > 0) {
    # Split signal regions that overlap with centromeres
    signal_minus_centromere <- setdiff(signal_gr, centromere_gr)
    
    # For segments that were split, we need to preserve the metadata
    # Find which original segments were affected
    affected_indices <- unique(queryHits(overlaps))
    
    # Create a list to store all final segments
    final_segments <- list()
    
    # Add unaffected segments
    unaffected_segments <- signal_gr[-affected_indices]
    if (length(unaffected_segments) > 0) {
      final_segments <- append(final_segments, list(unaffected_segments))
    }
    
    # Process affected segments
    for (i in affected_indices) {
      original_segment <- signal_gr[i]
      
      # Find the parts of this segment that don't overlap with centromeres
      non_centromeric_parts <- setdiff(original_segment, centromere_gr)
      
      if (length(non_centromeric_parts) > 0) {
        # Preserve original metadata for all split parts
        mcols(non_centromeric_parts) <- mcols(original_segment)
        final_segments <- append(final_segments, list(non_centromeric_parts))
      }
    }
    
    # Combine all segments
    if (length(final_segments) > 0) {
      # Flatten the list and combine GRanges objects
      all_ranges <- unlist(GRangesList(final_segments))
      filtered_signal <- all_ranges
    } else {
      filtered_signal <- GRanges()
    }
  } else {
    filtered_signal <- signal_gr
  }
  
  # Step 2: If callable_regions_df is provided, intersect with callable regions
  if (!is.null(callable_regions_df)) {
    callable_gr <- GRanges(
      seqnames = callable_regions_df$chr,
      ranges = IRanges(start = callable_regions_df$start, end = callable_regions_df$end)
    )
    
    # Find overlaps between filtered signal and callable regions
    overlaps_callable <- findOverlaps(filtered_signal, callable_gr)
    
    if (length(overlaps_callable) > 0) {
      # Get the overlapping parts while preserving metadata
      callable_parts <- list()
      
      for (i in unique(queryHits(overlaps_callable))) {
        original_segment <- filtered_signal[i]
        
        # Find all callable regions that overlap with this segment
        overlapping_callable_indices <- subjectHits(overlaps_callable)[queryHits(overlaps_callable) == i]
        overlapping_callable <- callable_gr[overlapping_callable_indices]
        
        # Get intersection of this segment with callable regions
        intersected_parts <- intersect(original_segment, overlapping_callable)
        
        if (length(intersected_parts) > 0) {
          # Preserve metadata from original segment
          mcols(intersected_parts) <- mcols(original_segment)
          callable_parts <- append(callable_parts, list(intersected_parts))
        }
      }
      
      # Combine all callable parts
      if (length(callable_parts) > 0) {
        filtered_signal <- unlist(GRangesList(callable_parts))
      } else {
        filtered_signal <- GRanges()
      }
    } else {
      filtered_signal <- GRanges()
    }
  }
  
  # Convert back to dataframe
  if (length(filtered_signal) > 0) {
    result_df <- data.frame(
      chr = as.character(seqnames(filtered_signal)),
      from = start(filtered_signal),
      to = end(filtered_signal),
      major = mcols(filtered_signal)$major,
      ratio = mcols(filtered_signal)$ratio,
      type = mcols(filtered_signal)$type,
      minor = mcols(filtered_signal)$minor,
      stringsAsFactors = FALSE
    )
  } else {
    # Return empty dataframe with correct structure
    result_df <- data.frame(
      chr = character(0),
      from = numeric(0),
      to = numeric(0),
      major = numeric(0),
      ratio = numeric(0),
      type = character(0),
      minor = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  return(result_df)
}
